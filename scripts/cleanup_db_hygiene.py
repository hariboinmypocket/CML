#!/usr/bin/env python3
"""Five independent data-hygiene fixes, each reported separately.

1. studies.data_type carried five values for three real platforms: 16S (181)
   alongside "16S rRNA" (29), and mNGS (93) alongside "shotgun metagenomics"
   (9). The species-vs-genus decision and `sync-abundances --data-type` both
   key on this column, so 38 studies silently fell out of any filter.

2. Missing values stored as real numbers. PRJEB6172 records absent BMI as 0.0
   (52 rows) while its real BMIs run 16.3-31.8, and absent age as 0 in the
   same rows. Two samples carry age -88 in a cohort otherwise aged 3-28, and
   one carries 119 where the next oldest is 71. A model reading these as
   measurements learns from a sentinel.

   Age 0 is only nulled in cohorts where it cannot be a newborn: PRJEB6172,
   PRJEB41297 and PRJNA1125836 have no sample under 16, 19 and 23 years
   respectively. PRJNA716780 is left alone -- it holds 47 samples under two
   years old and a minimum non-zero age of 0.1, so its single 0 is plausibly
   a real birth-day sample, not a sentinel.

3. Ten bacterial genera carry a plant or animal lineage from a cross-kingdom
   homonym, with the wrong ncbi_tax_id to match: Gordonia resolved to tax ID
   79255, a tea plant, across 159 abundance rows; Planococcus to a mealybug;
   Coxiella and Microcystis to molluscs. Re-resolved through
   TaxonomyClient.resolve_name with the domain constrained to Bacteria, which
   is the guard added for exactly this failure. Where NCBI cannot confirm a
   bacterial match the known-wrong lineage is cleared rather than kept, since
   a confidently wrong phylum is worse than a missing one.

4. 53 rows sit under phylum Thermodesulfobacteria because the synonym table
   folded NCBI's modern Thermodesulfobacteriota onto the classic thermophile
   phylum. Re-resolved from each row's tax ID now that the mapping is fixed.

5. sample_taxon_abundances.detection_threshold is NULL in all 2.58M rows.
   load_abundances did read a "detection_threshold" input column, but no
   source this pipeline ingests has ever supplied one, so the column only
   ever stored NULL. That read is removed alongside the column.

6. Ten species rows sit under phylum Arthropoda, which none of them are. Nine
   are MiMeDB "Bacillus <epithet>" rows that NCBI renamed in the 2020 Bacillus
   split (B. firmus -> Cytobacillus firmus, B. muralis -> Peribacillus
   muralis); the tenth is Victoria amazonica, a water lily whose tax ID is
   right and whose stored phylum is not.

Dry run by default; pass --apply to write. Pass --only=STEP (1-6, repeatable)
to run a subset.
"""
from __future__ import annotations

import csv
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.ncbi import TaxonomyClient
from gutdb.transform import normalize_phylum, parse_scientific_name, same_epithet, taxon_key

DRY_RUN = "--apply" not in sys.argv
EMAIL = os.environ.get("NCBI_EMAIL") or next(
    (a.split("=", 1)[1] for a in sys.argv if a.startswith("--email=")), None
)
ONLY = {a.split("=", 1)[1] for a in sys.argv if a.startswith("--only=")}

DATA_TYPE_MERGE = {"16S rRNA": "16S", "shotgun metagenomics": "mNGS"}

# Cohorts whose youngest real sample rules out a newborn, so age 0 is a
# sentinel. PRJNA716780 is deliberately absent: see the module docstring.
AGE_ZERO_SENTINEL_STUDIES = ("PRJEB6172", "PRJEB41297", "PRJNA1125836")

HOMONYM_GENERA = (
    "Gordonia", "Planococcus", "Microcystis", "Coxiella", "Yersinia",
    "Syntrophus", "Buchnera", "Schwartzia", "Moorella", "Lawsonia",
)
LINEAGE_FIELDS = ("superkingdom", "phylum", "class_name", "order_name", "family")
NCBI_RANK_FOR_FIELD = {
    "superkingdom": ("domain", "superkingdom"),
    "phylum": ("phylum",),
    "class_name": ("class",),
    "order_name": ("order",),
    "family": ("family",),
}

audit: list[dict] = []


def record(step: str, target: str, before: str, after: str, note: str = "") -> None:
    audit.append({"step": step, "target": target, "before": before,
                  "after": after, "note": note})


def wanted(step: int) -> bool:
    return not ONLY or str(step) in ONLY


def step1_data_type(cur) -> dict:
    stats = {"studies_updated": 0}
    for old, new in DATA_TYPE_MERGE.items():
        cur.execute("SELECT COUNT(*) FROM studies WHERE data_type = %s", (old,))
        n = cur.fetchone()[0]
        if not n:
            continue
        record("1_data_type", f"{n} studies", old, new, "same platform, two spellings")
        if not DRY_RUN:
            cur.execute("UPDATE studies SET data_type = %s WHERE data_type = %s", (new, old))
            stats["studies_updated"] += cur.rowcount
        else:
            stats["studies_updated"] += n
    return stats


def step2_sentinels(cur) -> dict:
    stats = {"bmi_nulled": 0, "age_nulled": 0}

    cur.execute("SELECT COUNT(*) FROM samples WHERE bmi = 0")
    n = cur.fetchone()[0]
    if n:
        record("2_sentinel", f"{n} samples", "bmi = 0", "NULL",
               "PRJEB6172 records absent BMI as 0; its real range is 16.3-31.8")
        if not DRY_RUN:
            cur.execute("UPDATE samples SET bmi = NULL WHERE bmi = 0")
            stats["bmi_nulled"] = cur.rowcount
        else:
            stats["bmi_nulled"] = n

    placeholders = ",".join(["%s"] * len(AGE_ZERO_SENTINEL_STUDIES))
    cur.execute(
        f"""SELECT COUNT(*) FROM samples s JOIN studies st ON st.id = s.study_id
            WHERE s.age_years = 0 AND st.project_accession IN ({placeholders})""",
        AGE_ZERO_SENTINEL_STUDIES,
    )
    n = cur.fetchone()[0]
    if n:
        record("2_sentinel", f"{n} samples", "age_years = 0", "NULL",
               "adult-only cohorts; PRJNA716780 excluded as a real infant cohort")
        if not DRY_RUN:
            cur.execute(
                f"""UPDATE samples s JOIN studies st ON st.id = s.study_id
                    SET s.age_years = NULL
                    WHERE s.age_years = 0 AND st.project_accession IN ({placeholders})""",
                AGE_ZERO_SENTINEL_STUDIES,
            )
            stats["age_nulled"] += cur.rowcount
        else:
            stats["age_nulled"] += n

    cur.execute("SELECT COUNT(*) FROM samples WHERE age_years < 0 OR age_years > 110")
    n = cur.fetchone()[0]
    if n:
        record("2_sentinel", f"{n} samples", "age_years out of [0,110]", "NULL",
               "-88 in a cohort aged 3-28; 119 where the next oldest is 71")
        if not DRY_RUN:
            cur.execute(
                "UPDATE samples SET age_years = NULL WHERE age_years < 0 OR age_years > 110")
            stats["age_nulled"] += cur.rowcount
        else:
            stats["age_nulled"] += n
    return stats


def step3_homonyms(cur, client) -> dict:
    stats = {"resolved": 0, "cleared": 0}
    placeholders = ",".join(["%s"] * len(HOMONYM_GENERA))
    cur.execute(
        f"""SELECT id, genus, scientific_name, ncbi_tax_id, phylum FROM taxa
            WHERE genus IN ({placeholders})
              AND phylum IN ('Streptophyta','Arthropoda','Mollusca','Chordata')""",
        HOMONYM_GENERA,
    )
    rows = cur.fetchall()
    for taxon_id, genus, sci, old_taxid, old_phylum in rows:
        taxid = client.resolve_name(genus, "Bacteria")
        time.sleep(0.4)
        if not taxid:
            record("3_homonym", f"taxon {taxon_id} {sci}",
                   f"taxid {old_taxid}, phylum {old_phylum}", "lineage cleared",
                   "NCBI could not confirm a bacterial match")
            if not DRY_RUN:
                cur.execute(
                    """UPDATE taxa SET ncbi_tax_id = NULL, superkingdom = NULL,
                       phylum = NULL, class_name = NULL, order_name = NULL, family = NULL
                       WHERE id = %s""",
                    (taxon_id,),
                )
                stats["cleared"] += 1
            else:
                stats["cleared"] += 1
            continue

        lineage = client.fetch_lineages([taxid]).get(taxid, {})
        time.sleep(0.4)
        values = {}
        for field, ranks in NCBI_RANK_FOR_FIELD.items():
            for rank in ranks:
                if lineage.get(rank):
                    values[field] = lineage[rank]
                    break
        if values.get("phylum"):
            values["phylum"] = normalize_phylum(values["phylum"])

        record("3_homonym", f"taxon {taxon_id} {sci}",
               f"taxid {old_taxid}, phylum {old_phylum}",
               f"taxid {taxid}, phylum {values.get('phylum','?')}",
               "re-resolved with the search constrained to Bacteria")
        if not DRY_RUN:
            assignments = ", ".join(f"`{f}` = %s" for f in values)
            cur.execute(
                f"UPDATE taxa SET ncbi_tax_id = %s{', ' + assignments if values else ''} WHERE id = %s",
                (int(taxid), *values.values(), taxon_id),
            )
        stats["resolved"] += 1
    return stats


def step4_phylum(cur, client) -> dict:
    stats = {"rephylumed": 0, "unchanged": 0, "no_taxid": 0}
    cur.execute(
        "SELECT id, scientific_name, ncbi_tax_id FROM taxa WHERE phylum = 'Thermodesulfobacteria'")
    rows = cur.fetchall()
    with_taxid = [(i, s, str(t)) for i, s, t in rows if t]
    stats["no_taxid"] = len(rows) - len(with_taxid)

    lineages: dict[str, dict] = {}
    ids = [t for _, _, t in with_taxid]
    for i in range(0, len(ids), 50):
        lineages.update(client.fetch_lineages(ids[i : i + 50]))
        time.sleep(0.4)

    for taxon_id, sci, taxid in with_taxid:
        raw = lineages.get(taxid, {}).get("phylum", "")
        new = normalize_phylum(raw) if raw else ""
        if not new or new == "Thermodesulfobacteria":
            stats["unchanged"] += 1
            continue
        record("4_phylum", f"taxon {taxon_id} {sci}", "Thermodesulfobacteria", new,
               "re-resolved from NCBI after the synonym fix")
        if not DRY_RUN:
            cur.execute("UPDATE taxa SET phylum = %s WHERE id = %s", (new, taxon_id))
        stats["rephylumed"] += 1
    return stats


def step5_drop_column(cur) -> dict:
    cur.execute(
        """SELECT COUNT(*) FROM information_schema.COLUMNS
           WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'sample_taxon_abundances'
             AND COLUMN_NAME = 'detection_threshold'""")
    if not cur.fetchone()[0]:
        return {"already_absent": 1}
    cur.execute(
        "SELECT COUNT(*) FROM sample_taxon_abundances WHERE detection_threshold IS NOT NULL")
    populated = cur.fetchone()[0]
    if populated:
        # Refuse rather than destroy data the survey said was absent.
        print(f"  !! detection_threshold has {populated} populated rows - NOT dropping")
        return {"refused_populated_rows": populated}
    record("5_drop_column", "sample_taxon_abundances.detection_threshold",
           "column present, 0 of 2.58M rows populated", "dropped",
           "no ingested source has ever supplied the column")
    if not DRY_RUN:
        cur.execute("ALTER TABLE sample_taxon_abundances DROP COLUMN detection_threshold")
    return {"column_dropped": 1}


def step6_arthropoda_species(cur, client) -> dict:
    """Species rows filed under phylum Arthropoda, which none of them are.

    Two unrelated causes sit in the same bucket:

    - Nine MiMeDB rows named "Bacillus <epithet>" with no tax ID. These are
      real organisms that NCBI renamed in the 2020 Bacillus split, so
      resolve_name refuses them: it will not claim "Bacillus firmus" is tax ID
      1399 when NCBI calls that Cytobacillus firmus. Identity is established
      here by the species epithet surviving the genus change, which is the
      same test the duplicate-tax-ID merge uses, and only then is the record
      adopted -- tax ID, lineage and NCBI's current name.
    - Victoria amazonica, a water lily, whose tax ID (85961) is correct and
      whose stored phylum is simply wrong. Its lineage is rewritten from its
      own tax ID; it stays a plant, since the MiMeDB reference library holds
      plants deliberately.

    Where the epithet does not survive, the genus is resolved for lineage only
    and ncbi_tax_id is left NULL. Writing a genus tax ID onto a species row is
    the defect that put the F. prausnitzii strain ID on the species.
    """
    stats = {"renamed_from_synonym": 0, "lineage_only": 0, "retaxid_lineage": 0,
             "deleted_stale_duplicate": 0, "needs_merge": 0, "unresolved": 0}
    cur.execute(
        "SELECT id, genus, species, scientific_name, ncbi_tax_id FROM taxa "
        "WHERE phylum = 'Arthropoda' ORDER BY id")
    rows = cur.fetchall()

    def lineage_values(record: dict) -> dict:
        values = {}
        for field, ranks in NCBI_RANK_FOR_FIELD.items():
            for rank in ranks:
                if record.get(rank):
                    values[field] = record[rank]
                    break
        if values.get("phylum"):
            values["phylum"] = normalize_phylum(values["phylum"])
        return values

    def write(taxon_id: int, values: dict, taxid: str | None) -> None:
        if DRY_RUN or not values:
            return
        assignments = ", ".join(f"`{f}` = %s" for f in values)
        params = list(values.values())
        if taxid is not None:
            assignments += ", ncbi_tax_id = %s"
            params.append(int(taxid))
        cur.execute(f"UPDATE taxa SET {assignments} WHERE id = %s", (*params, taxon_id))

    for taxon_id, genus, species, sci, existing_taxid in rows:
        if existing_taxid:
            record = client.fetch_lineages([str(existing_taxid)]).get(str(existing_taxid), {})
            time.sleep(0.35)
            values = lineage_values(record)
            if not values:
                stats["unresolved"] += 1
                continue
            record_arthropoda(taxon_id, sci, f"phylum Arthropoda, taxid {existing_taxid}",
                              f"phylum {values.get('phylum','?')}",
                              "tax ID was already correct; only the stored lineage was wrong")
            write(taxon_id, values, None)
            stats["retaxid_lineage"] += 1
            continue

        taxid, ncbi_name, record = client.resolve_synonym(sci, "Bacteria")
        time.sleep(0.35)
        _, ncbi_epithet = parse_scientific_name(ncbi_name) if ncbi_name else ("", "")
        if taxid and same_epithet(species, ncbi_epithet):
            values = lineage_values(record)
            genus_key, species_key = taxon_key(*parse_scientific_name(ncbi_name))
            cur.execute(
                "SELECT id FROM taxa WHERE genus_key = %s AND species_key = %s AND id <> %s",
                (genus_key, species_key, taxon_id))
            clash = cur.fetchone()
            if clash:
                # A row already holds NCBI's current name for this organism, so
                # this one is a stale duplicate under the superseded name.
                # Adopting the tax ID here would recreate exactly the duplicate
                # pair the merge removes, so an empty row is deleted outright
                # and a row carrying data is left for the merge to fold.
                cur.execute(
                    """SELECT (SELECT COUNT(*) FROM taxon_disease_associations a
                                WHERE a.taxon_id = %s),
                              (SELECT COUNT(*) FROM sample_taxon_abundances b
                                WHERE b.taxon_id = %s)""",
                    (taxon_id, taxon_id))
                n_assoc, n_abund = cur.fetchone()
                if n_assoc or n_abund:
                    record_arthropoda(
                        taxon_id, sci, "phylum Arthropoda, no taxid",
                        f"taxid {taxid}, phylum {values.get('phylum','?')}",
                        f"carries {n_assoc} assoc / {n_abund} abund; {ncbi_name!r} is "
                        f"row {clash[0]} - run merge_duplicate_taxa.py to fold them")
                    write(taxon_id, values, taxid)
                    stats["needs_merge"] += 1
                else:
                    record_arthropoda(
                        taxon_id, sci, "phylum Arthropoda, no taxid, no data", "deleted",
                        f"stale duplicate of row {clash[0]} {ncbi_name!r}, which already "
                        f"holds tax ID {taxid} and the correct lineage")
                    if not DRY_RUN:
                        cur.execute("DELETE FROM taxa WHERE id = %s", (taxon_id,))
                    stats["deleted_stale_duplicate"] += 1
                continue
            record_arthropoda(taxon_id, sci, "phylum Arthropoda, no taxid",
                              f"{ncbi_name}, taxid {taxid}, phylum {values.get('phylum','?')}",
                              "renamed: epithet survives the genus change, so same organism")
            if not DRY_RUN:
                new_genus, new_species = parse_scientific_name(ncbi_name)
                values_with_name = dict(values)
                cur.execute(
                    f"""UPDATE taxa SET genus = %s, species = %s, genus_key = %s,
                        species_key = %s, scientific_name = %s, ncbi_tax_id = %s,
                        {', '.join(f'`{f}` = %s' for f in values_with_name)}
                        WHERE id = %s""",
                    (new_genus, new_species, genus_key, species_key, ncbi_name,
                     int(taxid), *values_with_name.values(), taxon_id))
            stats["renamed_from_synonym"] += 1
            continue

        # Fall back to the genus for lineage only, never its tax ID.
        genus_taxid, genus_name, genus_record = client.resolve_synonym(genus, "Bacteria")
        time.sleep(0.35)
        if not genus_taxid or genus_name.casefold() != genus.casefold():
            record_arthropoda(taxon_id, sci, "phylum Arthropoda", "lineage cleared",
                              "neither the species nor the genus could be confirmed")
            if not DRY_RUN:
                cur.execute(
                    """UPDATE taxa SET superkingdom = NULL, phylum = NULL, class_name = NULL,
                       order_name = NULL, family = NULL WHERE id = %s""", (taxon_id,))
            stats["unresolved"] += 1
            continue
        values = lineage_values(genus_record)
        values.pop("family", None)  # a genus record's family is right; its species is unknown
        record_arthropoda(taxon_id, sci, "phylum Arthropoda, no taxid",
                          f"phylum {values.get('phylum','?')} (genus lineage, tax ID left NULL)",
                          "species epithet did not survive; genus lineage only")
        write(taxon_id, values, None)
        stats["lineage_only"] += 1
    return stats


def record_arthropoda(taxon_id: int, sci: str, before: str, after: str, note: str) -> None:
    record("6_arthropoda", f"taxon {taxon_id} {sci}", before, after, note)


def main() -> int:
    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    client = TaxonomyClient(email=EMAIL)
    results: dict[str, dict] = {}

    if wanted(1):
        results["1_data_type"] = step1_data_type(cur)
    if wanted(2):
        results["2_sentinels"] = step2_sentinels(cur)
    if wanted(3):
        results["3_homonyms"] = step3_homonyms(cur, client)
    if wanted(4):
        results["4_phylum"] = step4_phylum(cur, client)
    if wanted(5):
        results["5_drop_column"] = step5_drop_column(cur)
    if wanted(6):
        results["6_arthropoda"] = step6_arthropoda_species(cur, client)

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    else:
        conn.commit()

    for step, stats in results.items():
        print(f"\n{step}")
        for key, value in stats.items():
            print(f"    {key:28} {value}")

    # The per-change detail is printed above and nowhere else. It used to be
    # written to data/db_hygiene_audit.csv as well, which was removed: the file
    # recorded the sample-level fixes only in aggregate ("52 samples", "41
    # samples") and never which samples, so it could not be used to reverse or
    # verify anything the commit message did not already state. Worse, it was
    # rewritten on every run, so a `--only=2` pass clobbered the record of the
    # other four steps and left a file that read as though those had never
    # happened.
    if audit:
        print(f"\n{len(audit)} changes, detailed above")

    cur.close()
    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
