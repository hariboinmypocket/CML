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

Dry run by default; pass --apply to write. Pass --only STEP (1-5, repeatable)
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
from gutdb.transform import normalize_phylum

DRY_RUN = "--apply" not in sys.argv
AUDIT = Path("data/db_hygiene_audit.csv")
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

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    else:
        conn.commit()

    for step, stats in results.items():
        print(f"\n{step}")
        for key, value in stats.items():
            print(f"    {key:28} {value}")

    if audit:
        AUDIT.parent.mkdir(parents=True, exist_ok=True)
        with AUDIT.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(audit[0].keys()))
            writer.writeheader()
            writer.writerows(audit)
        print(f"\naudit written to {AUDIT} ({len(audit)} entries)")

    cur.close()
    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
