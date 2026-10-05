#!/usr/bin/env python3
"""Load MicrobiomeHD's genus-level q-values (Duvallet et al. 2017, file-S1).

Duvallet et al. (Nat Commun 8:1784, PMID 29209090) re-processed 28 case-control
16S studies through a single standardized pipeline -- 100% de novo OTUs, RDP
classifier, collapsed to genus, Kruskal-Wallis with Benjamini-Hochberg FDR --
and published the resulting q-values per genus per dataset. That makes it an
independent, consistently-processed evidence source to set against the
per-project LEfSe results this database holds from GMrepo.

file-S1 stores a SIGNED q-value: the magnitude is the FDR-corrected q, and the
sign says which arm the genus is higher in (positive = cases, negative =
controls). Direction is taken from the sign and the magnitude goes to the new
q_value column; effect_size stays NULL because S1 reports no effect size.
file-S5's log2 fold-changes are the matching effect sizes and are not loaded
here.

Only the 434 rows significant at |q| < 0.05 are loaded, out of 2,769 present in
the matrix. The remaining 2,335 are genuine "tested, not significant" results
and are informative -- but v_taxon_specificity counts COUNT(DISTINCT disease_id)
over every association row with no direction filter, so loading them as
direction='no_difference' would push any merely-tested genus toward 'broad' or
'pan_disease' and corrupt the specificity classification. They belong in a
separate table if they are ever wanted, not here.

Identity: the 30 datasets come from 28 papers (cdi_schubert/noncdi_schubert
share PMID 24803517; nash_zhu/ob_zhu share PMID 23055155), so each paper is one
studies row keyed PMID<pmid> to match this database's existing convention, and
each dataset is its own phenotype_comparison on that study. Every PMID and
accession was checked against the existing 314 studies: none is present, so
this adds no duplicate evidence.

Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import csv
import os
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.pipeline import (
    LoadStats,
    _finish_run,
    _start_run,
    _upsert_comparison,
    upsert_disease,
    upsert_study,
    upsert_taxon,
)
from gutdb.transform import clean_str, normalize_phylum

DRY_RUN = "--apply" not in sys.argv
# Upstream files are not redistributed here; point MICROBIOMEHD_DIR at a copy or
# fetch them into the default location (see FETCH_HINT below).
SOURCE = Path(os.environ.get("MICROBIOMEHD_DIR", "data/microbiomehd"))
QVALUES = SOURCE / "file-S1.qvalues.txt"
EFFECTS = SOURCE / "file-S5.effects.txt"
LITERATURE = SOURCE / "file-S4.literature_results.txt"
IDENTITY = SOURCE / "dataset_identity.csv"
AUDIT = Path("data/microbiomehd_s1_audit.csv")

RAW = "https://raw.githubusercontent.com/cduvallet/microbiomeHD/master"
FETCH_HINT = f"""
Required files are missing from {SOURCE}/. To fetch them:

  mkdir -p {SOURCE}
  curl -sS -o {SOURCE}/file-S1.qvalues.txt \\
    {RAW}/final/supp-files/file-S1.qvalues.txt
  curl -sS -o {SOURCE}/file-S5.effects.txt \\
    {RAW}/final/supp-files/file-S5.effects.txt
  curl -sS -o {SOURCE}/file-S4.literature_results.txt \\
    {RAW}/final/supp-files/file-S4.literature_results.txt
  curl -sS -o {SOURCE}/results_folders.yaml \\
    {RAW}/data/user_input/results_folders.yaml

dataset_identity.csv is derived from results_folders.yaml by resolving each
dataset's paper DOI to a PMID through NCBI; see scripts/README or rebuild it
with the snippet in this script's docstring.
"""

Q_THRESHOLD = 0.05
METHOD = "Kruskal-Wallis + BH FDR"
HEALTH = ("Health", "D006262")

# MicrobiomeHD's dataset-id prefix -> (disease name, MeSH id). Every name here
# was matched against the existing diseases table; only Hepatic Encephalopathy
# is new. 'edd' is enteric diarrheal disease and 'noncdi' is the non-C.difficile
# diarrhea arm of Schubert et al., so both map to Diarrhea.
DISEASE_FOR_CODE = {
    "art": ("Arthritis, Rheumatoid", "D001172"),
    "asd": ("Autism Spectrum Disorder", "D000067877"),
    "cdi": ("Clostridium Infections", "D003015"),
    "crc": ("Colorectal Neoplasms", "D015179"),
    "edd": ("Diarrhea", "D003967"),
    "hiv": ("HIV Infections", "D015658"),
    "ibd": ("Inflammatory Bowel Diseases", "D015212"),
    "liv": ("Hepatic Encephalopathy", "D006501"),
    "nash": ("Non-alcoholic Fatty Liver Disease", "D065626"),
    "noncdi": ("Diarrhea", "D003967"),
    "ob": ("Obesity", "D009765"),
    "par": ("Parkinson Disease", "D010300"),
    "t1d": ("Diabetes Mellitus, Type 1", "D003922"),
}

RANK_PREFIX = {"k__": "superkingdom", "p__": "phylum", "c__": "class_name",
               "o__": "order_name", "f__": "family", "g__": "genus"}
# RDP cluster labels and placeholder buckets, not organisms. Loading these as
# taxa would invent genera that do not exist.
NOT_A_GENUS = re.compile(
    r"(_incertae_sedis$|^Clostridium_[IVX]|_unclassified$|^unclassified|^$)", re.I)


def parse_lineage(label: str) -> dict[str, str]:
    """Greengenes-style 'k__X;p__Y;...;g__Z' -> our lineage column names."""
    out: dict[str, str] = {}
    for part in label.split(";"):
        part = part.strip()
        for prefix, column in RANK_PREFIX.items():
            if part.startswith(prefix):
                value = part[len(prefix):].strip().strip("[]")
                if value:
                    out[column] = normalize_phylum(value) if column == "phylum" else value
    return out


def load_effects() -> tuple[dict[str, list[str]], list[str], set[float]]:
    """Return (lineage -> row, datasets, sentinel values) from file-S5.

    S5 holds log2(mean_cases / mean_controls), but three of its values are
    placeholders rather than measurements: the table maximum stands in for a
    fold-change against a zero control mean, the table minimum for a zero case
    mean, and 0.0 for both means being zero. Those are division-by-zero markers,
    so storing the max as an effect size would assert a ~1300-fold enrichment
    the data cannot support. The sentinels are derived from the file instead of
    hardcoded, so a regenerated file with different extremes still works.
    """
    if not EFFECTS.exists():
        return {}, [], set()
    rows = list(csv.reader(EFFECTS.open(), delimiter="\t"))
    datasets = rows[0][1:]
    values: list[float] = []
    table: dict[str, list[str]] = {}
    for row in rows[1:]:
        table[row[0]] = row[1:]
        for raw in row[1:]:
            raw = raw.strip()
            if not raw:
                continue
            try:
                values.append(float(raw))
            except ValueError:
                pass
    sentinels = {max(values), min(values), 0.0} if values else set()
    return table, datasets, sentinels


TAXDUMP = Path(os.environ.get("GUTDB_TAXDUMP", "taxdump"))
_prokaryote_genera: set[str] | None = None


def taxdump_has_prokaryote_genus(name: str) -> bool:
    """True when `name` is a genus NCBI places under Bacteria or Archaea.

    Used to reject names that are not organisms before they become taxa rows.
    When no taxdump is available the check cannot run, so it passes rather than
    silently dropping valid findings -- the caller's report says which happened.
    """
    global _prokaryote_genera
    if _prokaryote_genera is None:
        names, nodes = TAXDUMP / "names.dmp", TAXDUMP / "nodes.dmp"
        if not (names.exists() and nodes.exists()):
            print(f"    (no taxdump at {TAXDUMP}/ - genus names not validated)")
            _prokaryote_genera = set()
            return True
        parent: dict[int, int] = {}
        rank: dict[int, str] = {}
        with nodes.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = line.split("\t|")
                parent[int(parts[0])] = int(parts[1])
                rank[int(parts[0])] = parts[2].strip()
        # Bacteria = 2, Archaea = 2157
        domains = {2, 2157}
        def is_prokaryote(taxid: int) -> bool:
            node, guard = taxid, 0
            while node and node != 1 and guard < 60:
                if node in domains:
                    return True
                node = parent.get(node, 0)
                guard += 1
            return False
        genera = set()
        with names.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = line.split("\t|")
                if len(parts) < 4 or parts[3].strip() != "scientific name":
                    continue
                taxid = int(parts[0])
                if rank.get(taxid) == "genus" and is_prokaryote(taxid):
                    genera.add(parts[1].strip().casefold())
        _prokaryote_genera = genera
        print(f"    taxdump: {len(genera):,} prokaryote genus names loaded")
    if not _prokaryote_genera:
        return True
    return clean_str(name).casefold() in _prokaryote_genera


# Control-arm labels used across S4's hand-curated group columns.
S4_CONTROL = {"control", "h", "healthy", "hc", "nt", "non-ibd", "nonibd"}


def load_literature(s1_datasets: list[str]) -> tuple[list[dict], Counter, dict[str, str]]:
    """Return (loadable rows, attrition counts, S4 study id -> S1 dataset id).

    file-S4 is the authors' hand-curated record of what each ORIGINAL paper
    reported, which makes it a third view: the publication's own claim, beside
    MicrobiomeHD's re-analysis (S1/S5) and this database's GMrepo evidence. It is
    notes rather than a matrix, though, and only a small part of it fits a
    genus-keyed association table:

      - 640 of 1,027 rows are not genus-level (OTU, species, phylum, family).
        taxa.genus is NOT NULL, so a family- or phylum-level finding has no key.
      - 512 rows belong to studies with no S1 counterpart. Their ids are
        abbreviated differently ("ra_littman", "ibd_hut", "mhe_zhang") and the
        repository documents no mapping, so they are left alone rather than
        matched on a shared author name -- attaching a paper's claim to the
        wrong study is worse than omitting it.
      - 497 rows report no parseable q value.

    What survives is 31 rows from 5 studies. That is 3% of the file, so this is
    a sample of the literature, not a summary of it; the attrition is returned
    so a caller can report it rather than imply coverage it does not have.
    """
    if not LITERATURE.exists():
        return [], Counter(), {}
    text = LITERATURE.read_bytes().decode("cp1252")  # not UTF-8 upstream
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))

    mapping: dict[str, str] = {}
    for sid in sorted({r["study"] for r in rows}):
        if sid in s1_datasets:
            mapping[sid] = sid
            continue
        code = sid.split("_")[0]
        exact = [d for d in s1_datasets if d.startswith(sid) and d.split("_")[0] == code]
        if len(exact) == 1:
            mapping[sid] = exact[0]

    attrition = Counter()
    collected: dict[tuple[str, str, str], dict] = {}
    for row in rows:
        # Species rows load too: taxa holds (genus, species), so a binomial
        # finding has a key. Ranks above genus do not -- taxa.genus is NOT NULL,
        # so a family- or phylum-level row cannot be represented without making
        # the column nullable, which is a schema change, not a cleanup.
        level = row["taxonomic_level"].casefold()
        if level not in ("genus", "species", "species_metagenomics"):
            attrition["rank_above_genus_or_otu"] += 1
            continue
        match = re.search(r"g__([A-Za-z0-9_\-]+)", row["taxa_or_feature"])
        if not match:
            attrition["genus_field_empty"] += 1
            continue
        epithet = ""
        if level != "genus":
            species_match = re.search(r"s__([A-Za-z0-9_\-]+)", row["taxa_or_feature"])
            if not species_match:
                attrition["species_field_empty"] += 1
                continue
            # Subspecies are collapsed to the species: "s__nucleatum;sb__polymorphum"
            # and "s__nucleatum;sb__nucleatum" are two findings about one organism
            # as far as a (genus, species) key is concerned, so the stronger q is
            # kept rather than one arbitrarily overwriting the other.
            epithet = species_match.group(1).casefold()
        if row["study"] not in mapping:
            attrition["study_unmappable"] += 1
            continue
        lower = row["group_lower"].strip().casefold()
        higher = row["group_higher"].strip().casefold()
        if (lower in S4_CONTROL) == (higher in S4_CONTROL):
            attrition["not_case_vs_control"] += 1
            continue
        try:
            qval = float(row["qval"])
        except ValueError:
            attrition["no_parseable_qval"] += 1
            continue
        if qval >= Q_THRESHOLD:
            attrition["not_significant"] += 1
            continue
        entry = {
            "dataset": mapping[row["study"]],
            "s4_study": row["study"],
            "genus": match.group(1),
            "species": epithet,
            "rank": "species" if epithet else "genus",
            "qval": qval,
            "direction": "depleted" if higher in S4_CONTROL else "enriched",
            "method": row["method"].strip() or "as reported",
            "multi_comp": row["multi_comp"].strip(),
        }
        key = (entry["dataset"], entry["genus"].casefold(), epithet)
        existing = collected.get(key)
        if existing is None:
            collected[key] = entry
        elif qval < existing["qval"]:
            collected[key] = entry
            attrition["subspecies_collapsed"] += 1
        else:
            attrition["subspecies_collapsed"] += 1
    return list(collected.values()), attrition, mapping


def load_identity() -> dict[str, dict[str, str]]:
    return {r["dataset"]: r for r in csv.DictReader(IDENTITY.open())}


def main() -> int:
    for required in (QVALUES, IDENTITY):
        if not required.exists():
            raise SystemExit(FETCH_HINT)
    identity = load_identity()
    s5_table, s5_datasets, s5_sentinels = load_effects()
    s5_index = {name: i for i, name in enumerate(s5_datasets)}
    if s5_table:
        print(f"S5 effects: {len(s5_table)} genera, sentinels "
              f"{sorted(round(v, 6) for v in s5_sentinels)}")
    else:
        print("S5 effects: file absent - effect sizes will be left NULL")
    rows = list(csv.reader(QVALUES.open(), delimiter="\t"))
    datasets = rows[0][1:]
    print(f"datasets in S1: {len(datasets)}   genera: {len(rows) - 1}")

    missing = [d for d in datasets if d not in identity]
    if missing:
        print(f"  !! no identity record for: {missing}")

    conn = connect(Settings.from_env(".env"))
    run_id = None if DRY_RUN else _start_run(
        conn, "MicrobiomeHD_S1_qvalues",
        "https://github.com/cduvallet/microbiomeHD final/supp-files/file-S1.qvalues.txt")

    stats = Counter()
    audit: list[dict] = []
    study_cache: dict[str, int] = {}
    comparison_cache: dict[str, int] = {}
    disease_cache: dict[str, int] = {}
    taxon_cache: dict[str, int] = {}

    def study_for(accession: str, dataset: str) -> int:
        if accession not in study_cache:
            study_cache[accession] = upsert_study(conn, {
                "project_id": accession,
                "title": f"MicrobiomeHD standardized re-analysis ({dataset})",
                "description": (
                    "Case-control 16S study re-processed by Duvallet et al. 2017 "
                    "(PMID 29209090) through one standardized pipeline: de novo OTUs, "
                    "RDP classifier, collapsed to genus, Kruskal-Wallis with "
                    "Benjamini-Hochberg FDR."),
                "data_type": "16S",
                "source_database": "MicrobiomeHD",
                "data_quality": "curated",
            })
        return study_cache[accession]

    def disease_id(name: str, mesh: str) -> int:
        if name not in disease_cache:
            disease_cache[name] = upsert_disease(conn, name, mesh)
        return disease_cache[name]

    health_id = None if DRY_RUN else disease_id(*HEALTH)

    for column, dataset in enumerate(datasets, start=1):
        code = dataset.split("_")[0]
        if code not in DISEASE_FOR_CODE:
            print(f"  !! unmapped disease code {code!r} ({dataset}) - skipped")
            stats["datasets_skipped"] += 1
            continue
        name, mesh = DISEASE_FOR_CODE[code]
        rec = identity.get(dataset, {})
        accession = f"PMID{rec['pmid']}" if rec.get("pmid") else rec.get("accession") or dataset

        significant = []
        for row in rows[1:]:
            raw = row[column].strip() if column < len(row) else ""
            if not raw:
                continue
            try:
                q = float(raw)
            except ValueError:
                continue
            stats["qvalues_present"] += 1
            if abs(q) >= Q_THRESHOLD:
                stats["not_significant"] += 1
                continue
            lineage = parse_lineage(row[0])
            genus = lineage.get("genus", "")
            if not genus or NOT_A_GENUS.search(genus):
                stats["skipped_not_a_genus"] += 1
                continue
            effect = None
            column_s5 = s5_index.get(dataset)
            if column_s5 is None:
                stats["effect_absent"] += 1
            else:
                cell_row = s5_table.get(row[0]) or []
                cell = cell_row[column_s5].strip() if column_s5 < len(cell_row) else ""
                if not cell:
                    stats["effect_absent"] += 1
                else:
                    try:
                        candidate = float(cell)
                    except ValueError:
                        candidate = None
                    if candidate is None:
                        stats["effect_absent"] += 1
                    elif candidate in s5_sentinels:
                        # division-by-zero placeholder, not a measurement
                        stats["effect_sentinel_skipped"] += 1
                    else:
                        effect = candidate
                        stats["effect_loaded"] += 1
            significant.append((genus, lineage, q, effect))

        if not significant:
            stats["datasets_with_no_hits"] += 1
            continue
        stats["datasets_loaded"] += 1

        if DRY_RUN:
            stats["associations"] += len(significant)
            for genus, _, q, effect in significant:
                audit.append({"dataset": dataset, "study": accession, "disease": name,
                              "genus": genus, "q_value": abs(q),
                              "log2_fold_change": "" if effect is None else effect,
                              "direction": "enriched" if q > 0 else "depleted",
                              "source_file": "S1/S5"})
            continue

        study_id = study_for(accession, dataset)

        case_id = disease_id(name, mesh)
        if dataset not in comparison_cache:
            comparison_cache[dataset] = _upsert_comparison(
                conn, study_id, accession, health_id, HEALTH[0], case_id, name,
                positive_id=case_id, negative_id=health_id,
                method=METHOD, rank="genus",
                notes=(f"MicrobiomeHD dataset {dataset}; signed q-value, positive = "
                       f"higher in cases. Source accession "
                       f"{rec.get('accession') or 'not reported'}."))
        comparison_id = comparison_cache[dataset]

        cursor = conn.cursor()
        for genus, lineage, q, effect in significant:
            if genus not in taxon_cache:
                taxon = {"genus": genus, "species": "", "taxonomic_rank": "genus",
                         "source_database": "MicrobiomeHD"}
                taxon.update({k: v for k, v in lineage.items() if k != "genus"})
                # upsert_taxon returns (id, EXISTED) -- the flag is true when the
                # row was already there, not when it was created.
                taxon_cache[genus], existed = upsert_taxon(conn, taxon, overwrite=False)
                stats["taxa_already_present" if existed else "taxa_created"] += 1
            # effect_type is part of uq_taxon_disease_evidence, so a row whose
            # effect_type changes between runs -- 'q_value' before S5 was loaded,
            # 'log2_fold_change' after -- inserts a second row for the same
            # finding instead of updating the first. Clear any other effect_type
            # from this source for this triple so one finding keeps one row.
            effect_type = "log2_fold_change" if effect is not None else "q_value"
            cursor.execute(
                """DELETE FROM taxon_disease_associations
                   WHERE taxon_id = %s AND disease_id = %s AND comparison_id = %s
                     AND source_database = 'MicrobiomeHD' AND effect_type <> %s""",
                (taxon_cache[genus], case_id, comparison_id, effect_type),
            )
            stats["superseded_rows_removed"] += cursor.rowcount
            cursor.execute(
                """
                INSERT INTO taxon_disease_associations
                    (taxon_id, disease_id, comparison_id, direction, effect_type,
                     effect_size, q_value, source_database)
                VALUES (%s, %s, %s, %s, %s, %s, %s, 'MicrobiomeHD')
                ON DUPLICATE KEY UPDATE
                    direction = VALUES(direction), q_value = VALUES(q_value),
                    effect_size = VALUES(effect_size)
                """,
                (taxon_cache[genus], case_id, comparison_id,
                 "enriched" if q > 0 else "depleted",
                 effect_type, effect, abs(q)),
            )
            stats["associations"] += 1
            audit.append({"dataset": dataset, "study": accession, "disease": name,
                          "genus": genus, "q_value": abs(q),
                          "log2_fold_change": "" if effect is None else effect,
                          "direction": "enriched" if q > 0 else "depleted",
                          "source_file": "S1/S5"})
        cursor.close()

    # --- file-S4: what the original publications reported ---------------------
    lit, attrition, _lit_map = load_literature(datasets)
    if lit:
        print(f"\nS4 literature: {len(lit)} loadable rows from "
              f"{len({r['dataset'] for r in lit})} studies")
        for key, n in attrition.most_common():
            stats[f"s4_excluded_{key}"] = n

        genus_ok: dict[str, bool] = {}
        for row in lit:
            genus = row["genus"]
            if genus not in genus_ok:
                # Validate against the local dump: the file contains at least one
                # typo ("Peptosreptococcus"), and loading it would invent a genus.
                # Not auto-corrected -- a silent name fix is how the homonym
                # lineages got into the curated table in the first place.
                genus_ok[genus] = taxdump_has_prokaryote_genus(genus)
                if not genus_ok[genus]:
                    print(f"    !! {genus!r} does not resolve to a prokaryote genus - skipped")
            if not genus_ok[genus]:
                stats["s4_genus_unresolved"] += 1
                continue
            code = row["dataset"].split("_")[0]
            name, mesh = DISEASE_FOR_CODE[code]
            rec = identity.get(row["dataset"], {})
            accession = f"PMID{rec['pmid']}" if rec.get("pmid") else rec.get("accession") or row["dataset"]
            stats["s4_associations"] += 1
            organism = f"{genus} {row['species']}".strip()
            audit.append({"dataset": row["dataset"], "study": accession, "disease": name,
                          "genus": organism, "q_value": row["qval"], "log2_fold_change": "",
                          "direction": row["direction"], "source_file": "S4"})
            if DRY_RUN:
                continue
            case_id = disease_id(name, mesh)
            key = f"{row['dataset']}|{row['method']}"
            if key not in comparison_cache:
                comparison_cache[key] = _upsert_comparison(
                    conn, study_for(accession, row["dataset"]), accession, health_id, HEALTH[0],
                    case_id, name, positive_id=case_id, negative_id=health_id,
                    method=f"as reported ({row['method']})", rank="genus",
                    notes=(f"As reported in the publication for MicrobiomeHD dataset "
                           f"{row['dataset']}, curated in file-S4; multiple-comparison "
                           f"correction: {row['multi_comp'] or 'not stated'}. Distinct "
                           f"from the standardized re-analysis comparison on this study."))
            cache_key = f"{genus}|{row['species']}"
            if cache_key not in taxon_cache:
                taxon_cache[cache_key], existed = upsert_taxon(conn, {
                    "genus": genus, "species": row["species"],
                    "taxonomic_rank": row["rank"],
                    "source_database": "MicrobiomeHD"}, overwrite=False)
                stats["taxa_already_present" if existed else "taxa_created"] += 1
            cursor = conn.cursor()
            cursor.execute(
                """
                INSERT INTO taxon_disease_associations
                    (taxon_id, disease_id, comparison_id, direction, effect_type,
                     effect_size, q_value, source_database)
                VALUES (%s, %s, %s, %s, 'reported', NULL, %s, 'MicrobiomeHD')
                ON DUPLICATE KEY UPDATE
                    direction = VALUES(direction), q_value = VALUES(q_value)
                """,
                (taxon_cache[cache_key], case_id, comparison_cache[key],
                 row["direction"], row["qval"]),
            )
            cursor.close()

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    else:
        conn.commit()
        _finish_run(conn, run_id, LoadStats(
            read=stats["qvalues_present"],
            inserted=stats["associations"],
            updated=0,
            skipped=stats["not_significant"] + stats["skipped_not_a_genus"],
        ))

    for key in ("datasets_loaded", "datasets_with_no_hits", "datasets_skipped",
                "qvalues_present", "not_significant", "skipped_not_a_genus",
                "associations", "effect_loaded", "effect_sentinel_skipped",
                "effect_absent", "superseded_rows_removed",
                "taxa_created", "taxa_already_present",
                "s4_associations", "s4_genus_unresolved"):
        print(f"    {key:26} {stats[key]}")
    for key in sorted(k for k in stats if k.startswith("s4_excluded_")):
        print(f"    {key:26} {stats[key]}")

    if audit:
        AUDIT.parent.mkdir(parents=True, exist_ok=True)
        with AUDIT.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=[
                "source_file", "dataset", "study", "disease", "genus",
                "q_value", "log2_fold_change", "direction"])
            writer.writeheader()
            writer.writerows(audit)
        print(f"\naudit written to {AUDIT} ({len(audit)} rows)")

    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
