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
from gutdb.transform import normalize_phylum

DRY_RUN = "--apply" not in sys.argv
# Upstream files are not redistributed here; point MICROBIOMEHD_DIR at a copy or
# fetch them into the default location (see FETCH_HINT below).
SOURCE = Path(os.environ.get("MICROBIOMEHD_DIR", "data/microbiomehd"))
QVALUES = SOURCE / "file-S1.qvalues.txt"
IDENTITY = SOURCE / "dataset_identity.csv"
AUDIT = Path("data/microbiomehd_s1_audit.csv")

RAW = "https://raw.githubusercontent.com/cduvallet/microbiomeHD/master"
FETCH_HINT = f"""
Required files are missing from {SOURCE}/. To fetch them:

  mkdir -p {SOURCE}
  curl -sS -o {SOURCE}/file-S1.qvalues.txt \\
    {RAW}/final/supp-files/file-S1.qvalues.txt
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


def load_identity() -> dict[str, dict[str, str]]:
    return {r["dataset"]: r for r in csv.DictReader(IDENTITY.open())}


def main() -> int:
    for required in (QVALUES, IDENTITY):
        if not required.exists():
            raise SystemExit(FETCH_HINT)
    identity = load_identity()
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
            significant.append((genus, lineage, q))

        if not significant:
            stats["datasets_with_no_hits"] += 1
            continue
        stats["datasets_loaded"] += 1

        if DRY_RUN:
            stats["associations"] += len(significant)
            for genus, _, q in significant:
                audit.append({"dataset": dataset, "study": accession, "disease": name,
                              "genus": genus, "q_value": abs(q),
                              "direction": "enriched" if q > 0 else "depleted"})
            continue

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
        study_id = study_cache[accession]

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
        for genus, lineage, q in significant:
            if genus not in taxon_cache:
                taxon = {"genus": genus, "species": "", "taxonomic_rank": "genus",
                         "source_database": "MicrobiomeHD"}
                taxon.update({k: v for k, v in lineage.items() if k != "genus"})
                # upsert_taxon returns (id, EXISTED) -- the flag is true when the
                # row was already there, not when it was created.
                taxon_cache[genus], existed = upsert_taxon(conn, taxon, overwrite=False)
                stats["taxa_already_present" if existed else "taxa_created"] += 1
            cursor.execute(
                """
                INSERT INTO taxon_disease_associations
                    (taxon_id, disease_id, comparison_id, direction, effect_type,
                     effect_size, q_value, source_database)
                VALUES (%s, %s, %s, %s, 'q_value', NULL, %s, 'MicrobiomeHD')
                ON DUPLICATE KEY UPDATE
                    direction = VALUES(direction), q_value = VALUES(q_value)
                """,
                (taxon_cache[genus], case_id, comparison_id,
                 "enriched" if q > 0 else "depleted", abs(q)),
            )
            stats["associations"] += 1
            audit.append({"dataset": dataset, "study": accession, "disease": name,
                          "genus": genus, "q_value": abs(q),
                          "direction": "enriched" if q > 0 else "depleted"})
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
                "associations", "taxa_created", "taxa_already_present"):
        print(f"    {key:26} {stats[key]}")

    if audit:
        AUDIT.parent.mkdir(parents=True, exist_ok=True)
        with AUDIT.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(audit[0].keys()))
            writer.writeheader()
            writer.writerows(audit)
        print(f"\naudit written to {AUDIT} ({len(audit)} rows)")

    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
