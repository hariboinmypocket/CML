#!/usr/bin/env python3
"""Load the BioMapAI ME/CFS cohort (PRJNA1125469) — Xiong et al., Nat Med 2025.

https://doi.org/10.1038/s41591-025-03788-3
Individual-level data published openly at github.com/ohlab/BioMapAI

479 runs, 249 subjects, with age/sex complete plus diet, antibiotic and
probiotic exposure. Longitudinal: most subjects contributed 2-3 timepoints,
so subject_id is populated - treating repeated measures as independent
observations would inflate the effective sample size and leak an individual
across train/test splits.

Abundances published as 135 species per sample, summing to a mean of 0.89
rather than 1.0 because the authors kept only those species. They are
aggregated to genus (matching this database's single-rank convention) and
renormalized, with the dropped share recorded in unclassified_fraction so
poorly-represented samples stay identifiable and filterable.

Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from gutdb.config import Settings
from gutdb.db import connect
from gutdb.pipeline import upsert_taxon
from gutdb.transform import parse_scientific_name

BASE = Path("/private/tmp/claude-501/-Users-lucyshaa-Desktop-CML/"
            "0a6c9862-7ad6-4ee9-a30e-151a04129642/scratchpad/biomapai/Omics_Dataset")
META = BASE / "Metadata.csv"
ABUND = BASE / "Specie Abundance.csv"
ENA = Path("/tmp/bio_ena.tsv")
PROJECT = "PRJNA1125469"
TITLE = "AI-driven multi-omics modeling of myalgic encephalomyelitis/chronic fatigue syndrome"
DRY_RUN = "--apply" not in sys.argv

DISEASE = {"MECFS": "Fatigue Syndrome, Chronic", "Control": "Health"}


def main() -> int:
    meta_rows = list(csv.DictReader(META.open(encoding="utf-8-sig")))
    idx_col = list(meta_rows[0].keys())[0]
    meta = {r[idx_col]: r for r in meta_rows}
    print(f"metadata samples: {len(meta)}")

    lib2run: dict[str, str] = {}
    with ENA.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            lib = re.sub(r"_lib$", "", (row.get("library_name") or "").strip())
            if lib:
                lib2run[lib] = row["run_accession"]
    print(f"SRA runs: {len(lib2run)}")

    matrix = list(csv.reader(ABUND.open(encoding="utf-8-sig")))
    header, species_rows = matrix[0], matrix[1:]
    col_of = {name: i for i, name in enumerate(header)}

    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    cur.execute("SELECT id, name FROM diseases WHERE name IN (%s, %s)",
                (DISEASE["MECFS"], DISEASE["Control"]))
    disease_id = {name: did for did, name in cur.fetchall()}

    cur.execute("SELECT id FROM studies WHERE project_accession=%s", (PROJECT,))
    row = cur.fetchone()
    if row:
        study_id = row[0]
    elif DRY_RUN:
        study_id = None
    else:
        cur.execute(
            """INSERT INTO studies (project_accession, title, data_type,
                                    source_database, data_quality)
               VALUES (%s, %s, %s, %s, %s)""",
            (PROJECT, TITLE, "mNGS", "PubMed", "qualified"))
        study_id = cur.lastrowid
        conn.commit()

    stats = {"samples": 0, "already_present": 0, "no_run": 0,
             "abundance_rows": 0, "low_retention": 0, "no_abundance_col": 0}

    for sample_key, m in meta.items():
        run = lib2run.get(sample_key)
        if not run:
            stats["no_run"] += 1
            continue
        cur.execute("SELECT id FROM samples WHERE run_accession=%s", (run,))
        if cur.fetchone():
            stats["already_present"] += 1
            continue

        subject = re.sub(r"_tp\d+$", "", sample_key)
        tp = (re.search(r"_(tp\d+)$", sample_key) or [None, None])[1]
        group = (m.get("study_ptorhc") or "").strip()
        sex_raw = (m.get("gender") or "").strip().casefold()
        sex = "male" if sex_raw == "male" else ("female" if sex_raw == "female" else None)
        try:
            age = float(m.get("age") or "")
        except ValueError:
            age = None
        if age is not None and not (0 <= age <= 122):
            age = None

        # aggregate the published species to genus, then renormalize
        genus_totals: dict[str, float] = {}
        gross = 0.0
        ci = col_of.get(sample_key)
        if ci is None:
            stats["no_abundance_col"] += 1
            continue
        for srow in species_rows:
            try:
                value = float(srow[ci])
            except (ValueError, IndexError):
                continue
            if value <= 0:
                continue
            gross += value
            genus, _ = parse_scientific_name(srow[0].replace("_", " "))
            if genus:
                genus_totals[genus] = genus_totals.get(genus, 0.0) + value
        if not genus_totals or gross <= 0:
            continue
        if gross < 0.7:
            stats["low_retention"] += 1

        stats["samples"] += 1
        if DRY_RUN:
            stats["abundance_rows"] += len(genus_totals)
            continue

        cur.execute(
            """INSERT INTO samples (study_id, disease_id, run_accession, subject_id,
                                    timepoint, sex, age_years, body_site,
                                    unclassified_fraction)
               VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)""",
            (study_id, disease_id.get(DISEASE.get(group, "")), run, subject, tp,
             sex, age, "feces", round(1.0 - gross, 6)))
        sid = cur.lastrowid

        total = sum(genus_totals.values())
        for genus, value in genus_totals.items():
            taxon_id, _ = upsert_taxon(conn, {
                "genus": genus, "species": "", "taxonomic_rank": "genus",
                "source_database": "PubMed"}, overwrite=False)
            cur.execute(
                """INSERT INTO sample_taxon_abundances (sample_id, taxon_id, relative_abundance)
                   VALUES (%s, %s, %s)
                   ON DUPLICATE KEY UPDATE relative_abundance = VALUES(relative_abundance)""",
                (sid, taxon_id, round(value / total, 8)))
            stats["abundance_rows"] += 1

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    else:
        conn.commit()

    for key, value in stats.items():
        print(f"  {key:18} {value}")
    cur.close()
    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
