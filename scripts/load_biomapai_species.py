#!/usr/bin/env python3
"""Add the species layer for the BioMapAI cohort (PRJNA1125469).

These samples were first loaded genus-only, matching the database's single-rank
convention at the time. Now that a species layer exists alongside genus, the
published species resolution can be restored from the same local file — GMrepo
does not hold this study, so it is the only source.

Species and genus layers are independent: each sums to 1 within its own rank.
The genus rows already present are left untouched.

The authors published 135 species per sample summing to ~0.89 rather than 1.0,
so the same renormalization used for the genus layer applies here, and the
dropped share is already recorded in samples.unclassified_fraction.

Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.pipeline import upsert_taxon
from gutdb.transform import parse_scientific_name

ABUND = Path("/private/tmp/claude-501/-Users-lucyshaa-Desktop-CML/"
             "0a6c9862-7ad6-4ee9-a30e-151a04129642/scratchpad/biomapai/"
             "Omics_Dataset/Specie Abundance.csv")
ENA = Path("/tmp/bio_ena.tsv")
DRY_RUN = "--apply" not in sys.argv


def main() -> int:
    import re

    lib2run: dict[str, str] = {}
    with ENA.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            lib = re.sub(r"_lib$", "", (row.get("library_name") or "").strip())
            if lib:
                lib2run[lib] = row["run_accession"]

    matrix = list(csv.reader(ABUND.open(encoding="utf-8-sig")))
    header, species_rows = matrix[0], matrix[1:]
    print(f"species in file: {len(species_rows)}  sample columns: {len(header)-1}")

    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    cur.execute("SELECT run_accession, id FROM samples WHERE study_id = "
                "(SELECT id FROM studies WHERE project_accession='PRJNA1125469')")
    run2id = {r[0]: r[1] for r in cur.fetchall()}
    print(f"BioMapAI samples in DB: {len(run2id)}")

    stats = {"samples": 0, "rows": 0, "no_run": 0, "already_has_species": 0}
    taxon_cache: dict[tuple[str, str], int] = {}

    for col_index in range(1, len(header)):
        sample_key = header[col_index]
        run = lib2run.get(sample_key)
        if not run or run not in run2id:
            stats["no_run"] += 1
            continue
        sample_id = run2id[run]

        cur.execute(
            """SELECT 1 FROM sample_taxon_abundances a JOIN taxa t ON t.id=a.taxon_id
               WHERE a.sample_id=%s AND t.taxonomic_rank='species' LIMIT 1""",
            (sample_id,),
        )
        if cur.fetchone():
            stats["already_has_species"] += 1
            continue

        values: dict[tuple[str, str], float] = {}
        total = 0.0
        for srow in species_rows:
            try:
                value = float(srow[col_index])
            except (ValueError, IndexError):
                continue
            if value <= 0:
                continue
            genus, species = parse_scientific_name(srow[0].replace("_", " "))
            if not genus or not species:
                continue  # not a binomial; belongs to the genus layer
            key = (genus, species)
            values[key] = values.get(key, 0.0) + value
            total += value
        if not values or total <= 0:
            continue

        stats["samples"] += 1
        if DRY_RUN:
            stats["rows"] += len(values)
            continue

        for (genus, species), value in values.items():
            if (genus, species) not in taxon_cache:
                tid, _ = upsert_taxon(conn, {
                    "genus": genus, "species": species, "taxonomic_rank": "species",
                    "source_database": "PubMed"}, overwrite=False)
                taxon_cache[(genus, species)] = tid
            cur.execute(
                """INSERT INTO sample_taxon_abundances (sample_id, taxon_id, relative_abundance)
                   VALUES (%s,%s,%s)
                   ON DUPLICATE KEY UPDATE relative_abundance = VALUES(relative_abundance)""",
                (sample_id, taxon_cache[(genus, species)], round(value / total, 8)),
            )
            stats["rows"] += 1

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    else:
        conn.commit()
    for key, value in stats.items():
        print(f"  {key:20} {value}")
    cur.close()
    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
