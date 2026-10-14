#!/usr/bin/env python3
"""Repair the one taxa row that concatenated two genera across a phylum boundary.

MiMeDB record MMDBm0412197 is strain HER1398. Its `name` and `species` columns
both say "Sphingobacterium paucimobilis" and its `ncbi_tax_id` (1346330) agrees,
but its `genus`/`klass`/`phylum` columns say Sphingomonas / Gammaproteobacteria /
Proteobacteria -- MiMeDB filed a Bacteroidetes organism under a Proteobacteria
lineage. _collapse_mimedb builds the stored name as genus + species, so the two
disagreeing fields were concatenated into "Sphingomonas sphingobacterium
paucimobilis", a name for an organism that does not exist, carrying the lineage
of one genus and the tax ID of another.

The strain explains how MiMeDB got there: HER1398 was deposited as Sphingomonas
(then Pseudomonas) paucimobilis, and White and Suttle re-sequenced it in 2013 and
moved it to Sphingobacterium, keeping the epithet. MiMeDB took the new name and
the old lineage. NCBI resolves 1346330 as strain HER1398 under species 1385985
"Sphingobacterium paucimobilis", genus Sphingobacterium, phylum Bacteroidota --
which normalize_phylum folds to Bacteroidetes, the spelling used here.

So the organism meant is Sphingobacterium paucimobilis, attested by the source
record's name, its species and its tax ID; only the lineage columns, already
known to be corrupt, point at Sphingomonas. The row is renamed to NCBI's current
name for the tax ID with the Bacteroidetes lineage, which is what
merge_duplicate_taxa.py does for this corruption class ("a parser artifact that
kept both the old and new genus") whenever a row sharing the tax ID exists to
merge into. Here none does -- 1346330 is unique to this row -- so that script
never considered it, and the same canonicalization is applied on its own.

Not merged into Sphingomonas paucimobilis (row 1684, tax 13689): that asserts
the misidentification the 2013 re-sequencing corrected, moves Bacteroidetes
traits onto a Proteobacterium, and merges across two tax IDs, which
merge_duplicate_taxa.py refuses by design. Not deleted either: the row has no
dependents, so deleting is safe, but the source record is real and exported, its
traits came from it, and nothing else in the schema records MiMeDB's microbe_id,
so deletion would drop a reference organism rather than correct it.

The tax ID moves from the strain node 1346330 to its species node 1385985. A
rank='species' row carrying a strain tax ID is the mismatch OVERRIDES in
merge_duplicate_taxa.py exists to undo (F. prausnitzii, L. salivarius, both
"retaxid"), and it would leave this row looking like a rename candidate forever,
since NCBI's name for 1346330 keeps the strain designation. Both IDs resolve to
the same phylum, so the phylum reconciliation is satisfied either way. Note that
NCBI quotes the authority for 1385985, meaning the species name is not validly
published; it is kept because the tax ID is this database's canonical key.

Traits are left as ingested -- MiMeDB's own claims about this strain.

Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect

DRY_RUN = "--apply" not in sys.argv
SQLITE_PATH = "gut_microbiome.sqlite"

# The corrupt row, identified by what makes it corrupt rather than by id.
CORRUPT_KEY = ("sphingomonas", "sphingobacterium paucimobilis")
CORRUPT_TAXID = 1346330

# NCBI's lineage for species 1385985, in the spellings already used by the
# Sphingobacterium rows in this table (rows 364, 8751, 2012115).
REPAIR = {
    "scientific_name": "Sphingobacterium paucimobilis",
    "genus": "Sphingobacterium",
    "species": "paucimobilis",
    "genus_key": "sphingobacterium",
    "species_key": "paucimobilis",
    "superkingdom": "Bacteria",
    "phylum": "Bacteroidetes",
    "class_name": "Sphingobacteriia",
    "order_name": "Sphingobacteriales",
    "family": "Sphingobacteriaceae",
    "ncbi_tax_id": 1385985,
}
COLUMNS = tuple(REPAIR)

# Dependents are checked rather than repointed: this repair changes which
# organism the row denotes, so evidence attached to it would have been filed
# against the wrong name and would need re-adjudication, not moving.
DEPENDENTS = (
    ("sample_taxon_abundances", "taxon_id"),
    ("taxon_disease_associations", "taxon_id"),
    ("taxon_disease_adjudication", "taxon_id"),
)


def repair(store: str, cursor, tables: set[str]) -> int:
    """Repair the row in one store. Returns 1 if written (or would be)."""
    ph = "%s" if store == "mysql" else "?"
    cursor.execute(
        f"SELECT id, scientific_name, phylum, ncbi_tax_id FROM taxa "
        f"WHERE genus_key={ph} AND species_key={ph}",
        CORRUPT_KEY,
    )
    found = cursor.fetchall()
    if not found:
        cursor.execute(
            f"SELECT id FROM taxa WHERE genus_key={ph} AND species_key={ph}",
            (REPAIR["genus_key"], REPAIR["species_key"]),
        )
        already = cursor.fetchone()
        print(f"  corrupt row absent"
              + (f"; repaired row present as id {already[0]}" if already
                 else "; nothing matching either key"))
        return 0

    written = 0
    for taxon_id, name, phylum, taxid in found:
        print(f"  id {taxon_id}: {name!r}  phylum={phylum}  tax_id={taxid}")
        if taxid != CORRUPT_TAXID:
            print(f"    SKIPPED: expected tax id {CORRUPT_TAXID}, found {taxid}. "
                  "The row is not the one this script documents.")
            continue

        blocked = False
        for table, column in DEPENDENTS:
            if table not in tables:
                continue
            cursor.execute(
                f"SELECT COUNT(*) FROM `{table}` WHERE `{column}`={ph}", (taxon_id,))
            count = cursor.fetchone()[0]
            print(f"    {table:<28} {count}")
            if count:
                blocked = True
        if blocked:
            print("    SKIPPED: row has dependents. Renaming it would relabel "
                  "existing evidence; re-adjudicate instead.")
            continue

        cursor.execute(
            f"SELECT id FROM taxa WHERE genus_key={ph} AND species_key={ph} "
            f"AND id<>{ph}",
            (REPAIR["genus_key"], REPAIR["species_key"], taxon_id),
        )
        clash = cursor.fetchone()
        if clash:
            print(f"    SKIPPED: id {clash[0]} already holds "
                  f"({REPAIR['genus_key']}, {REPAIR['species_key']}); "
                  "this is a merge, not a rename.")
            continue

        cursor.execute(
            f"SELECT id, scientific_name FROM taxa WHERE ncbi_tax_id={ph} AND id<>{ph}",
            (REPAIR["ncbi_tax_id"], taxon_id),
        )
        taxid_clash = cursor.fetchone()
        if taxid_clash:
            print(f"    SKIPPED: id {taxid_clash[0]} ({taxid_clash[1]!r}) already "
                  f"carries tax id {REPAIR['ncbi_tax_id']}.")
            continue

        print(f"    -> {REPAIR['scientific_name']!r}  "
              f"phylum={REPAIR['phylum']}  tax_id={REPAIR['ncbi_tax_id']}")
        if not DRY_RUN:
            assignments = ", ".join(f"{column}={ph}" for column in COLUMNS)
            cursor.execute(
                f"UPDATE taxa SET {assignments} WHERE id={ph}",
                tuple(REPAIR[column] for column in COLUMNS) + (taxon_id,),
            )
        written += 1
    return written


def main() -> int:
    total = 0

    print("mysql")
    conn = connect(Settings.from_env(".env"))
    cursor = conn.cursor()
    cursor.execute("SHOW TABLES")
    tables = {row[0] for row in cursor.fetchall()}
    total += repair("mysql", cursor, tables)
    if not DRY_RUN:
        conn.commit()
    cursor.close()
    conn.close()

    print("\nsqlite")
    if Path(SQLITE_PATH).exists():
        lite = sqlite3.connect(SQLITE_PATH)
        cursor = lite.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = {row[0] for row in cursor.fetchall()}
        total += repair("sqlite", cursor, tables)
        if not DRY_RUN:
            lite.commit()
        cursor.close()
        lite.close()
    else:
        print(f"  (no SQLite mirror at {SQLITE_PATH}; skipped)")

    if DRY_RUN:
        print(f"\n*** DRY RUN - nothing written ({total} row(s) would change). "
              "Re-run with --apply ***")
    else:
        print(f"\nrows repaired: {total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
