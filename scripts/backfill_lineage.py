#!/usr/bin/env python3
"""Fill in lineage ranks that were never populated, from the local NCBI dump.

1,417 taxa carry no phylum and 1,563 no superkingdom, and 1,415 of those are not
dead rows -- they have associations or abundance behind them. Waltera alone
accounts for 7,944 abundance rows. Most are recently described gut organisms that
arrived from GMrepo with a name and a tax ID but no lineage, so anything grouped
by phylum silently drops them.

ONLY EMPTY FIELDS ARE FILLED. A populated value is never overwritten, and that is
the whole safety property of this script rather than a nicety: several lineages
in this database were repaired by hand after cross-kingdom homonyms put Bosea
under Streptophyta and Bogoriella under Ascomycota, and a backfill that trusted
the dump over the stored value would quietly undo that work the next time someone
ran it.

The phylum goes through normalize_phylum, because NCBI's scientific name for
these is the 2021 ICNP one while this schema keeps the classic name so a phylum
groups as one value.

Idempotent. Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import os
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from fix_curated_csv_lineage import load_taxdump

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.transform import normalize_phylum, strip_nomenclatural_code

DRY_RUN = "--apply" not in sys.argv

# Database column -> NCBI rank. genus and species are deliberately absent: they
# are keys here, and rewriting them would change identity rather than fill a gap.
COLUMN_RANK = {
    "superkingdom": ("superkingdom", "domain"),
    "phylum": ("phylum",),
    "class_name": ("class",),
    "order_name": ("order",),
    "family": ("family",),
}


def main() -> int:
    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    columns = list(COLUMN_RANK)
    blank = " OR ".join(f"({c} IS NULL OR {c} = '')" for c in columns)
    cur.execute(
        f"""SELECT id, scientific_name, ncbi_tax_id, {', '.join(columns)}
            FROM taxa WHERE ncbi_tax_id IS NOT NULL AND ({blank})""")
    rows = cur.fetchall()
    print(f"taxa with a tax ID and at least one empty lineage rank: {len(rows):,}")
    if not rows:
        conn.close()
        return 0

    print("reading the NCBI dump...")
    _, lineages = load_taxdump()

    filled: Counter = Counter()
    unresolved = 0
    updates: list[tuple] = []
    for row in rows:
        taxon_id, name, taxid = row[0], row[1], int(row[2])
        current = dict(zip(columns, row[3:]))
        lineage = lineages.get(taxid)
        if not lineage:
            unresolved += 1
            continue
        changes = {}
        for column, ranks in COLUMN_RANK.items():
            if (current[column] or "").strip():
                continue                      # populated: never touched
            value = next((lineage.get(r) for r in ranks if lineage.get(r)), "")
            if not value:
                continue
            value = (normalize_phylum(value) if column == "phylum"
                     else strip_nomenclatural_code(value))
            if value:
                changes[column] = value
                filled[column] += 1
        if changes:
            updates.append((taxon_id, name, changes))

    print(f"   resolvable in the dump: {len(rows) - unresolved:,}"
          f"   unresolvable: {unresolved:,}")
    print(f"   taxa that would gain at least one rank: {len(updates):,}\n")
    print(f"   {'column':16} {'filled':>8}")
    for column in columns:
        print(f"   {column:16} {filled[column]:>8,}")

    if updates:
        print("\n   examples:")
        for taxon_id, name, changes in updates[:5]:
            shown = ", ".join(f"{k}={v}" for k, v in list(changes.items())[:3])
            print(f"      {str(name)[:36]:38} {shown}")

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
        conn.close()
        return 0

    for taxon_id, _name, changes in updates:
        sets = ", ".join(f"{c} = %s" for c in changes)
        cur.execute(f"UPDATE taxa SET {sets} WHERE id = %s",
                    (*changes.values(), taxon_id))
    conn.commit()
    print(f"\nupdated {len(updates):,} taxa")

    cur.execute("SELECT COUNT(*) FROM taxa WHERE phylum IS NULL OR phylum = ''")
    print(f"   taxa still without a phylum: {cur.fetchone()[0]:,}")
    cur.execute("SELECT COUNT(*) FROM taxa WHERE superkingdom IS NULL OR superkingdom = ''")
    print(f"   taxa still without a superkingdom: {cur.fetchone()[0]:,}")
    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
