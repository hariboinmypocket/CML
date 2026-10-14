#!/usr/bin/env python3
"""Fold post-2021 phylum names onto the classic name, everywhere they are stored.

Seven phyla were present in `taxa` under BOTH their classic and their 2021 ICNP
name at once -- 1,365 rows saying Firmicutes beside 97 saying Bacillota, 424
Bacteroidetes beside 12 Bacteroidota, and so on across 3,636 taxa. Nothing
errors on that. It just means `GROUP BY phylum` silently reports one phylum as
two, which is why PHYLUM_SYNONYMS exists and why the Firmicutes/Bacteroidetes
ratio, the most quoted number in this field, was being computed over an
incomplete set.

The names came in through loaders that read NCBI's lineage directly. NCBI's
SCIENTIFIC name for these is now the ICNP one -- taxid 1239 is Bacillota, with
Firmicutes recorded as a synonym -- so anything resolving a lineage from the dump
and storing the answer verbatim writes the new name. 119 of the 131 affected rows
came from scripts/load_muller2022.py, which now calls normalize_phylum at
creation; this repairs what is already stored.

Only renames are folded, never reclassifications. The test is the tax ID:
Euryarchaeota and Methanobacteriota are both taxid 28890, so one is a renaming of
the other, while Thermoplasmatota (2283796) and Halobacteriota (3100166) are
separate nodes that GTDB split out of the old Euryarchaeota. Moving taxa between
phyla is a taxonomic decision and belongs in normalize_reclassified_phyla.py's
APPROVED set, where a human signs off on it.

Run this BEFORE scripts/normalize_reclassified_phyla.py. That script compares the
stored phylum against NCBI's, but it normalizes only NCBI's side, so a stored
value carrying a post-2021 name or a typo is reported as a disagreement it cannot
act on. Folding first leaves that script looking at genuine moves between phyla,
which is the only thing it should be asked to judge.

Both stores are folded: MySQL and the SQLite mirror, which is built from the CSVs
in data/ and otherwise keeps whatever names those carried.

Idempotent. Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import os
import sqlite3
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.transform import normalize_phylum

DRY_RUN = "--apply" not in sys.argv
SQLITE = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")


def fold_store(label: str, cur, placeholder: str) -> int:
    """Fold one store's phylum column. Returns the number of rows rewritten."""
    cur.execute("SELECT phylum, COUNT(*) FROM taxa "
                "WHERE phylum IS NOT NULL AND phylum <> '' GROUP BY phylum")
    stored = cur.fetchall()
    counts = dict(stored)
    changes = [(p, normalize_phylum(p), n) for p, n in stored
               if normalize_phylum(p) and normalize_phylum(p) != p]
    if not changes:
        print(f"   {label}: already folded")
        return 0
    total = 0
    for old, new, n in sorted(changes, key=lambda c: -c[2]):
        total += n
        print(f"   {label}: {old:26} -> {new:24} {n:>5} rows "
              f"(joins {counts.get(new, 0)})")
        if not DRY_RUN:
            cur.execute(f"UPDATE taxa SET phylum = {placeholder} "
                        f"WHERE phylum = {placeholder}", (new, old))
    return total


def main() -> int:
    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()

    cur.execute("""SELECT phylum, COUNT(*) FROM taxa
                   WHERE phylum IS NOT NULL AND phylum <> '' GROUP BY phylum""")
    stored = cur.fetchall()

    changes: list[tuple[str, str, int]] = []
    for phylum, n in stored:
        folded = normalize_phylum(phylum)
        if folded and folded != phylum:
            changes.append((phylum, folded, n))

    if not changes:
        print("every stored phylum is already the classic name; nothing to fold")
        conn.close()
        return 0

    counts = dict(stored)
    print(f"{'stored name':26} {'folds to':24} {'rows':>6}  joins an existing group of")
    total = 0
    for old, new, n in sorted(changes, key=lambda c: -c[2]):
        total += n
        existing = counts.get(new, 0)
        print(f"   {old:23} {new:24} {n:>6}  {existing:>6}")
    print(f"\n{len(changes)} names, {total} rows")

    merged = sum(1 for old, new, _ in changes if counts.get(new))
    print(f"{merged} of them merge into a group that already exists, which is the "
          f"split being repaired")

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
        conn.close()
        return 0

    applied: Counter = Counter()
    for old, new, _ in changes:
        cur.execute("UPDATE taxa SET phylum = %s WHERE phylum = %s", (new, old))
        applied[old] = cur.rowcount
    conn.commit()
    print(f"\nmysql: updated {sum(applied.values())} rows")

    if Path(SQLITE).exists():
        print("\nsqlite mirror:")
        sq = sqlite3.connect(SQLITE)
        n = fold_store("sqlite", sq.cursor(), "?")
        if not DRY_RUN:
            sq.commit()
            print(f"   sqlite: updated {n} rows")
        sq.close()
    else:
        print(f"\nsqlite mirror {SQLITE} not present; skipped")

    cur.execute("""SELECT phylum, COUNT(*) FROM taxa
                   WHERE phylum IS NOT NULL AND phylum <> '' GROUP BY phylum""")
    remaining = [p for p, _ in cur.fetchall() if normalize_phylum(p) != p]
    print(f"stored names still needing folding: {len(remaining)}"
          + (f" -> {remaining}" if remaining else "  (clean)"))
    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
