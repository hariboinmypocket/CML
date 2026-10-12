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

Idempotent. Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.transform import normalize_phylum

DRY_RUN = "--apply" not in sys.argv


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
    print(f"\nupdated {sum(applied.values())} rows")

    cur.execute("""SELECT phylum, COUNT(*) FROM taxa
                   WHERE phylum IS NOT NULL AND phylum <> '' GROUP BY phylum""")
    remaining = [p for p, _ in cur.fetchall() if normalize_phylum(p) != p]
    print(f"stored names still needing folding: {len(remaining)}"
          + (f" -> {remaining}" if remaining else "  (clean)"))
    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
