#!/usr/bin/env python3
"""Stop the mirror keeping a METHOD name in effect_type.

effect_type records what the stored number IS -- 'LDA', 'log2_fold_change', or a
sentinel when there is no number -- not how it was produced. resolve_effect_type()
in gutdb/transform.py exists to make that conversion, and MySQL goes through it:
no MySQL row holds a method name there.

The mirror holds 123 rows reading 'LEfSe', which is a method. It is built by
scripts/build_gutdb_sqlite.py from the CSVs and never called the resolver, so the
raw declared method landed in the column -- the wrong-axis defect this schema
keeps producing, in the one store that skipped the shared function.

    122 have no effect_size -> 'reported'
      1 carries LDA 4.31    -> 'LDA'

Those are the resolver's own answers, not a rule invented here, and they bring the
mirror's vocabulary to MySQL's: LDA and reported, with no method names.

The MySQL-only values (log2_fold_change, consensus, q_value) are absent from the
mirror because the MicrobiomeHD tables they come from are not part of its build;
that is a membership difference, not a vocabulary one, and nothing here touches
it.

Idempotent. Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import os
import sqlite3
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.transform import EFFECT_TYPE_METHODS, resolve_effect_type

DRY_RUN = "--apply" not in sys.argv


def main() -> int:
    path = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")
    if not Path(path).exists():
        print(f"sqlite mirror {path} not present; nothing to do")
        return 0

    conn = sqlite3.connect(path)
    cur = conn.cursor()
    cur.execute("""SELECT effect_type, COUNT(*) FROM taxon_disease_associations
                   WHERE effect_type IS NOT NULL GROUP BY effect_type""")
    before = dict(cur.fetchall())
    offenders = [v for v in before
                 if str(v).strip().casefold() in EFFECT_TYPE_METHODS]
    if not offenders:
        print("no effect_type names a method; nothing to do")
        print(f"   mirror effect_type: {before}")
        conn.close()
        return 0

    print(f"effect_type values that name a method: {offenders}")
    counts: Counter = Counter()
    for declared in offenders:
        # Resolve per row, because the answer depends on whether that row has a
        # magnitude -- which is the whole point of the function.
        cur.execute("""SELECT id, effect_size FROM taxon_disease_associations
                       WHERE effect_type = ?""", (declared,))
        rows = cur.fetchall()
        updates = []
        for row_id, effect_size in rows:
            resolved = resolve_effect_type(declared, effect_size)
            counts[f"{declared} -> {resolved}"] += 1
            updates.append((resolved, row_id))
        if not DRY_RUN:
            cur.executemany(
                "UPDATE taxon_disease_associations SET effect_type = ? WHERE id = ?",
                updates)
    for key, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"   {key:34} {n:>5}")
    if not DRY_RUN:
        conn.commit()

    cur.execute("""SELECT effect_type, COUNT(*) FROM taxon_disease_associations
                   WHERE effect_type IS NOT NULL GROUP BY effect_type
                   ORDER BY COUNT(*) DESC""")
    print(f"\nmirror effect_type after: {dict(cur.fetchall())}")
    conn.close()

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
