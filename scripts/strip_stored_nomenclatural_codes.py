#!/usr/bin/env python3
"""Strip (SeqCode)/(ICNP) annotations already written into stored lineage ranks.

NCBI annotates names published outside the ICNP with their nomenclatural code,
so the dump carries "Methylomirabilota (SeqCode)" and siblings at the class,
order and family ranks. transform.strip_nomenclatural_code drops the annotation
on read, but it is reached only through normalize_phylum, and the taxdump fill in
cleanup_gutdb.py normalized the phylum column alone -- so every other rank stored
the annotation verbatim. The SQLite mirror holds 22 such values across phylum,
class_name, order_name and family; MySQL holds none.

The code is metadata about how a name was published, not part of the name. Left
attached it splits a rank into two values the moment anything arrives under the
bare name, which is the same defect PHYLUM_SYNONYMS exists to prevent.

cleanup_gutdb.py now strips at every rank, so the annotation stops arriving. This
script cleans what is already stored. It is deliberately separate from a mirror
rebuild: rebuilding from the CSVs alone drops the enrichment the live mirror has
accumulated (superkingdom 5,278 -> 0, phylum 5,268 -> 2,270, gram_stain
2,278 -> 0), so the stored values are repaired in place instead.

Only values that actually change are written, so a re-run is a no-op.

Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import os
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.transform import strip_nomenclatural_code

DRY_RUN = "--apply" not in sys.argv
SQLITE_PATH = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")

# Every lineage rank the taxdump fill writes, plus genus, which it does not
# annotate today but would be stored the same way if NCBI ever did.
COLUMNS = ("superkingdom", "phylum", "class_name", "order_name", "family", "genus")


def strip(store: str, cursor) -> int:
    ph = "%s" if store == "mysql" else "?"
    changed = 0
    for column in COLUMNS:
        cursor.execute(
            f"SELECT id, `{column}` FROM taxa "
            f"WHERE `{column}` LIKE {ph} OR `{column}` LIKE {ph}",
            ("%(SeqCode)%", "%(ICNP)%"),
        )
        for taxon_id, value in cursor.fetchall():
            stripped = strip_nomenclatural_code(value)
            if stripped == value or not stripped:
                continue
            print(f"    id {taxon_id:<8} {column:<13} {value!r} -> {stripped!r}")
            if not DRY_RUN:
                cursor.execute(
                    f"UPDATE taxa SET `{column}`={ph} WHERE id={ph}",
                    (stripped, taxon_id),
                )
            changed += 1
    if not changed:
        print("    nothing to strip")
    return changed


def main() -> int:
    total = 0

    print("mysql")
    conn = connect(Settings.from_env(".env"))
    cursor = conn.cursor()
    total += strip("mysql", cursor)
    if not DRY_RUN:
        conn.commit()
    cursor.close()
    conn.close()

    print("\nsqlite")
    if Path(SQLITE_PATH).exists():
        lite = sqlite3.connect(SQLITE_PATH)
        cursor = lite.cursor()
        total += strip("sqlite", cursor)
        if not DRY_RUN:
            lite.commit()
        cursor.close()
        lite.close()
    else:
        print(f"  (no SQLite mirror at {SQLITE_PATH}; skipped)")

    if DRY_RUN:
        print(f"\n*** DRY RUN - nothing written ({total} value(s) would change). "
              "Re-run with --apply ***")
    else:
        print(f"\nvalues stripped: {total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
