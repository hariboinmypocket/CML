#!/usr/bin/env python3
"""Fill empty lineage cells from the NCBI dump, in both stores.

Only EMPTY cells, and that restriction is the whole point. Stored lineage names
that disagree with the dump are already analysed in
reconcile_lineage_ranks.py, which separates a rename (the stored name is an
older spelling of the same node -- `Actinobacteria` and `Actinomycetes` are both
class 1760) from a reclassification (the name resolves to a different node), and
deliberately leaves the renames alone: adopting them would rewrite hundreds of
rows to a name the rest of the table does not use, which is a vocabulary
decision rather than a repair. Nothing here reopens that. A cell with no value
has no competing name, so filling it decides nothing.

100 such cells exist: 49 in MySQL and 51 in the mirror, across superkingdom,
phylum, class, order and family. Some became fillable only after the 114
one-sided tax IDs were filled -- a row with no tax ID has no lineage to look up,
so that repair made this one possible.

cleanup_gutdb.py already fills empty cells from the dump and, like the rest of
that script, reaches only the mirror. This reaches both.

A superkingdom value is written only when it is one of the four top-level names.
`Taxdump.lineage` resolves superkingdom through ("superkingdom", "domain",
"kingdom"), so for a virus it can return a kingdom -- `Heunggongvirae` for a row
whose correct superkingdom is `Viruses`. No empty cell hits that case today, and
the guard is here so that stays true.

Preview by default; --apply writes both stores and an audit.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect

COLUMNS = ["superkingdom", "phylum", "class_name", "order_name", "family"]
TOP_LEVEL = {"bacteria", "archaea", "eukaryota", "viruses"}


def load_taxdump(directory: str):
    """Reuse the canonical reader rather than parse the dump again here."""
    path = Path(__file__).resolve().parent / "fix_concatenated_genus_taxa.py"
    spec = importlib.util.spec_from_file_location("fix_concatenated_genus_taxa", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.Taxdump(directory)


def rows_of(connection, sqlite_mode):
    sql = ("SELECT id, scientific_name, ncbi_tax_id, "
           + ", ".join(COLUMNS) + " FROM taxa WHERE ncbi_tax_id IS NOT NULL")
    if sqlite_mode:
        return connection.execute(sql).fetchall()
    cursor = connection.cursor()
    try:
        cursor.execute(sql)
        return cursor.fetchall()
    finally:
        cursor.close()


def plan(connection, sqlite_mode, dump):
    fills, refusals = [], []
    for row in rows_of(connection, sqlite_mode):
        lineage = dump.lineage(int(row[2]))
        for index, column in enumerate(COLUMNS):
            if row[3 + index] not in (None, ""):
                continue
            value = lineage.get(column)
            if not value:
                continue
            if column == "superkingdom" and value.strip().casefold() not in TOP_LEVEL:
                refusals.append(f"{row[1]!r}: dump gives superkingdom {value!r}, "
                                f"which is not a top-level name -- the rank "
                                f"fallback has returned a kingdom")
                continue
            fills.append((row[0], row[1], column, value))
    return fills, refusals


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env-file", default=".env")
    parser.add_argument("--sqlite", default="gut_microbiome.sqlite")
    parser.add_argument("--taxdump", default="taxdump")
    parser.add_argument("--audit-dir", type=Path, default=Path("data"))
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    dump = load_taxdump(args.taxdump)
    mysql = connect(Settings.from_env(args.env_file))
    lite = sqlite3.connect(args.sqlite, timeout=60)
    audit = []
    try:
        for label, connection, sqlite_mode in (("mysql", mysql, False),
                                               ("sqlite", lite, True)):
            fills, refusals = plan(connection, sqlite_mode, dump)
            by_column: dict[str, int] = {}
            for _, _, column, _ in fills:
                by_column[column] = by_column.get(column, 0) + 1
            print(f"--- {label}: {len(fills)} cells to fill {by_column or ''} ---")
            for message in refusals:
                print(f"  REFUSED {message}")
            for row_id, name, column, value in fills[:6]:
                print(f"  {name!r}: {column} <- {value!r}")
            if len(fills) > 6:
                print(f"  ... and {len(fills) - 6} more")
            audit.extend([label, row_id, name, column, value]
                         for row_id, name, column, value in fills)
            if args.apply:
                for row_id, name, column, value in fills:
                    sql = (f"UPDATE taxa SET {column} = ? WHERE id = ?" if sqlite_mode
                           else f"UPDATE taxa SET {column} = %s WHERE id = %s")
                    if sqlite_mode:
                        lite.execute(sql, (value, row_id))
                    else:
                        cursor = mysql.cursor()
                        try:
                            cursor.execute(sql, (value, row_id))
                        finally:
                            cursor.close()
                connection.commit()
        if args.apply:
            args.audit_dir.mkdir(parents=True, exist_ok=True)
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            path = args.audit_dir / f"lineage_fill_audit_{stamp}.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["store", "taxon_id", "scientific_name",
                                 "column", "filled_with"])
                writer.writerows(audit)
            print(f"\nApplied. Audit: {path.resolve()}")
        else:
            print("\nPreview only. Use --apply to write both stores.")
    finally:
        lite.close()
        mysql.close()


if __name__ == "__main__":
    main()
