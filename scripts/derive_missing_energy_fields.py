#!/usr/bin/env python3
"""Derive energy_mode and primary_food_source where they are still 'N/A'.

Both columns are derived, not sourced: `_collapse_mimedb` calls
`derive_energy_mode` and `derive_primary_food_source` while loading MiMeDB
rows, and both default to 'N/A' for anything that arrives another way. 'N/A' is
a real member of the energy_mode enum, so this is not a stray sentinel to
clean -- it is the column saying "not known", on rows where the row's own
traits say otherwise.

Measured rather than assumed: 477 cells in MySQL and 2,535 in the mirror hold
'N/A' while the project's own derivation returns a definite answer from
oxygen_requirement, metabolism and energy_source -- values already in the same
row. `aerobe` gives respirator, `anaerobe` fermenter, `facultative anaerobe`
mixed.

Nothing is invented and nothing is overwritten. The functions are the project's,
called with the row's own cells, and only a cell reading 'N/A' is written. A
derivation that returns 'N/A' leaves the cell alone, so re-running changes
nothing.

The mirror can derive five times as many as MySQL, and that is not a bug here:
the mirror holds the trait columns that cleanup_gutdb.py back-filled into it
alone, so it has the evidence MySQL is missing. Filling MySQL's traits from
those is a separate and still-open question -- see the README -- and this does
not prejudge it. It derives from whatever each store actually holds.

Preview by default; --apply writes both stores and an audit.
"""
from __future__ import annotations

import argparse
import csv
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.transform import derive_energy_mode, derive_primary_food_source

# There is no `activity` column on taxa, so the derivations get "" for it --
# the same value the loader passes when a source row has no activity field.
SELECT = ("SELECT id, scientific_name, oxygen_requirement, metabolism, "
          "energy_source, energy_mode, primary_food_source FROM taxa "
          "WHERE energy_mode = 'N/A' OR primary_food_source = 'N/A'")
ENERGY_VALUES = {"fermenter", "respirator", "mixed"}


def rows_of(connection, sqlite_mode):
    if sqlite_mode:
        return connection.execute(SELECT).fetchall()
    cursor = connection.cursor()
    try:
        cursor.execute(SELECT)
        return cursor.fetchall()
    finally:
        cursor.close()


def plan(connection, sqlite_mode):
    fills, refusals = [], []
    for row_id, name, oxygen, metabolism, source, mode, food in rows_of(
            connection, sqlite_mode):
        if mode == "N/A":
            derived = derive_energy_mode(oxygen, metabolism, source, "")
            if derived and derived != "N/A":
                if derived not in ENERGY_VALUES:
                    refusals.append(f"{name!r}: derived energy_mode {derived!r} is "
                                    f"not in the enum")
                else:
                    fills.append((row_id, name, "energy_mode", derived))
        if food == "N/A":
            derived = derive_primary_food_source(metabolism, "", source)
            if derived and derived != "N/A":
                fills.append((row_id, name, "primary_food_source", derived))
    return fills, refusals


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env-file", default=".env")
    parser.add_argument("--sqlite", default="gut_microbiome.sqlite")
    parser.add_argument("--audit-dir", type=Path, default=Path("data"))
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    mysql = connect(Settings.from_env(args.env_file))
    lite = sqlite3.connect(args.sqlite, timeout=60)
    audit = []
    try:
        for label, connection, sqlite_mode in (("mysql", mysql, False),
                                               ("sqlite", lite, True)):
            fills, refusals = plan(connection, sqlite_mode)
            counts: dict[tuple[str, str], int] = {}
            for _, _, column, value in fills:
                counts[(column, value)] = counts.get((column, value), 0) + 1
            print(f"--- {label}: {len(fills)} cells derivable ---")
            for (column, value), n in sorted(counts.items()):
                print(f"  {column} <- {value!r}: {n}")
            for message in refusals:
                print(f"  REFUSED {message}")
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
            path = args.audit_dir / f"energy_field_derivation_audit_{stamp}.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["store", "taxon_id", "scientific_name",
                                 "column", "derived"])
                writer.writerows(audit)
            print(f"\nApplied. Audit: {path.resolve()}")
        else:
            print("\nPreview only. Use --apply to write both stores.")
    finally:
        lite.close()
        mysql.close()


if __name__ == "__main__":
    main()
