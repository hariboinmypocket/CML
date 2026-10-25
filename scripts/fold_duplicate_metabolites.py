#!/usr/bin/env python3
"""Remove metabolite rows that duplicate another row of the same compound.

Two imports populated `metabolites`: one naming compounds with a KEGG prefix
(`C00042_Succinate`) and a later one using plain names (`glycolic acid`). Where
the later import re-added a compound already present, its KEGG id is now held
by several rows and `sample_metabolite_levels` is split across them, so a query
filtering on kegg_id either double-counts one compound or pools several.

ONLY THE INDEFENSIBLE CASE IS FOLDED: rows sharing both kegg_id AND an
identical name. A shared kegg_id alone is not duplication, because KEGG ids for
lipid classes name a class rather than a molecule -- C00157 is held by three
distinct phosphatidylcholine species (16:0/18:2, 16:0/18:1, 18:0/18:1) and
C04230 by two GPC species, which are different molecules that must not be
merged. Deciding the remaining cases needs chemistry per pair (p-Hydroxybenzoate
and 4-hydroxybenzoate are one compound under two HMDB ids; Citramalate and
Methylsuccinate may not be), so this does not guess at them. They stay visible
through v_metabolite_kegg_ambiguity.

Two conditions have to hold, and the first is the one doing the work. A
candidate must share both kegg_id and an identical name; then the guard
additionally requires every level row under the duplicate to be matched, in the
same sample, by a byte-identical value under the kept row.

VALUE IDENTITY ALONE WOULD NOT BE SAFE, which was measured rather than assumed.
Applying the guard by hand to C04230 -- `1-palmitoyl-GPC (16:0)` against
`_1_OLEOYLGLYCEROPHOSPHOCHOLIN`, a palmitoyl and an oleoyl species, so different
molecules -- reports all 131 level rows byte-identical and would fold them.
C02614, Citramalate against Methylsuccinate, reports all 648 identical and
would fold too. Identical numbers therefore do not establish that two rows are
the same compound. Anyone widening `candidates()` to "same kegg_id" and leaning
on the guard would merge distinct lipid species without a single refusal, so the
name match is a safety condition and not a convenience.

Preview by default; --apply writes an audit and deletes.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect

TOLERANCE = 1e-9


def candidates(connection) -> list[tuple[str, str, list[int]]]:
    cursor = connection.cursor()
    try:
        cursor.execute("""SELECT kegg_id, name, GROUP_CONCAT(id ORDER BY id)
                          FROM metabolites WHERE kegg_id IS NOT NULL
                          GROUP BY kegg_id, name HAVING COUNT(*) > 1""")
        return [(k, n, [int(x) for x in ids.split(",")])
                for k, n, ids in cursor.fetchall()]
    finally:
        cursor.close()


def redundant(connection, keep: int, drop: int) -> tuple[bool, str, int]:
    """Is every level under `drop` already stored identically under `keep`?"""
    cursor = connection.cursor()
    try:
        cursor.execute("SELECT COUNT(*) FROM sample_metabolite_levels WHERE metabolite_id = %s",
                       (drop,))
        total = cursor.fetchone()[0]
        cursor.execute("""SELECT COUNT(*), SUM(ABS(d.level - k.level) < %s)
                          FROM sample_metabolite_levels d
                          JOIN sample_metabolite_levels k ON k.sample_id = d.sample_id
                           AND k.metabolite_id = %s
                          WHERE d.metabolite_id = %s""", (TOLERANCE, keep, drop))
        matched, identical = cursor.fetchone()
        identical = int(identical or 0)
        if matched != total:
            return False, f"{total - matched} of {total} level rows have no counterpart", total
        if identical != total:
            return False, f"{total - identical} of {total} counterparts differ in value", total
        return True, f"all {total} level rows duplicate the kept row exactly", total
    finally:
        cursor.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env-file", default=".env")
    parser.add_argument("--audit-dir", type=Path, default=Path("data"))
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    connection = connect(Settings.from_env(args.env_file))
    try:
        plan, refused = [], []
        for kegg, name, ids in candidates(connection):
            keep, drops = ids[0], ids[1:]
            for drop in drops:
                ok, why, rows = redundant(connection, keep, drop)
                (plan if ok else refused).append((kegg, name, keep, drop, rows, why))
        for kegg, name, keep, drop, rows, why in plan:
            print(f"fold {kegg} {name!r}: drop id {drop} into {keep} -- {why}")
        for kegg, name, keep, drop, rows, why in refused:
            print(f"REFUSED {kegg} {name!r}: id {drop} vs {keep} -- {why}")
        if not plan:
            print("Nothing is provably redundant.")
            return
        print(f"\n{len(plan)} metabolite row(s) and "
              f"{sum(r for *_, r, _ in plan):,} duplicate level row(s) would go")
        if not args.apply:
            print("Preview only. Use --apply to fold.")
            return

        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        audit = args.audit_dir / f"metabolite_fold_audit_{stamp}.csv"
        cursor = connection.cursor()
        try:
            with audit.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["kegg_id", "name", "kept_metabolite_id",
                                 "dropped_metabolite_id", "sample_id", "level"])
                for kegg, name, keep, drop, rows, why in plan:
                    cursor.execute("""SELECT sample_id, level FROM sample_metabolite_levels
                                      WHERE metabolite_id = %s ORDER BY sample_id""", (drop,))
                    for sample_id, level in cursor.fetchall():
                        writer.writerow([kegg, name, keep, drop, sample_id, level])
            for kegg, name, keep, drop, rows, why in plan:
                cursor.execute("DELETE FROM sample_metabolite_levels WHERE metabolite_id = %s",
                               (drop,))
                cursor.execute("DELETE FROM metabolites WHERE id = %s", (drop,))
            connection.commit()
        finally:
            cursor.close()
        print(f"Folded. Audit: {audit.resolve()}")
    finally:
        connection.close()


if __name__ == "__main__":
    main()
