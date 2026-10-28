#!/usr/bin/env python3
"""Fill a tax ID that one store has and the other is missing.

120 rows exist in both stores under the same (genus_key, species_key) while
only one of them carries an ncbi_tax_id: 89 where MySQL has the ID and the
mirror has NULL, 31 the other way. Such a row can be joined to NCBI in one
store and not in the other -- the same organism, differently reachable.

NOTHING WATCHED THIS, INCLUDING THE TEST WRITTEN TO WATCH FOR EXACTLY THIS KIND
OF THING. tests/test_cross_store_tax_ids.py selects `WHERE ncbi_tax_id IS NOT
NULL`, so a key whose ID is missing on one side never enters both maps and is
never compared. The gap was in the guard, not only in the data.

The donor store's ID is not copied on trust. For each row the stored name must
appear among NCBI's names for that ID -- scientific name, synonym, equivalent
name or `includes` -- and the ID must be a live node. 119 of the 120 pass that
check; the one that does not is refused and named. A fill is also refused if it
would put two rows of the same abundance layer on one tax ID, which is what
taxa_sharing_a_tax_id_within_one_rank exists to catch.

ncbi_rank is written from nodes.dmp in the same statement. Filling the ID alone
would leave the rank NULL, which is the marker this repo uses for "points at an
ID NCBI does not have" -- so a half-filled row would masquerade as a retired one.

Preview by default; --apply writes both stores and an audit.
"""
from __future__ import annotations

import argparse
import csv
import re
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings, require_local
from gutdb.db import connect

IDENTITY = {"scientific name", "synonym", "equivalent name", "genbank synonym",
            "includes"}
SELECT_TAXA = ("SELECT genus_key, species_key, scientific_name, ncbi_tax_id, id, "
               "taxonomic_rank FROM taxa")


def norm(text):
    text = (text or "").strip().lower().replace("[", "").replace("]", "")
    return re.sub(r"\s+", " ", text)


def read_dump(taxdump: Path, wanted: set[int]):
    ranks: dict[int, str] = {}
    with open(require_local(taxdump / "nodes.dmp"), encoding="utf-8",
              errors="replace") as handle:
        for line in handle:
            parts = line.split("\t|")
            tax_id = int(parts[0])
            if tax_id in wanted:
                ranks[tax_id] = parts[2].strip()
    names: dict[int, set[str]] = {}
    with open(require_local(taxdump / "names.dmp"), encoding="utf-8",
              errors="replace") as handle:
        for line in handle:
            parts = line.split("\t|")
            if len(parts) < 4:
                continue
            tax_id = int(parts[0])
            if tax_id in wanted and parts[3].strip() in IDENTITY:
                names.setdefault(tax_id, set()).add(norm(parts[1]))
    return ranks, names


def load(connection, sqlite_mode):
    if sqlite_mode:
        rows = connection.execute(SELECT_TAXA).fetchall()
    else:
        cursor = connection.cursor()
        try:
            cursor.execute(SELECT_TAXA)
            rows = cursor.fetchall()
        finally:
            cursor.close()
    return {(a, b): {"name": n, "tax_id": t, "id": i, "layer": r}
            for a, b, n, t, i, r in rows}


def analyse(mysql, lite, taxdump: Path):
    """(fills, refusals) for the one-sided tax IDs, as main() decides them.

    Shared with tests/test_cross_store_tax_ids.py so the test cannot drift from
    the rule the script applies.
    """
    stores = {"mysql": load(mysql, False), "sqlite": load(lite, True)}
    shared = set(stores["mysql"]) & set(stores["sqlite"])
    gaps = []
    for key in sorted(shared):
        a, b = stores["mysql"][key], stores["sqlite"][key]
        if a["tax_id"] is None and b["tax_id"] is not None:
            gaps.append(("mysql", key, int(b["tax_id"])))
        elif b["tax_id"] is None and a["tax_id"] is not None:
            gaps.append(("sqlite", key, int(a["tax_id"])))
    ranks, names = read_dump(taxdump, {t for _, _, t in gaps})
    occupied = {label: {(int(v["tax_id"]), v["layer"])
                        for v in rows.values() if v["tax_id"] is not None}
                for label, rows in stores.items()}
    fills, refusals = [], []
    for target, key, tax_id in gaps:
        row = stores[target][key]
        if tax_id not in ranks:
            refusals.append((key, f"{row['name']!r}: donor ID {tax_id} is not a "
                             f"live node, so it would fill in a broken join"))
            continue
        if norm(row["name"]) not in names.get(tax_id, set()):
            refusals.append((key, f"{row['name']!r}: NCBI does not list that name "
                             f"for {tax_id}, so the donor ID is not confirmed"))
            continue
        if (tax_id, row["layer"]) in occupied[target]:
            refusals.append((key, f"{row['name']!r}: {target} already has a "
                             f"{row['layer']}-layer row on {tax_id}"))
            continue
        fills.append((target, key, row, tax_id, ranks[tax_id]))
        occupied[target].add((tax_id, row["layer"]))
    return fills, refusals, len(gaps)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env-file", default=".env")
    parser.add_argument("--sqlite", default="gut_microbiome.sqlite")
    parser.add_argument("--taxdump", type=Path, default=Path("taxdump"))
    parser.add_argument("--audit-dir", type=Path, default=Path("data"))
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    mysql = connect(Settings.from_env(args.env_file))
    lite = sqlite3.connect(args.sqlite, timeout=60)
    try:
        fills, refusals, total = analyse(mysql, lite, args.taxdump)
        print(f"one-sided tax IDs on shared keys: {total}")
        print(f"  fillable: {len(fills)}   refused: {len(refusals)}\n")
        from collections import Counter
        print("  by store:", dict(Counter(t for t, _, _, _, _ in fills)))
        for _key, message in refusals:
            print(f"  REFUSED {message}")

        if args.apply and fills:
            for target, key, row, tax_id, rank in fills:
                sql = ("UPDATE taxa SET ncbi_tax_id = ?, ncbi_rank = ? WHERE id = ?"
                       if target == "sqlite" else
                       "UPDATE taxa SET ncbi_tax_id = %s, ncbi_rank = %s WHERE id = %s")
                if target == "sqlite":
                    lite.execute(sql, (tax_id, rank, row["id"]))
                else:
                    cursor = mysql.cursor()
                    try:
                        cursor.execute(sql, (tax_id, rank, row["id"]))
                    finally:
                        cursor.close()
            mysql.commit()
            lite.commit()
            args.audit_dir.mkdir(parents=True, exist_ok=True)
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            path = args.audit_dir / f"one_sided_taxid_audit_{stamp}.csv"
            with path.open("w", newline="", encoding="utf-8") as fh:
                writer = csv.writer(fh)
                writer.writerow(["filled_store", "taxon_id", "scientific_name",
                                 "tax_id", "ncbi_rank"])
                writer.writerows([[t, r["id"], r["name"], x, k]
                                  for t, _, r, x, k in fills])
            print(f"\nApplied. Audit: {path.resolve()}")
        elif not args.apply:
            print("\nPreview only. Use --apply to write both stores.")
    finally:
        lite.close()
        mysql.close()


if __name__ == "__main__":
    main()
