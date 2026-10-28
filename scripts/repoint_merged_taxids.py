#!/usr/bin/env python3
"""Repoint rows whose tax ID NCBI has retired, using merged.dmp.

38 taxa carried a tax ID absent from nodes.dmp, with 30,950 abundance rows and
186 associations hanging off them -- Waltera alone had 7,944. A retired ID joins
to nothing, which is what Klebsiella pneumoniae was doing with the deleted
1328388. `ncbi_rank IS NULL AND ncbi_tax_id IS NOT NULL` marks the class
exactly: 38 rows, no escapes and no false positives, because the rank backfill
cannot answer for an ID the dump has never heard of.

RESOLVED FROM merged.dmp, NOT FROM NAMES. Every one of the 38 also resolved to a
single candidate by name lookup, and the two methods agreed everywhere they
could be compared -- but only 15 matched a current scientific name and the other
23 matched a synonym or an `includes`, which proves a relationship and not an
identity. Tax 668 lists `Vibrio fischeri` among the synonyms of `Aliivibrio
fischeri`, and reading that as identity is how the Vibrio direction was got
wrong once already. merged.dmp states the replacement outright, so nothing here
is inferred from a name.

A REPOINT IS REFUSED WHERE THE TARGET IS ALREADY OCCUPIED. 13 of the 38 merge
onto an ID another row in the same abundance layer already holds, because the
database carries both the old and the new name as separate rows: Allisonella
(4,629 abundance rows) beside Dialister, Anaeroglobus (1,471) beside
Megasphaera, Raoultella beside Klebsiella, and Anaeroglobus and Colibacter both
merge onto Megasphaera. Those are organism merges, not ID corrections: the
per-(sample, rank) shares must be summed where both rows appear in one sample or
the sum to 1 breaks, which is the unwinding this repo has done once before by
hand. This script repoints the 25 that are unambiguous and names the rest.

ncbi_rank is set from nodes.dmp in the same pass, because a repoint that leaves
the rank describing the retired ID is the half-finished repair this project
keeps finding.

Preview by default; --apply writes both stores and an audit.
"""
from __future__ import annotations

import argparse
import csv
import sqlite3
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings, require_local
from gutdb.db import connect

SELECT_DEAD = ("SELECT id, scientific_name, ncbi_tax_id, taxonomic_rank, "
               "genus_key, species_key FROM taxa "
               "WHERE ncbi_tax_id IS NOT NULL AND ncbi_rank IS NULL")
SELECT_OCCUPIED = ("SELECT ncbi_tax_id, taxonomic_rank, id, scientific_name FROM taxa "
                   "WHERE ncbi_tax_id IS NOT NULL")


def read_merged(path: Path) -> dict[int, int]:
    out: dict[int, int] = {}
    with open(require_local(path), encoding="utf-8") as handle:
        for line in handle:
            parts = line.split("\t|")
            if len(parts) >= 2:
                out[int(parts[0])] = int(parts[1])
    return out


def read_ranks(path: Path) -> dict[int, str]:
    out: dict[int, str] = {}
    with open(require_local(path), encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split("\t|")
            out[int(parts[0])] = parts[2].strip()
    return out


class Store:
    def __init__(self, label, connection, placeholder):
        self.label, self.connection, self.ph = label, connection, placeholder

    def rows(self, sql, args=()):
        sql = sql.replace("%s", self.ph)
        if self.label == "sqlite":
            return self.connection.execute(sql, args).fetchall()
        cur = self.connection.cursor()
        try:
            cur.execute(sql, args)
            return cur.fetchall()
        finally:
            cur.close()

    def execute(self, sql, args=()):
        sql = sql.replace("%s", self.ph)
        if self.label == "sqlite":
            self.connection.execute(sql, args)
            return
        cur = self.connection.cursor()
        try:
            cur.execute(sql, args)
        finally:
            cur.close()


def plan(store, merged, ranks):
    dead = store.rows(SELECT_DEAD)
    occupied: dict[tuple[int, str], list] = {}
    for tax_id, layer, row_id, name in store.rows(SELECT_OCCUPIED):
        occupied.setdefault((int(tax_id), layer), []).append((row_id, name))

    wanted = Counter()
    for _, _, tax_id, _, _, _ in dead:
        if int(tax_id) in merged:
            wanted[merged[int(tax_id)]] += 1

    repoints, refusals = [], []
    for row_id, name, tax_id, layer, genus_key, species_key in dead:
        key = (genus_key, species_key)
        tax_id = int(tax_id)
        target = merged.get(tax_id)
        if target is None:
            refusals.append((key, f"{name!r}: tax ID {tax_id} is in neither "
                            f"merged.dmp nor a live node, so nothing states a "
                            f"replacement"))
            continue
        rank = ranks.get(target)
        if rank is None:
            refusals.append((key, f"{name!r}: merged.dmp sends {tax_id} to "
                            f"{target}, which is not a live node either"))
            continue
        others = [o for o in occupied.get((target, layer), []) if o[0] != row_id]
        if others:
            refusals.append((key, f"{name!r}: {tax_id} merges onto {target}, "
                            f"already held in the {layer} layer by "
                            f"{others[0][1]!r} -- an organism merge, whose shares "
                            f"must be summed, not a repoint"))
            continue
        if wanted[target] > 1:
            refusals.append((key, f"{name!r}: {tax_id} merges onto {target}, and "
                            f"so does another retired row -- they are one organism, "
                            f"to be folded rather than both repointed"))
            continue
        repoints.append((key, row_id, name, tax_id, target, rank))
    return repoints, refusals


def blocked_keys(plans: dict) -> dict[tuple[str, str], str]:
    """Keys no store may repoint, keyed by (genus_key, species_key).

    A key refused by ANY store is refused everywhere. The stores hold different
    subsets, so a target that collides in MySQL can look free in the mirror --
    13 refusals against 11 on the first run. Repointing on each store's own
    verdict would leave the two disagreeing about the organism, which is this
    project's recurring defect, manufactured by the script meant to fix one.

    Shared with tests/test_retired_taxids.py so the test cannot drift from the
    rule the script actually applies.
    """
    blocked: dict[tuple[str, str], str] = {}
    for label, (_, refusals) in plans.items():
        for key, message in refusals:
            blocked.setdefault(key, f"[{label}] {message}")
    return blocked


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env-file", default=".env")
    parser.add_argument("--sqlite", default="gut_microbiome.sqlite")
    parser.add_argument("--taxdump", type=Path, default=Path("taxdump"))
    parser.add_argument("--audit-dir", type=Path, default=Path("data"))
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    merged = read_merged(args.taxdump / "merged.dmp")
    ranks = read_ranks(args.taxdump / "nodes.dmp")
    print(f"merged.dmp: {len(merged):,} retired ids; nodes.dmp: {len(ranks):,} live\n")

    mysql = connect(Settings.from_env(args.env_file))
    lite = sqlite3.connect(args.sqlite, timeout=60)
    audit = []
    try:
        stores = [Store("mysql", mysql, "%s"), Store("sqlite", lite, "?")]
        plans = {store.label: plan(store, merged, ranks) for store in stores}

        blocked = blocked_keys(plans)

        for store in stores:
            repoints, _ = plans[store.label]
            allowed = [r for r in repoints if r[0] not in blocked]
            print(f"--- {store.label}: {len(allowed)} repointing, "
                  f"{len(repoints) - len(allowed)} held back by another store ---")
            for key, row_id, name, old, new, rank in allowed:
                print(f"  repoint {name!r}: {old} -> {new} [{rank}]")
                audit.append([store.label, row_id, name, old, new, rank])
            for key, row_id, name, old, new, rank in repoints:
                if key in blocked:
                    print(f"  HELD   {name!r}: safe here but {blocked[key]}")
            if args.apply:
                for key, row_id, name, old, new, rank in allowed:
                    store.execute("UPDATE taxa SET ncbi_tax_id = %s, ncbi_rank = %s "
                                  "WHERE id = %s", (new, rank, row_id))
                store.connection.commit()
            print()

        print("refused, with the reason recorded once per organism:")
        for key, message in sorted(blocked.items()):
            print(f"  {message}")
        print()
        if args.apply:
            args.audit_dir.mkdir(parents=True, exist_ok=True)
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            path = args.audit_dir / f"merged_taxid_repoint_audit_{stamp}.csv"
            with path.open("w", newline="", encoding="utf-8") as fh:
                writer = csv.writer(fh)
                writer.writerow(["store", "taxon_id", "scientific_name",
                                 "retired_tax_id", "replacement_tax_id", "ncbi_rank"])
                writer.writerows(audit)
            print(f"Applied. Audit: {path.resolve()}")
        else:
            print("Preview only. Use --apply to write both stores.")
    finally:
        lite.close()
        mysql.close()


if __name__ == "__main__":
    main()
