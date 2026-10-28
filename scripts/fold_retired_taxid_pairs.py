#!/usr/bin/env python3
"""Fold the retired-tax-ID rows whose merge needs no share arithmetic.

repoint_merged_taxids.py repointed the 25 retired IDs whose replacement was
free and refused 13 whose replacement is already held, in the same abundance
layer, by the row for the same organism under its newer name. Those 13 are not
one problem. Measured per pair -- samples holding both rows, and abundance and
association counts on each side -- they fall into four cases, and three of them
need no arithmetic at all:

  retired row holds no evidence        2  Sphingobium japonicum, Streptococcus sp. i-p16
  the occupant holds no evidence       2  Merdimmobilis hominis, Anaeroarcus
  both hold evidence, never co-occur   2  Defluviimonas, Massilibacillus
  the two co-occur in samples          7  Allisonella, Anaeroglobus, Raoultella, ...

THE FOURTH CASE IS REFUSED AND THAT IS THE POINT. Allisonella and Dialister
share 2,869 samples, Anaeroglobus and Megasphaera 1,195, Raoultella and
Klebsiella 310. Folding those means SUMMING the per-(sample, rank) shares
wherever both appear, because each is a real share of one sample's mass and the
layer must still sum to 1 -- the unwinding this repo has done once before by
hand. It also collapses two ML features into one, which changes any model
already trained on them. That is a data decision, so this states it and stops.

The other three cases move no mass. Where only one side holds evidence the
other row is a dead name and is dropped. Where both hold evidence but no sample
holds both, each sample still ends with exactly one row for the organism, so the
layer's sum is untouched -- only which taxon_id carries it changes.

Traits and the tax ID are carried onto the survivor before anything is deleted,
because a row being dropped here is often the only one holding them -- that was
the near-miss in merge_doubled_token_taxa.py. Rows are matched on
(genus_key, species_key), never on id, and an action is applied only if every
store that holds the row classifies it the same way: the mirror is a subset, so
deciding per store is how a one-store repair gets made by accident.

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

from gutdb.config import Settings, require_local
from gutdb.db import connect

CARRIED = [
    "ncbi_tax_id", "ncbi_rank", "gram_stain", "oxygen_requirement",
    "ph_preference", "sporulation", "shape", "cell_arrangement", "mobility",
    "flagella_presence", "number_of_membranes", "biotic_relationship",
    "habitat", "temperature_range", "optimal_temperature", "metabolism",
    "energy_source", "human_pathogen",
]
DROP_RETIRED = "retired row holds no evidence"
DROP_OCCUPANT = "the occupant holds no evidence"
MOVE_ROWS = "both hold evidence, never co-occur"
REFUSE_SUM = "the two co-occur, so shares would have to be summed"


class Store:
    def __init__(self, label, connection, placeholder):
        self.label, self.connection, self.ph = label, connection, placeholder

    def rows(self, sql, args=()):
        sql = sql.replace("%s", self.ph)
        if self.label == "sqlite":
            cur = self.connection.execute(sql, args)
            cols = [d[0] for d in cur.description]
            return [dict(zip(cols, r)) for r in cur.fetchall()]
        cur = self.connection.cursor(dictionary=True)
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

    def scalar(self, sql, args=()):
        row = self.rows(sql, args)[0]
        return list(row.values())[0]

    def has(self, table):
        if self.label == "sqlite":
            return bool(self.connection.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
                (table,)).fetchone())
        return bool(self.rows("SHOW TABLES LIKE %s", (table,)))


def read_merged(path: Path) -> dict[int, int]:
    out: dict[int, int] = {}
    with open(require_local(path), encoding="utf-8") as handle:
        for line in handle:
            parts = line.split("\t|")
            if len(parts) >= 2:
                out[int(parts[0])] = int(parts[1])
    return out


def evidence(store, taxon_id):
    abundance = store.scalar(
        "SELECT COUNT(*) FROM sample_taxon_abundances WHERE taxon_id = %s", (taxon_id,))
    assoc = 0
    if store.has("taxon_disease_associations"):
        assoc = store.scalar(
            "SELECT COUNT(*) FROM taxon_disease_associations WHERE taxon_id = %s",
            (taxon_id,))
    return int(abundance), int(assoc)


def classify(store, merged):
    """One decision per retired row, with the numbers it was decided on."""
    out = {}
    retired = store.rows(
        "SELECT * FROM taxa WHERE ncbi_tax_id IS NOT NULL AND ncbi_rank IS NULL")
    for row in retired:
        target = merged.get(int(row["ncbi_tax_id"]))
        if target is None:
            continue
        occupant = store.rows(
            "SELECT * FROM taxa WHERE ncbi_tax_id = %s AND taxonomic_rank = %s "
            "AND id <> %s", (target, row["taxonomic_rank"], row["id"]))
        if not occupant:
            continue
        occupant = occupant[0]
        ab_r, as_r = evidence(store, row["id"])
        ab_o, as_o = evidence(store, occupant["id"])
        shared = store.scalar(
            "SELECT COUNT(*) FROM (SELECT sample_id FROM sample_taxon_abundances "
            "WHERE taxon_id = %s INTERSECT SELECT sample_id FROM "
            "sample_taxon_abundances WHERE taxon_id = %s) both_rows",
            (row["id"], occupant["id"]))
        if ab_r == 0 and as_r == 0:
            case = DROP_RETIRED
        elif ab_o == 0 and as_o == 0:
            case = DROP_OCCUPANT
        elif shared == 0 and not (as_r and as_o):
            # Moving rows is only safe when nothing can collide on arrival:
            # no sample holds both (checked above), and associations cannot
            # both exist, since two association rows for one taxon and one
            # disease would conflict on arrival. Neither current pair has any.
            case = MOVE_ROWS
        else:
            case = REFUSE_SUM
        out[(row["genus_key"], row["species_key"])] = {
            "case": case, "row": row, "occupant": occupant, "target": target,
            "numbers": f"co-occur={shared} retired={ab_r}/{as_r} occupant={ab_o}/{as_o}",
        }
    return out


def carry(store, source, survivor, apply):
    missing = {c: source[c] for c in CARRIED
               if source.get(c) not in (None, "") and survivor.get(c) in (None, "")}
    if missing and apply:
        assignments = ", ".join(f"{c} = %s" for c in missing)
        store.execute(f"UPDATE taxa SET {assignments} WHERE id = %s",
                      (*missing.values(), survivor["id"]))
    return missing


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env-file", default=".env")
    parser.add_argument("--sqlite", default="gut_microbiome.sqlite")
    parser.add_argument("--taxdump", type=Path, default=Path("taxdump"))
    parser.add_argument("--audit-dir", type=Path, default=Path("data"))
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    merged = read_merged(args.taxdump / "merged.dmp")
    mysql = connect(Settings.from_env(args.env_file))
    lite = sqlite3.connect(args.sqlite, timeout=60)
    audit = []
    try:
        stores = [Store("mysql", mysql, "%s"), Store("sqlite", lite, "?")]
        decisions = {s.label: classify(s, merged) for s in stores}

        # An action runs only where every store holding the row agrees on it.
        agreed, disputed = {}, {}
        keys = set().union(*(set(d) for d in decisions.values()))
        for key in keys:
            cases = {label: d[key]["case"] for label, d in decisions.items() if key in d}
            if len(set(cases.values())) == 1:
                agreed[key] = next(iter(cases.values()))
            else:
                disputed[key] = cases

        for key, cases in sorted(disputed.items()):
            print(f"  DISPUTED {key}: stores classify it differently: {cases}")

        for store in stores:
            print(f"--- {store.label} ---")
            for key, plan in sorted(decisions[store.label].items()):
                case = agreed.get(key)
                name = plan["row"]["scientific_name"]
                other = plan["occupant"]["scientific_name"]
                if case is None:
                    continue
                if case == REFUSE_SUM:
                    print(f"  REFUSED {name!r} + {other!r}: {plan['numbers']} -- "
                          f"shares must be summed, which is a data decision")
                    continue
                print(f"  {case}: {name!r} + {other!r}  ({plan['numbers']})")
                if case == DROP_RETIRED:
                    survivor, doomed = plan["occupant"], plan["row"]
                elif case == DROP_OCCUPANT:
                    survivor, doomed = plan["row"], plan["occupant"]
                else:
                    survivor, doomed = plan["occupant"], plan["row"]
                # The survivor may not exist in every store -- the mirror holds
                # 54 keys MySQL does not, and lacks some MySQL rows. Folding
                # here then leaves the other store's copy on a retired ID, so
                # say so: repoint_merged_taxids.py closes it on a second pass,
                # since this store no longer holds the row to veto it. Refusing
                # outright would be worse, leaving a broken join behind.
                elsewhere = [other.label for other in stores
                             if other.label != store.label
                             and other.rows("SELECT 1 FROM taxa WHERE genus_key = %s "
                                            "AND species_key = %s", key)
                             and not other.rows(
                                 "SELECT 1 FROM taxa WHERE genus_key = %s AND "
                                 "species_key = %s",
                                 (survivor["genus_key"], survivor["species_key"]))]
                if elsewhere:
                    print(f"       NOTE {', '.join(elsewhere)} holds this row but "
                          f"not {survivor['scientific_name']!r}; re-run "
                          f"repoint_merged_taxids.py afterwards to converge it")
                moved = carry(store, doomed, survivor, args.apply)
                if moved:
                    print(f"       carrying onto {survivor['scientific_name']!r}: {moved}")
                if case == MOVE_ROWS and args.apply:
                    store.execute("UPDATE sample_taxon_abundances SET taxon_id = %s "
                                  "WHERE taxon_id = %s", (survivor["id"], doomed["id"]))
                    if store.has("taxon_disease_associations"):
                        store.execute("UPDATE taxon_disease_associations SET taxon_id = %s "
                                      "WHERE taxon_id = %s", (survivor["id"], doomed["id"]))
                if case == DROP_OCCUPANT and args.apply:
                    # Free the live ID first, then move the survivor onto it.
                    store.execute("DELETE FROM taxa WHERE id = %s", (doomed["id"],))
                    store.execute("UPDATE taxa SET ncbi_tax_id = %s WHERE id = %s",
                                  (plan["target"], survivor["id"]))
                elif args.apply:
                    store.execute("DELETE FROM taxa WHERE id = %s", (doomed["id"],))
                audit.append([store.label, case, doomed["id"], doomed["scientific_name"],
                              survivor["id"], survivor["scientific_name"],
                              plan["numbers"]])
            if args.apply:
                store.connection.commit()
            print()

        if args.apply:
            args.audit_dir.mkdir(parents=True, exist_ok=True)
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            path = args.audit_dir / f"retired_pair_fold_audit_{stamp}.csv"
            with path.open("w", newline="", encoding="utf-8") as fh:
                writer = csv.writer(fh)
                writer.writerow(["store", "case", "dropped_id", "dropped_name",
                                 "survivor_id", "survivor_name", "numbers"])
                writer.writerows(audit)
            print(f"Applied. Audit: {path.resolve()}")
        else:
            print("Preview only. Use --apply to write both stores.")
    finally:
        lite.close()
        mysql.close()


if __name__ == "__main__":
    main()
