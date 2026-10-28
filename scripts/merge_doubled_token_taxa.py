#!/usr/bin/env python3
"""Fold a taxa row whose name says one word twice into its well-formed twin.

Three rows read like `Unclassified burkholderiales burkholderiales bacterium`,
a name no organism bears. The cause is not a broken field: `genus` is
`Unclassified Burkholderiales` and `species` is `Burkholderiales bacterium`,
which is NCBI's name for an unnamed member of that order, and
tests/test_transform.py deliberately pins that epithet rather than stripping
its leading token -- a rank name is not a genus. The doubling comes from
`upsert_taxon` joining the two with f"{genus} {species}".

Fixing the join alone would make things worse. `_taxon_exists` keys on
species_key, so the well-formed row for the same organism already exists with
species_key 'bacterium' (ids 2089, 2091, 1863, and those are the only ones the
SQLite mirror carries); renaming without converging the key would leave two
rows wearing one name. So this is a duplicate pair to fold, not a string to
repair.

IT IS A MERGE, NOT A DELETE, AND THAT DISTINCTION IS THE POINT. The malformed
row holds what the survivor lacks: the only ncbi_tax_id those three organisms
have (469610, 457421, 658089) plus gram_stain and metabolism. Deleting it
outright would destroy that, which is the failure mode this repo has already
recorded once in a repair script whose premise had expired. Every column the
survivor is missing is copied first, into BOTH stores, because the survivor
exists in both and leaving the mirror behind is this project's recurring defect.

A row is folded only when it carries no abundance and no associations. If it
carries evidence it is not a stale duplicate, and this refuses rather than
choosing where that evidence should go.

Preview by default; --apply writes an audit and the changes.
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

CARRIED = [
    "ncbi_tax_id", "ncbi_rank", "gram_stain", "oxygen_requirement",
    "ph_preference", "sporulation", "shape", "cell_arrangement", "mobility",
    "flagella_presence", "number_of_membranes", "biotic_relationship",
    "habitat", "temperature_range", "optimal_temperature", "metabolism",
    "energy_source", "human_pathogen",
]
EVIDENCE = [("sample_taxon_abundances", "taxon_id"),
            ("taxon_disease_associations", "taxon_id")]


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

    def table_exists(self, name):
        if self.label == "sqlite":
            return bool(self.connection.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
                (name,)).fetchone())
        return bool(self.rows("SHOW TABLES LIKE %s", (name,)))


def doubled(store):
    """Rows whose species field repeats the final word of a multi-word genus."""
    out = []
    for row in store.rows("SELECT * FROM taxa WHERE genus IS NOT NULL "
                          "AND species IS NOT NULL AND species <> ''"):
        genus_words = (row["genus"] or "").split()
        species_words = (row["species"] or "").split()
        if len(genus_words) < 2 or len(species_words) < 2:
            continue
        if species_words[0].casefold() == genus_words[-1].casefold():
            out.append((row, " ".join(species_words[1:])))
    return out


def evidence_for(store, taxon_id):
    counts = {}
    for table, column in EVIDENCE:
        if store.table_exists(table):
            counts[table] = store.rows(
                f"SELECT COUNT(*) AS n FROM {table} WHERE {column} = %s",
                (taxon_id,))[0]["n"]
    return {k: v for k, v in counts.items() if v}


def plan(store):
    folds, refusals = [], []
    for row, reduced in doubled(store):
        twin = store.rows("SELECT * FROM taxa WHERE genus_key = %s AND species_key = %s",
                          (row["genus_key"], reduced.casefold()))
        if not twin:
            refusals.append(f"{row['scientific_name']!r}: no row carries the "
                            f"reduced key ({row['genus_key']!r}, {reduced.casefold()!r})")
            continue
        held = evidence_for(store, row["id"])
        if held:
            refusals.append(f"{row['scientific_name']!r}: carries {held}, so it is "
                            f"not a stale duplicate")
            continue
        carry = {c: row[c] for c in CARRIED
                 if row.get(c) not in (None, "") and twin[0].get(c) in (None, "")}
        folds.append({"row": row, "twin": twin[0], "carry": carry})
    return folds, refusals


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
        stores = [Store("mysql", mysql, "%s"), Store("sqlite", lite, "?")]
        # The survivor must be enriched in every store that has it, even where
        # the malformed row never reached -- the mirror holds only the twins.
        carried_by_key = {}
        for store in stores:
            folds, refusals = plan(store)
            print(f"--- {store.label} ---")
            for fold in folds:
                key = (fold["twin"]["genus_key"], fold["twin"]["species_key"])
                carried_by_key.setdefault(key, {}).update(fold["carry"])
                print(f"  fold {fold['row']['scientific_name']!r} (id "
                      f"{fold['row']['id']}) into {fold['twin']['scientific_name']!r} "
                      f"(id {fold['twin']['id']})")
                print(f"       carrying {fold['carry'] or 'nothing the twin lacks'}")
                audit.append([store.label, "fold", fold["row"]["id"],
                              fold["row"]["scientific_name"], fold["twin"]["id"],
                              fold["twin"]["scientific_name"],
                              ";".join(f"{k}={v}" for k, v in fold["carry"].items())])
            for message in refusals:
                print(f"  REFUSED {message}")
            if not folds and not refusals:
                print("  nothing to do")

        for store in stores:
            for key, carry in carried_by_key.items():
                twin = store.rows("SELECT * FROM taxa WHERE genus_key = %s "
                                  "AND species_key = %s", key)
                if not twin:
                    continue
                missing = {c: v for c, v in carry.items() if twin[0].get(c) in (None, "")}
                if missing:
                    print(f"  {store.label}: enrich {twin[0]['scientific_name']!r} "
                          f"(id {twin[0]['id']}) with {missing}")
                    audit.append([store.label, "enrich", "", "", twin[0]["id"],
                                  twin[0]["scientific_name"],
                                  ";".join(f"{k}={v}" for k, v in missing.items())])
                    if args.apply:
                        assignments = ", ".join(f"{c} = %s" for c in missing)
                        store.execute(f"UPDATE taxa SET {assignments} WHERE id = %s",
                                      (*missing.values(), twin[0]["id"]))

        if args.apply:
            for store in stores:
                folds, _ = plan(store)
                for fold in folds:
                    store.execute("DELETE FROM taxa WHERE id = %s", (fold["row"]["id"],))
                store.connection.commit()
            args.audit_dir.mkdir(parents=True, exist_ok=True)
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            path = args.audit_dir / f"doubled_token_merge_audit_{stamp}.csv"
            with path.open("w", newline="", encoding="utf-8") as fh:
                writer = csv.writer(fh)
                writer.writerow(["store", "action", "folded_id", "folded_name",
                                 "survivor_id", "survivor_name", "carried"])
                writer.writerows(audit)
            print(f"\nApplied. Audit: {path.resolve()}")
        else:
            print("\nPreview only. Use --apply to write both stores.")
    finally:
        lite.close()
        mysql.close()


if __name__ == "__main__":
    main()
