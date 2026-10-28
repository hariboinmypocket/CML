#!/usr/bin/env python3
"""Correct tax IDs that name a different organism than the row they sit on.

Found by taxa_sharing_a_tax_id_within_one_rank, which reported 0 in MySQL and 2
in the SQLite mirror. The mirror was not the broken one -- it was the only store
whose data made the collision visible.

VERIFIED AGAINST taxdump/names.dmp AND nodes.dmp, not against intuition:

  29494  is `Vibrio furnissii`, rank species.
  903510 is `Vibrio furnissii NCTC 11218`, rank STRAIN.
  668    is `Aliivibrio fischeri`, rank species, and lists `Vibrio fischeri`
         among its synonyms -- so it is the same organism under an older name,
         which is the tax ID identity test this project uses: same ID means a
         rename, a different ID would mean a reclassification.

So `Vibrio fischeri` was carrying furnissii's species ID in BOTH stores, and
MySQL's `Vibrio furnissii` was pointing at a strain rather than at the species.
MySQL avoided colliding only because its two rows were wrong in two different
ways; the mirror had furnissii right, so its fischeri row collided with it and
the counter fired.

The third correction is a rename left half-applied. Tax ID 1624 is
`Ligilactobacillus salivarius`, the 2021 ICNP name for `Lactobacillus
salivarius`. The source CSV carries both spellings as separate rows
(`Lactobacillus,Salivarius` and `Ligilactobacillus,salivarius`); MySQL folded
them and the mirror did not, so the mirror holds a dead old-name row sharing the
tax ID with the live one.

Rows are matched on (genus_key, species_key), never on id, because the mirror
maintains its own taxa.id space and an id-keyed repair silently hits a different
organism there. Evidence is checked before anything is removed: a row is only
deleted if it carries no abundance and no associations, and the script refuses
rather than choosing if it does.

Preview by default; --apply writes an audit and the corrections.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path
import sqlite3
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect


def _ranks_from_dump():
    """Reuse the canonical nodes.dmp reader rather than parse it again here."""
    import importlib.util
    path = Path(__file__).resolve().parent / "backfill_ncbi_rank.py"
    spec = importlib.util.spec_from_file_location("backfill_ncbi_rank", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.ranks_from_dump()

# (genus_key, species_key, wrong_tax_id, correct_tax_id, why)
RETAG = [
    ("vibrio", "fischeri", 29494, 668,
     "29494 is Vibrio furnissii; 668 is Aliivibrio fischeri, which lists "
     "Vibrio fischeri as a synonym"),
    ("vibrio", "furnissii", 903510, 29494,
     "903510 is the strain Vibrio furnissii NCTC 11218; the species is 29494"),
    ("oxalobacter", "formigenes", 556268, 847,
     "556268 is Oxalobacter paraformigenes (Chmiel et al. 2022), a species "
     "split out of formigenes that lists only the strain Oxalobacter "
     "formigenes HOxBLS under `includes`; the species is 847"),
    ("klebsiella", "pneumoniae", 1328388, 573,
     "1328388 is absent from nodes.dmp -- NCBI deleted it -- so 2,681 "
     "abundance rows joined to nothing; 573 is Klebsiella pneumoniae"),
    ("eubacterium", "saphenum", 592031, 51123,
     "592031 is the strain Eubacterium saphenum ATCC 49989; the species is "
     "51123. The genus-layer row that already carries 51123 is a different "
     "layer, so the within-rank collision counter is unaffected"),
]
# (genus_key, species_key, tax_id, why) -- removed only when evidence-free.
DROP_STALE_SYNONYM = [
    ("lactobacillus", "salivarius", 1624,
     "dead pre-2021 name sharing tax ID 1624 with Ligilactobacillus "
     "salivarius, which carries the evidence"),
]


def rows_for(store, genus_key, species_key):
    sql = ("SELECT id, scientific_name, ncbi_tax_id FROM taxa "
           "WHERE genus_key = %s AND species_key = %s")
    return store.query(sql, (genus_key, species_key))


class MySQLStore:
    label = "mysql"
    placeholder = "%s"

    def __init__(self, connection):
        self.connection = connection

    def query(self, sql, args=()):
        cursor = self.connection.cursor()
        try:
            cursor.execute(sql, args)
            return cursor.fetchall()
        finally:
            cursor.close()

    def execute(self, sql, args=()):
        cursor = self.connection.cursor()
        try:
            cursor.execute(sql, args)
        finally:
            cursor.close()

    def commit(self):
        self.connection.commit()


class SqliteStore:
    label = "sqlite"
    placeholder = "?"

    def __init__(self, connection):
        self.connection = connection

    def query(self, sql, args=()):
        return self.connection.execute(sql.replace("%s", "?"), args).fetchall()

    def execute(self, sql, args=()):
        self.connection.execute(sql.replace("%s", "?"), args)

    def commit(self):
        self.connection.commit()


def evidence(store, taxon_id) -> tuple[int, int]:
    abundance = store.query(
        "SELECT COUNT(*) FROM sample_taxon_abundances WHERE taxon_id = %s", (taxon_id,))[0][0]
    associations = store.query(
        "SELECT COUNT(*) FROM taxon_disease_associations WHERE taxon_id = %s", (taxon_id,))[0][0]
    return int(abundance), int(associations)


def plan_for(store, ranks):
    retag, drop, refused, reranked = [], [], [], []
    for genus_key, species_key, wrong, correct, why in RETAG:
        for taxon_id, name, current in rows_for(store, genus_key, species_key):
            if current == wrong:
                retag.append((taxon_id, name, wrong, correct, why))
            elif current != correct:
                refused.append(f"{name}: tax ID is {current}, expected {wrong} or {correct}")
            # ncbi_rank describes whatever id the row carried, so a retag leaves
            # it stale -- MySQL kept 'strain' after Vibrio furnissii moved off
            # the strain id. Checked even when the id needed no change, so
            # re-running converges instead of depending on call order.
            stored = store.query(
                "SELECT ncbi_rank FROM taxa WHERE id = %s", (taxon_id,))[0][0]
            expected = ranks.get(correct)
            if expected and stored != expected:
                reranked.append((taxon_id, name, stored, expected, correct))
    for genus_key, species_key, tax_id, why in DROP_STALE_SYNONYM:
        for taxon_id, name, current in rows_for(store, genus_key, species_key):
            if current != tax_id:
                continue
            ab, assoc = evidence(store, taxon_id)
            if ab or assoc:
                refused.append(f"{name}: carries {ab} abundance rows and {assoc} "
                               f"associations, so it is not a dead synonym")
            else:
                drop.append((taxon_id, name, tax_id, why))
    return retag, drop, refused, reranked


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env-file", default=".env")
    parser.add_argument("--sqlite", default="gut_microbiome.sqlite")
    parser.add_argument("--audit-dir", type=Path, default=Path("data"))
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    ranks = _ranks_from_dump()
    print(f"ranks read from the dump: {len(ranks):,}\n")
    mysql = connect(Settings.from_env(args.env_file))
    lite = sqlite3.connect(args.sqlite, timeout=60)
    audit_rows = []
    try:
        for store in (MySQLStore(mysql), SqliteStore(lite)):
            retag, drop, refused, reranked = plan_for(store, ranks)
            print(f"--- {store.label} ---")
            for taxon_id, name, wrong, correct, why in retag:
                print(f"  retag {name!r}: {wrong} -> {correct}  ({why})")
                audit_rows.append([store.label, "retag", taxon_id, name, wrong, correct, why])
            for taxon_id, name, tax_id, why in drop:
                print(f"  drop  {name!r} (tax ID {tax_id}, no evidence)  ({why})")
                audit_rows.append([store.label, "drop", taxon_id, name, tax_id, "", why])
            for taxon_id, name, stored, expected, tax_id in reranked:
                print(f"  rerank {name!r}: ncbi_rank {stored!r} -> {expected!r} "
                      f"(what taxdump gives {tax_id})")
                audit_rows.append([store.label, "rerank", taxon_id, name,
                                   stored, expected, f"taxdump rank of {tax_id}"])
            for message in refused:
                print(f"  REFUSED {message}")
            if not retag and not drop and not reranked:
                print("  nothing to do")
            if not args.apply:
                continue
            for taxon_id, name, wrong, correct, why in retag:
                store.execute("UPDATE taxa SET ncbi_tax_id = %s WHERE id = %s",
                              (correct, taxon_id))
            for taxon_id, name, tax_id, why in drop:
                store.execute("DELETE FROM taxa WHERE id = %s", (taxon_id,))
            for taxon_id, name, stored, expected, tax_id in reranked:
                store.execute("UPDATE taxa SET ncbi_rank = %s WHERE id = %s",
                              (expected, taxon_id))
            store.commit()
        if not args.apply:
            print("\nPreview only. Use --apply to correct both stores.")
            return
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        audit = args.audit_dir / f"tax_id_correction_audit_{stamp}.csv"
        with audit.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["store", "action", "taxon_id", "scientific_name",
                             "tax_id_before", "tax_id_after", "reason"])
            writer.writerows(audit_rows)
        print(f"\nApplied. Audit: {audit.resolve()}")
    finally:
        mysql.close()
        lite.close()


if __name__ == "__main__":
    main()
