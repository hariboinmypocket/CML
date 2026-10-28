#!/usr/bin/env python3
"""Restore the underscores NCBI puts in strain designators, in MySQL.

The two stores key the same organism two different ways. MySQL holds
`Bacteroides sp. 1 1 14` and the SQLite mirror holds `Bacteroides sp. 1_1_14`,
which is what NCBI calls tax ID 469585 -- so the mirror and NCBI agree and
MySQL is the outlier. 24 rows are affected, and the difference is case as well
as punctuation: NCBI writes `3_1_33FAA` and `2_1_33B` where MySQL has
`3 1 33faa` and `2 1 33b`.

THE COST IS NOT COSMETIC. (genus_key, species_key) is how every two-store repair
in this repo matches a row, because the mirror keeps its own taxa.id space. A row
keyed differently in each store is therefore invisible to all of them, and to
tests/test_cross_store_tax_ids.py, which compares only keys present in both: 24
of the 54 keys the mirror holds and MySQL does not are these variants, so those
rows have never been compared across the stores at all.

No transform is to blame. clean_str collapses whitespace and touches no
underscore, so the spelling arrived this way from the source, which means a
reload can bring it back. That is why the invariant is asserted in the test
rather than assumed to stay fixed.

Each row is checked individually against names.dmp and changed only where the
stored name differs from NCBI's scientific name for its own tax ID by nothing
but spaces, underscores and case -- never by a blanket space-to-underscore
substitution, which would corrupt the many legitimate spaces in names like
`sp. oral taxon 102`. The name, the species field and species_key are set
together from NCBI's spelling so they cannot drift apart, and a row is skipped
if the corrected key is already taken.

Preview by default; --apply writes MySQL and an audit.
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
from gutdb.transform import parse_scientific_name


def squash(text: str) -> str:
    """Spaces, underscores and case made irrelevant -- nothing else."""
    return re.sub(r"[ _]+", " ", text.strip().lower())


def ncbi_names(path: Path, wanted: set[int]) -> dict[int, str]:
    out: dict[int, str] = {}
    with open(require_local(path), encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split("\t|")
            if len(parts) >= 4 and parts[3].strip() == "scientific name":
                tax_id = int(parts[0])
                if tax_id in wanted:
                    out[tax_id] = parts[1].strip()
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env-file", default=".env")
    parser.add_argument("--taxdump", type=Path, default=Path("taxdump"))
    parser.add_argument("--sqlite", default="gut_microbiome.sqlite")
    parser.add_argument("--audit-dir", type=Path, default=Path("data"))
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    connection = connect(Settings.from_env(args.env_file))
    cursor = connection.cursor()
    try:
        # Every row, not just those with a tax ID: the rows that need the
        # mirror fallback are precisely the ones with none.
        cursor.execute("SELECT id, scientific_name, ncbi_tax_id, genus, species, "
                       "genus_key, species_key FROM taxa")
        rows = cursor.fetchall()
        names = ncbi_names(args.taxdump / "names.dmp",
                           {int(r[2]) for r in rows if r[2] is not None})
        taken = {(r[5], r[6]) for r in rows}

        fixes, refusals = [], []
        for row_id, stored, tax_id, genus, species, genus_key, species_key in rows:
            if tax_id is None:
                continue
            ncbi = names.get(int(tax_id))
            if not ncbi or stored.strip() == ncbi:
                continue
            if squash(stored) != squash(ncbi) or "_" not in ncbi:
                continue
            # NCBI's punctuation, this project's capitalisation. parse_
            # scientific_name lowercases the epithet, and scientific_name is
            # genus + species everywhere else -- 4 rows out of 5,829 break that
            # relationship today. Storing NCBI's `3_1_33FAA` while species held
            # `3_1_33faa` would take that to 42 for a cosmetic gain, so only the
            # underscores are adopted. The rows that differ from NCBI by case
            # alone are therefore already correct, and drop out here.
            new_genus, new_species = parse_scientific_name(ncbi)
            rendered = f"{new_genus} {new_species}".strip()
            if rendered == stored.strip():
                continue
            new_key = (new_genus.casefold(), new_species.casefold())
            if new_key != (genus_key, species_key) and new_key in taken:
                refusals.append(f"{stored!r}: corrected key {new_key} is already "
                                f"held by another row")
                continue
            fixes.append((row_id, stored, rendered, new_genus, new_species, new_key))

        # A row with no tax ID has nothing in the dump to check against, and
        # three such rows are still keyed one way here and another in the
        # mirror. The mirror's spelling is a weaker basis than NCBI's -- it is
        # the other store, not the authority -- so it is used only where the
        # mirror's key is this key with underscores restored, the row carries no
        # evidence, and no NCBI node exists to ask. The basis is recorded per
        # row in the audit rather than blurred into the NCBI-confirmed ones.
        lite = sqlite3.connect(args.sqlite, timeout=60)
        try:
            mirror = {(a, b): n for a, b, n in lite.execute(
                "SELECT genus_key, species_key, scientific_name FROM taxa").fetchall()}
        finally:
            lite.close()
        basis = {f[0]: "ncbi" for f in fixes}
        for row_id, stored, tax_id, genus, species, genus_key, species_key in rows:
            if tax_id is not None or any(f[0] == row_id for f in fixes):
                continue
            candidates = [k for k in mirror
                          if k[0] == genus_key and k not in taken
                          and squash(k[1]) == squash(species_key) and "_" in k[1]]
            if len(candidates) != 1:
                continue
            key = candidates[0]
            new_species = key[1]
            rendered = f"{genus} {new_species}".strip()
            if rendered == stored.strip():
                continue
            fixes.append((row_id, stored, rendered, genus, new_species, key))
            basis[row_id] = "mirror spelling; no NCBI node to confirm"

        print(f"rows to respell: {len(fixes)}, refused: {len(refusals)}\n")
        for row_id, stored, corrected, new_genus, new_species, new_key in fixes:
            print(f"  id={row_id}: {stored!r} -> {corrected!r}  [{basis[row_id]}]")
        for message in refusals:
            print(f"  REFUSED {message}")

        if args.apply and fixes:
            for row_id, stored, corrected, new_genus, new_species, new_key in fixes:
                cursor.execute(
                    "UPDATE taxa SET scientific_name = %s, genus = %s, species = %s, "
                    "genus_key = %s, species_key = %s WHERE id = %s",
                    (corrected, new_genus, new_species, new_key[0], new_key[1], row_id))
            connection.commit()
            args.audit_dir.mkdir(parents=True, exist_ok=True)
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            path = args.audit_dir / f"strain_underscore_audit_{stamp}.csv"
            with path.open("w", newline="", encoding="utf-8") as fh:
                writer = csv.writer(fh)
                writer.writerow(["taxon_id", "was", "now", "genus", "species", "basis"])
                writer.writerows([[f[0], f[1], f[2], f[3], f[4], basis[f[0]]]
                                  for f in fixes])
            print(f"\nApplied to MySQL. Audit: {path.resolve()}")
        elif not args.apply:
            print("\nPreview only. Use --apply to write MySQL.")
    finally:
        cursor.close()
        connection.close()


if __name__ == "__main__":
    main()
