#!/usr/bin/env python3
"""Correct trait values that contradict their genus, where that is a source slip.

gutdb/traits.py measures which traits a genus constrains and counts the rows
that dissent. This script resolves them, and it does not resolve them by
folding the minority onto the majority: that would be the mistake the phylum
work already names -- "an unapproved disagreement may be a curated value that is
deliberately better", so it is surfaced for a human rather than normalized away.

Three outcomes, and a row has to be put in one of them deliberately:

  * listed in FOLD below -- a confirmed source error, corrected;
  * listed in traits.APPROVED_DISSENTERS -- looked at, genuinely dissenting, so
    it never reaches this script and never trips the check;
  * neither -- reported only, and the check stays non-zero until someone picks.

The third case is the one that matters. A check that forced a repair would have
quietly rewritten Clostridium symbiosum and Clostridium innocuum, which are
gram-positive organisms that stain gram-negative, on 4,143 and 4,336 abundance
rows between them.

Keyed on (genus_key, species_key), never on taxa.id: the mirror numbers its rows
independently, so Aspergillus nidulans is 725 in MySQL and 724 there.

Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.traits import (APPROVED_DISSENTERS, GENUS_CONSTRAINED_TRAITS,
                          NON_ANSWERS, dissents_from_genus,
                          load_genus_consensus)

DRY_RUN = "--apply" not in sys.argv
SQLITE_PATH = "gut_microbiome.sqlite"

# (trait, genus_key, species_key) -> why this one is a source error.
FOLD: dict[tuple[str, str, str], str] = {
    # MiMeDB is internally inconsistent here rather than describing an unusual
    # organism: 28 of its 29 Aspergillus records say positive and this one says
    # negative, and across its fungi the split is 141 positive to 18 negative.
    # Fungal cell walls are conventionally reported gram-positive, A. nidulans
    # is not a documented exception to its own genus, and the row carries no
    # abundances or associations, so nothing downstream moves.
    ("gram_stain", "aspergillus", "nidulans"):
        "MiMeDB inconsistency: 28 of 29 Aspergillus records say positive",
}


def resolve(store: str, cursor) -> tuple[int, int]:
    ph = "%s" if store == "mysql" else "?"
    consensus = load_genus_consensus(cursor, ph)
    folded = reported = 0

    for trait in GENUS_CONSTRAINED_TRAITS:
        cursor.execute(
            f"SELECT id, scientific_name, genus_key, species_key, {trait} "
            f"FROM taxa WHERE genus_key IS NOT NULL AND TRIM(genus_key) <> '' "
            f"  AND {trait} IS NOT NULL"
        )
        for taxon_id, name, genus_key, species_key, stored in cursor.fetchall():
            agreed = dissents_from_genus(consensus, genus_key, trait, stored)
            if agreed is None:
                continue
            key = (trait, (genus_key or "").lower(), (species_key or "").lower())
            if key in APPROVED_DISSENTERS:
                # Excluded from the check too; printed so the list stays visible.
                print(f"  {trait}: {name} -- dissent approved, left alone")
                print(f"      {APPROVED_DISSENTERS[key]}")
                continue
            current = str(stored).strip().lower()
            if key in FOLD:
                print(f"  {trait}: {name} (id {taxon_id}) "
                      f"{current!r} -> {agreed!r}")
                print(f"      {FOLD[key]}")
                if not DRY_RUN:
                    cursor.execute(
                        f"UPDATE taxa SET {trait} = {ph} WHERE id = {ph}",
                        (agreed, taxon_id))
                folded += 1
                continue
            print(f"  {trait}: {name} (id {taxon_id}) dissents "
                  f"({current!r} against {agreed!r}) -- REPORTED ONLY")
            print("      Decide it: add to FOLD if the source is wrong, or to "
                  "traits.APPROVED_DISSENTERS if the row is right.")
            reported += 1
    return folded, reported


def main() -> int:
    totals = [0, 0]

    print("mysql")
    conn = connect(Settings.from_env(".env"))
    cursor = conn.cursor()
    folded, reported = resolve("mysql", cursor)
    totals[0] += folded
    totals[1] += reported
    if not DRY_RUN:
        conn.commit()
    cursor.close()
    conn.close()

    print("\nsqlite")
    if Path(SQLITE_PATH).exists():
        lite = sqlite3.connect(SQLITE_PATH)
        cursor = lite.cursor()
        folded, reported = resolve("sqlite", cursor)
        totals[0] += folded
        totals[1] += reported
        if not DRY_RUN:
            lite.commit()
        cursor.close()
        lite.close()
    else:
        print(f"  (no SQLite mirror at {SQLITE_PATH}; skipped)")

    print(f"\n  corrected      : {totals[0]}")
    print(f"  reported only  : {totals[1]}")
    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    # A reported-only row leaves the matching validate counter non-zero, which
    # the pre-push gate refuses, so exit non-zero rather than looking clean.
    return 1 if totals[1] else 0


if __name__ == "__main__":
    raise SystemExit(main())
