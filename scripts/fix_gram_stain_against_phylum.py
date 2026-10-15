#!/usr/bin/env python3
"""Correct seven gram_stain values that contradict their phylum.

Gram stain describes the cell envelope, so it is near-constant within a phylum.
Seven rows dissented, and every one of them is a taxon whose gram character is
textbook: Bacteroidetes and the Vibrio/Yersinia/Yokenella enterobacteria are
gram-negative. They were not a discovery, they were a bad source value -- MiMeDB
records 'positive' for all of them, 8 source rows across the five
Proteobacteria plus MMDBm0000056/57 for Bacteroides vulgatus.

Two of the seven matter. The MiMeDB trait backfill moved that Bacteroides
vulgatus value onto Phocaeicola vulgatus (row 1569, 10,263 abundance rows, 30
associations), and Parabacteroides distasonis (row 216, 9,731 abundances) holds
the same error from its own load. They were the only two positives among 109
Bacteroidetes. The other five carry 0-2 abundance rows between them; they are
corrected because leaving a known-wrong value to surface later when the row gains
evidence is how the first two happened.

Fixed as an explicit list rather than by folding the minority onto the phylum
majority. The majority is the right way to *detect* this -- validate's
gram_stain_against_phylum check does exactly that -- but it is the wrong way to
decide it, because a phylum can genuinely hold both: Firmicutes contains the
gram-negative Negativicutes, and a blanket fold would rewrite them. Each row
below is named with the organism that settles it.

The source CSV is not corrected, so a MiMeDB reload reintroduces the five rows
sourced from it. validate will catch that; this script re-run will fix it.

Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect

DRY_RUN = "--apply" not in sys.argv

# (taxon id, expected name, expected phylum, wrong value, correct value, why)
REPAIRS = (
    (1569, "Phocaeicola vulgatus", "Bacteroidetes", "positive", "negative",
     "Bacteroidetes are gram-negative; value arrived from MiMeDB's "
     "Bacteroides vulgatus records via the trait backfill"),
    (216, "Parabacteroides distasonis", "Bacteroidetes", "positive", "negative",
     "Bacteroidetes are gram-negative; error predates and is independent of "
     "the concatenated-genus repair"),
    (1907, "Vibrio fischeri", "Proteobacteria", "positive", "negative",
     "Vibrio are gram-negative"),
    (1911, "Vibrio furnissii", "Proteobacteria", "positive", "negative",
     "Vibrio are gram-negative"),
    (1914, "Vibrio mimicus", "Proteobacteria", "positive", "negative",
     "Vibrio are gram-negative"),
    (1998, "Yersinia rohdei", "Proteobacteria", "positive", "negative",
     "Yersinia are gram-negative"),
    (2000, "Yokenella regensburgei", "Proteobacteria", "positive", "negative",
     "Yokenella is an enterobacterium, gram-negative"),
)


def main() -> int:
    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()

    written = skipped = 0
    for taxon_id, name, phylum, wrong, correct, why in REPAIRS:
        cur.execute(
            "SELECT scientific_name, phylum, gram_stain FROM taxa WHERE id = %s",
            (taxon_id,))
        row = cur.fetchone()
        if not row:
            print(f"  id {taxon_id}: SKIPPED, no such row")
            skipped += 1
            continue
        stored_name, stored_phylum, stored_gram = row
        current = (stored_gram or "").strip().lower()

        if current == correct:
            print(f"  id {taxon_id}: already {correct!r} ({stored_name})")
            continue
        # Guard on the row still being the one described above. The database
        # changes hourly, and a row renamed or reclassified since this list was
        # written is no longer the organism whose gram character settles it.
        if stored_name != name:
            print(f"  id {taxon_id}: SKIPPED, expected {name!r}, found {stored_name!r}")
            skipped += 1
            continue
        if (stored_phylum or "").strip() != phylum:
            print(f"  id {taxon_id}: SKIPPED, expected phylum {phylum!r}, "
                  f"found {stored_phylum!r}")
            skipped += 1
            continue
        if current != wrong:
            print(f"  id {taxon_id}: SKIPPED, expected {wrong!r}, found {current!r}")
            skipped += 1
            continue

        print(f"  id {taxon_id}: {stored_name} ({phylum}) "
              f"{current!r} -> {correct!r}")
        print(f"      {why}")
        if not DRY_RUN:
            cur.execute("UPDATE taxa SET gram_stain = %s WHERE id = %s",
                        (correct, taxon_id))
        written += 1

    if not DRY_RUN:
        conn.commit()

    print(f"\n  corrected : {written}")
    print(f"  skipped   : {skipped}")
    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    cur.close()
    conn.close()
    return 1 if skipped else 0


if __name__ == "__main__":
    raise SystemExit(main())
