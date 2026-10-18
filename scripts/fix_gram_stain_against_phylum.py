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

import os
import sqlite3
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
    cur.close()
    conn.close()

    m_written, m_skipped = repair_mirror()
    written += m_written
    skipped += m_skipped

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    return 1 if skipped else 0


def repair_mirror() -> tuple[int, int]:
    """Apply the same corrections to the SQLite mirror, keyed on NAME not id.

    The mirror held six of these seven when the check was finally able to run
    against it -- Parabacteroides distasonis among them, with 9,291 abundance
    rows there -- because this repair reached MySQL only. That is the defect
    shape this project keeps producing, and the reason validate now runs its
    portable checks against both stores.

    KEYED ON scientific_name, NOT ON THE IDS ABOVE. The mirror has its own
    taxa.id space: Vibrio fischeri is 1907 in MySQL and 1904 here, so reusing
    those ids would silently correct whatever organism happens to hold them. The
    name is what identifies the row whose gram character the list settles, and
    the phylum and wrong-value guards still have to pass.
    """
    path = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")
    if not Path(path).exists():
        print(f"\n  sqlite mirror {path} not present; skipped")
        return 0, 0

    print("\nsqlite mirror:")
    connection = sqlite3.connect(path)
    cursor = connection.cursor()
    written = skipped = 0
    for _mysql_id, name, phylum, wrong, correct, why in REPAIRS:
        cursor.execute(
            "SELECT id, phylum, gram_stain FROM taxa WHERE scientific_name = ?",
            (name,))
        rows = cursor.fetchall()
        if not rows:
            print(f"  {name}: not in the mirror")
            continue
        if len(rows) > 1:
            print(f"  {name}: SKIPPED, {len(rows)} rows share this name here")
            skipped += 1
            continue
        row_id, stored_phylum, stored_gram = rows[0]
        current = (stored_gram or "").strip().lower()
        if current == correct:
            print(f"  {name}: already {correct!r}")
            continue
        if (stored_phylum or "").strip() != phylum:
            print(f"  {name}: SKIPPED, expected phylum {phylum!r}, "
                  f"found {stored_phylum!r}")
            skipped += 1
            continue
        if not current:
            # The mirror never carried this trait for this row. Nothing to
            # correct and nothing wrong, so it is not counted as a skip -- the
            # exit code is reserved for a row that needs a human.
            print(f"  {name}: no gram stain recorded here, nothing to correct")
            continue
        if current != wrong:
            print(f"  {name}: SKIPPED, expected {wrong!r}, found {current!r}")
            skipped += 1
            continue
        print(f"  {name} ({phylum}) id {row_id}: {current!r} -> {correct!r}")
        if not DRY_RUN:
            cursor.execute("UPDATE taxa SET gram_stain = ? WHERE id = ?",
                           (correct, row_id))
        written += 1
    if not DRY_RUN:
        connection.commit()
    connection.close()
    print(f"  corrected : {written}")
    print(f"  skipped   : {skipped}")
    return written, skipped


if __name__ == "__main__":
    raise SystemExit(main())
