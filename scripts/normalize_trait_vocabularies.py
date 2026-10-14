#!/usr/bin/env python3
"""Apply the taxa trait vocabularies to MySQL, which never had them.

scripts/cleanup_gutdb.py has normalized these columns since it was written, but
it runs only against the SQLite mirror. So the two stores held different
vocabularies for the same rows, and this is the third time that has happened --
the (SeqCode) annotations and body_site were the first two.

    column               MySQL distinct   mirror distinct
    oxygen_requirement               21                10
    gram_stain                        7                 4
    temperature_range                 6                 5
    sporulation                      32                 2

oxygen_requirement is the one that costs something. 'anaerobe' (576),
'Anaerobic' (50) and 'Obligate anaerobe' (35) were three groups rather than one,
and 'aerobe' (486), 'Aerobic' (161), 'obligate aerobe' (8) and 'Obligate
aerobic' (1) were four. Oxygen tolerance is among the few mechanistic traits for
a gut community -- the lumen is anaerobic, and a shift toward aerotolerant taxa
is a recognised dysbiosis signature -- so splitting it across case variants
quietly weakens any model or query that groups on it.

sporulation was worse than inconsistent: it held three different kinds of value
at once.

    0 / 1                                     429   a yes/no flag as integers
    Sporulating / Nonsporulating              183   the same flag as words
    Sporulating (ascospores); reproduction...  216   prose about spore type
    Singles / Pairs - Singles / Chains ...     128   CELL ARRANGEMENT, another axis

The last group is the real defect: that text answers how cells sit next to each
other, which `cell_arrangement` already exists for. Those values are moved there
rather than deleted, and only when cell_arrangement is empty, so a real reading
is never overwritten by a value recovered from the wrong column. In practice only
2 of the 128 need the copy -- for the other 126 cell_arrangement already holds
the same reading, which means nothing was ever stored only in the wrong column
and clearing it loses nothing.

Nothing unrecognised is touched. A value no rule matches keeps its cell, because
the alternative is discarding a category nobody anticipated.

Vocabularies live in gutdb/transform.py so both stores use one implementation.
Idempotent. Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.transform import (FREE_TEXT_TRAITS, TRAIT_VOCABULARIES, canon_trait,
                             split_sporulation)

DRY_RUN = "--apply" not in sys.argv
AUDIT = Path("data/trait_vocabulary_audit.csv")
AUDIT_FIELDS = ("column", "stored", "normalized", "n_taxa", "action")


def normalize_column(cur, column, mapping, rows_out) -> int:
    """Fold one column onto its vocabulary. Returns taxa rewritten."""
    cur.execute(f"SELECT {column}, COUNT(*) FROM taxa "
                f"WHERE {column} IS NOT NULL AND TRIM({column}) <> '' "
                f"GROUP BY {column}")
    changed = 0
    for stored, n in cur.fetchall():
        tidy = canon_trait(stored, mapping)
        # Compare cleaned to cleaned: an integer-looking cell reads back
        # differently from its stored text and would look changed every run.
        if tidy == canon_trait(stored):
            continue
        rows_out.append({"column": column, "stored": stored, "normalized": tidy,
                         "n_taxa": n, "action": "folded"})
        print(f"   {column:20} {str(stored)[:40]:40} -> {tidy or 'NULL':24} {n:>5}")
        changed += n
        if not DRY_RUN:
            cur.execute(f"UPDATE taxa SET {column} = %s WHERE {column} = %s",
                        (tidy or None, stored))
    return changed


def resplit_sporulation(cur, rows_out) -> Counter:
    """Move cell-arrangement text out of sporulation; reduce the rest to yes/no."""
    counts: Counter = Counter()
    cur.execute("""SELECT sporulation, COUNT(*),
                          SUM(cell_arrangement IS NULL OR TRIM(cell_arrangement) = '')
                   FROM taxa
                   WHERE sporulation IS NOT NULL AND TRIM(sporulation) <> ''
                   GROUP BY sporulation""")
    for stored, n, n_free in cur.fetchall():
        spore, arrangement = split_sporulation(stored)
        if not spore and not arrangement:
            print(f"   sporulation {str(stored)[:44]:44} unrecognised, left alone "
                  f"({n})")
            rows_out.append({"column": "sporulation", "stored": stored,
                             "normalized": "", "n_taxa": n, "action": "unmatched"})
            counts["unmatched"] += n
            continue
        if arrangement:
            # n_free is the honest part of this. For all but a couple of these
            # rows cell_arrangement ALREADY holds the same reading, so the move
            # copies nothing and only clears the wrong column -- which means the
            # information was never only in sporulation.
            print(f"   sporulation {str(stored)[:30]:30} -> cell_arrangement: "
                  f"{n} taxa, {n_free} need the copy, "
                  f"{n - n_free} already recorded there")
            rows_out.append({"column": "sporulation", "stored": stored,
                             "normalized": arrangement, "n_taxa": n,
                             "action": "moved_to_cell_arrangement"})
            counts["moved"] += n
            if not DRY_RUN:
                # Only where cell_arrangement is empty: a value recovered from
                # the wrong column must never overwrite a real reading.
                cur.execute("""UPDATE taxa SET cell_arrangement = %s
                               WHERE sporulation = %s
                                 AND (cell_arrangement IS NULL
                                      OR TRIM(cell_arrangement) = '')""",
                            (arrangement, stored))
                cur.execute("UPDATE taxa SET sporulation = NULL "
                            "WHERE sporulation = %s", (stored,))
            continue
        if canon_trait(stored) != spore:
            print(f"   sporulation {str(stored)[:44]:44} -> {spore:5} {n:>5}")
            rows_out.append({"column": "sporulation", "stored": stored,
                             "normalized": spore, "n_taxa": n, "action": "folded"})
            counts["folded"] += n
            if not DRY_RUN:
                cur.execute("UPDATE taxa SET sporulation = %s WHERE sporulation = %s",
                            (spore, stored))
    return counts


def main() -> int:
    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    rows: list[dict] = []

    before = {}
    columns = list(TRAIT_VOCABULARIES) + list(FREE_TEXT_TRAITS) + ["sporulation"]
    for column in columns:
        cur.execute(f"SELECT COUNT(DISTINCT {column}) FROM taxa "
                    f"WHERE {column} IS NOT NULL AND TRIM({column}) <> ''")
        before[column] = cur.fetchone()[0]

    print("closed vocabularies:")
    total = 0
    for column, mapping in TRAIT_VOCABULARIES.items():
        total += normalize_column(cur, column, mapping, rows)

    print("\nfree-text columns (whitespace and case only, no vocabulary):")
    for column in FREE_TEXT_TRAITS:
        total += normalize_column(cur, column, None, rows)

    print("\nsporulation, which holds three different kinds of value:")
    spore_counts = resplit_sporulation(cur, rows)

    if not DRY_RUN:
        conn.commit()

    print("\ndistinct values per column, before -> after:")
    for column in columns:
        cur.execute(f"SELECT COUNT(DISTINCT {column}) FROM taxa "
                    f"WHERE {column} IS NOT NULL AND TRIM({column}) <> ''")
        after = cur.fetchone()[0]
        flag = "" if after == before[column] else "   <-"
        print(f"   {column:22} {before[column]:>4} -> {after:>4}{flag}")

    print(f"\ntaxa rewritten: {total}; sporulation {dict(spore_counts)}")
    # shape keeps variants like 'Yeast-like cell/globose to elongate' beside
    # 'Yeast-like cell, globose to elongate'. They differ by punctuation inside
    # free prose, not by case or whitespace, and there is no closed vocabulary for
    # shape to fold them onto -- inventing one would be a taxonomic judgement
    # rather than a repair, so they are reported and left.
    cur.execute("""SELECT COUNT(DISTINCT shape) FROM taxa
                   WHERE shape LIKE 'Yeast-like%'""")
    print(f"shape: {cur.fetchone()[0]} distinct yeast-cell descriptions remain; "
          f"free prose with no vocabulary to fold onto")

    AUDIT.parent.mkdir(parents=True, exist_ok=True)
    with AUDIT.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=AUDIT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {AUDIT} ({len(rows)} rows)")
    conn.close()

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
