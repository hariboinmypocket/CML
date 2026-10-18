#!/usr/bin/env python3
"""Repair the mirror's three mangled bracket names, without moving any mass.

NCBI writes a misapplied name in square brackets: "[Clostridium] methylpentosum"
means the organism is not a Clostridium but keeps the epithet. The mirror holds
three of those with the OPENING bracket lost, which is a parsing artifact, not a
taxonomic statement:

    id 2247  'Clostridium] methylpentosum'  genus rank  372 abundance rows
    id 2729  'Bacteroides] pectinophilus'   genus rank  142
    id 4533  'Eubacterium] saphenum'        genus rank   19

Both scientific_name and genus hold the whole mangled string, so the row is a
genus-layer feature named after a species.

WHY THIS DOES NOT DO WHAT THE MYSQL SCRIPT DID. fix_genus_layer_purity.py
repointed such rows to their leading name token -- the "genus" -- summing into
that genus where the sample already held it. The brackets are precisely the
warning against that. NCBI assigns NO GENUS to any of these three; they are
incertae sedis, and all three are Bacillota:

    84026   [Clostridium] methylpentosum   family Oscillospiraceae, no genus
    384638  [Bacteroides] pectinophilus    no family, no genus
    51123   [Eubacterium] saphenum         family Anaerovoracaceae, no genus

So repointing put a Firmicute under the Bacteroidetes genus Bacteroides -- 142
rows averaging 1.5% of a sample's genus composition and reaching 11.6%, in the
most-measured genus in this database. Replicating that here to make the stores
match would be propagating an error, not fixing one.

WHAT THIS DOES INSTEAD: restores the bracket so the name says what NCBI says,
and leaves the abundance exactly where it is. Mass is untouched, so no rank group
changes and no renormalisation is needed. The row stays a genus-layer feature
that is honestly named and, via taxa.ncbi_rank, honestly described as something
NCBI ranks as a species with no genus -- which v_abundance_genus already exposes
as coarser_than_genus for features like it.

The same tax ID still appears on a species-rank row, so these three remain
duplicate-tax-ID groups. That is correct: one is genus-layer mass, the other a
species-layer taxon, and merging across ranks is what broke 663 rank groups when
merge_duplicate_taxa_mirror.py first ran. They stay listed in
data/taxa_merge_mirror_audit.csv as spans_two_ranks.

Mirror only -- MySQL has no such rows. Idempotent. Dry run by default.
"""
from __future__ import annotations

import os
import re
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

DRY_RUN = "--apply" not in sys.argv
SQLITE = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")

# A name with a closing bracket and no opening one: 'Clostridium] methylpentosum'.
MANGLED = re.compile(r"^\s*([A-Z][A-Za-z-]*)\]\s*(.+?)\s*$")


def repair_name(value: str) -> str | None:
    """'Clostridium] methylpentosum' -> '[Clostridium] methylpentosum'."""
    if not value or "]" not in value or "[" in value:
        return None
    m = MANGLED.match(value)
    if not m:
        return None
    return f"[{m.group(1)}] {m.group(2)}"


def main() -> int:
    if not Path(SQLITE).exists():
        print(f"sqlite mirror {SQLITE} not present; nothing to do")
        return 0

    conn = sqlite3.connect(SQLITE)
    cur = conn.cursor()
    cur.execute("""SELECT id, scientific_name, genus, species, taxonomic_rank,
                          ncbi_tax_id
                   FROM taxa
                   WHERE scientific_name LIKE '%]%' OR genus LIKE '%]%'
                      OR species LIKE '%]%'""")
    rows = cur.fetchall()
    candidates = [r for r in rows if repair_name(r[1]) or repair_name(r[2] or "")]
    print(f"taxa with a closing bracket and no opening one: {len(candidates)}")
    if not candidates:
        print("nothing to do")
        conn.close()
        return 0

    touched = 0
    for taxon_id, sci, genus, species, rank, taxid in candidates:
        new_sci = repair_name(sci) or sci
        new_genus = repair_name(genus or "") or genus
        cur.execute("""SELECT COUNT(*) FROM sample_taxon_abundances
                       WHERE taxon_id = ?""", (taxon_id,))
        n_ab = cur.fetchone()[0]
        print(f"   id={taxon_id} rank={rank} taxid={taxid} abundances={n_ab}")
        print(f"      scientific_name {sci!r} -> {new_sci!r}")
        if new_genus != genus:
            print(f"      genus           {genus!r} -> {new_genus!r}")
        print(f"      abundance rows: unchanged (mass stays in the genus layer)")
        touched += 1
        if not DRY_RUN:
            cur.execute("""UPDATE taxa SET scientific_name = ?, genus = ?
                           WHERE id = ?""", (new_sci, new_genus, taxon_id))

    if not DRY_RUN:
        conn.commit()

    cur.execute("""SELECT COUNT(*) FROM (
                     SELECT a.sample_id, t.taxonomic_rank,
                            SUM(a.relative_abundance) total
                     FROM sample_taxon_abundances a
                     JOIN taxa t ON t.id = a.taxon_id
                     GROUP BY a.sample_id, t.taxonomic_rank
                     HAVING ABS(total - 1.0) > 0.01) bad""")
    print(f"\nrank groups not summing to 1: {cur.fetchone()[0]} "
          f"(must be 0 -- no mass was moved)")
    cur.execute("""SELECT COUNT(*) FROM taxa
                   WHERE (scientific_name LIKE '%]%' OR genus LIKE '%]%')
                     AND scientific_name NOT LIKE '[%'""")
    print(f"mangled bracket names remaining: {cur.fetchone()[0]}")
    conn.close()

    print(f"\nnames repaired: {touched}")
    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
