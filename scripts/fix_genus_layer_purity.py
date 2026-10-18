#!/usr/bin/env python3
"""Move stray species-rank rows out of the genus abundance layer.

DO NOT RUN THIS WITHOUT READING THE NEXT PARAGRAPH. Its premise expired.

It was written when every abundance in the database was a genus-level layer, so
"taxonomic_rank <> 'genus'" identified a stray exactly. The species layer was
loaded afterwards (data/gmrepo_species_abundances.csv, August 2026) and now holds
1,277,809 legitimate rows, every one of which that selector matches. A dry run
today reports "rows moved to genus: 1277809" -- it would collapse the entire
species layer into the genus layer and add species mass to genus values. The
script is kept for its history and for the narrower job below, and it now
REFUSES to run when the species layer is present rather than relying on whoever
reads the output to notice.

The original defect: GMrepo's genus array contains binomials such as
"[Clostridium] methylpentosum", which parse to genus+epithet and were stored at
species rank while their abundance belonged to the genus layer.

THE REPOINT IT PERFORMED IS ALSO QUESTIONABLE, which the brackets should have
warned about. NCBI's square brackets mean the name is MISAPPLIED, and for all
three of these organisms NCBI assigns no genus at all -- they are incertae sedis,
and all three are Bacillota:

    [Clostridium] methylpentosum  84026   family Oscillospiraceae, no genus
    [Bacteroides] pectinophilus  384638   no family, no genus
    [Eubacterium] saphenum        51123   family Anaerovoracaceae, no genus

So repointing to the leading name token put a Firmicute's mass under the
Bacteroidetes genus Bacteroides. Measured in the SQLite mirror, which still holds
the un-repointed rows: 142 rows averaging 1.5% of a sample's genus composition
and reaching 11.6%. That is not a rounding error, and Bacteroides is the
most-measured genus here. Do not extend this behaviour to another store; see
scripts/fix_mirror_bracket_names.py, which cleans the names and leaves the mass
where it is.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.pipeline import upsert_taxon

DRY_RUN = "--apply" not in sys.argv


def main() -> int:
    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()

    # The selector below matches every species-layer row. Refuse rather than
    # trust the operator to read the count.
    cur.execute("""SELECT COUNT(*) FROM sample_taxon_abundances a
                   JOIN taxa t ON t.id = a.taxon_id
                   WHERE t.taxonomic_rank = 'species'""")
    species_rows = cur.fetchone()[0]
    if species_rows > 1000:
        print(f"REFUSING TO RUN: the species layer holds {species_rows:,} rows.")
        print("This script's selector (taxonomic_rank <> 'genus') matches all of")
        print("them, so it would collapse the species layer into the genus layer.")
        print("Its premise -- that every abundance is a genus-level layer -- has")
        print("not held since the species layer was loaded. See the module")
        print("docstring.")
        cur.close()
        conn.close()
        return 1

    cur.execute(
        """SELECT a.sample_id, a.taxon_id, a.relative_abundance, t.genus, t.scientific_name
           FROM sample_taxon_abundances a
           JOIN taxa t ON t.id = a.taxon_id
           WHERE t.taxonomic_rank <> 'genus'"""
    )
    strays = cur.fetchall()
    print(f"abundance rows at a non-genus rank: {len(strays)}")
    if not strays:
        print("nothing to do")
        return 0

    by_name: dict[str, int] = {}
    for _s, _t, _v, _g, name in strays:
        by_name[name] = by_name.get(name, 0) + 1
    for name, count in sorted(by_name.items(), key=lambda kv: -kv[1])[:10]:
        print(f"   {count:5d}  {name}")

    moved = summed = 0
    genus_cache: dict[str, int] = {}

    for sample_id, taxon_id, value, genus, _name in strays:
        if not genus:
            continue
        if DRY_RUN:
            moved += 1
            continue
        if genus not in genus_cache:
            gid, _ = upsert_taxon(
                conn,
                {"genus": genus, "species": "", "taxonomic_rank": "genus",
                 "source_database": "GMrepo"},
                overwrite=False,
            )
            genus_cache[genus] = gid
        target = genus_cache[genus]

        cur.execute(
            "SELECT relative_abundance FROM sample_taxon_abundances "
            "WHERE sample_id=%s AND taxon_id=%s", (sample_id, target),
        )
        existing = cur.fetchone()
        if existing:
            cur.execute(
                "UPDATE sample_taxon_abundances SET relative_abundance=%s "
                "WHERE sample_id=%s AND taxon_id=%s",
                (float(existing[0]) + float(value), sample_id, target),
            )
            summed += 1
        else:
            cur.execute(
                "INSERT INTO sample_taxon_abundances (sample_id, taxon_id, relative_abundance) "
                "VALUES (%s,%s,%s)", (sample_id, target, value),
            )
        cur.execute(
            "DELETE FROM sample_taxon_abundances WHERE sample_id=%s AND taxon_id=%s",
            (sample_id, taxon_id),
        )
        moved += 1

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    else:
        conn.commit()
    print(f"  rows moved to genus : {moved}")
    print(f"  of which summed into an existing genus row : {summed}")
    cur.close()
    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
