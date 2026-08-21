#!/usr/bin/env python3
"""Move stray species-rank rows out of the genus abundance layer.

Every abundance currently stored was loaded as a genus-level layer (GMrepo,
MGnify and BioMapAI alike). A handful still point at species-rank taxa because
GMrepo's genus array contains binomials such as "[Clostridium] methylpentosum",
which parse to genus+epithet. The value is correct; only the taxon it hangs off
is at the wrong rank, which breaks the per-(sample, rank) sum-to-1 invariant.

Each such row is repointed to its genus, summing where the sample already holds
that genus. The species taxon itself is left in place — it is legitimate as a
taxonomic entity and may carry disease associations; only the abundance rows
move. Dry run by default.
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
