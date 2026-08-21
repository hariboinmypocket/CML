#!/usr/bin/env python3
"""Repair taxa split under a fabricated genus called "Candidatus".

"Candidatus" is a nomenclatural status prefix for uncultured organisms, not a
genus. A two-word name such as "Candidatus Soleaferrea" is a bare genus, but the
parser previously only handled the three-word form, so these became
genus="Candidatus" / species="soleaferrea" — a genus that does not exist, with
unrelated organisms collapsed beneath it and an incorrect species rank that also
broke the per-rank abundance invariant.

Where the corrected name already exists the rows are merged and dependents
repointed; otherwise the row is renamed in place. Dry run by default.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.transform import normalize_genus

DRY_RUN = "--apply" not in sys.argv


def main() -> int:
    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()

    cur.execute("SELECT COUNT(*) FROM taxa")
    before = cur.fetchone()[0]
    cur.execute(
        "SELECT id, scientific_name, genus, species, taxonomic_rank "
        "FROM taxa WHERE genus = 'Candidatus' ORDER BY id"
    )
    targets = cur.fetchall()
    print(f"taxa with fabricated genus 'Candidatus': {len(targets)}")

    merged = renamed = repointed_ab = repointed_as = 0

    for taxon_id, sci_name, _genus, species, _rank in targets:
        if not species:
            continue
        correct_genus = normalize_genus(f"Candidatus {species}")
        genus_key, species_key = correct_genus.casefold(), ""

        cur.execute(
            "SELECT id FROM taxa WHERE genus_key=%s AND species_key=%s AND id<>%s",
            (genus_key, species_key, taxon_id),
        )
        existing = cur.fetchone()

        if existing:
            canonical = existing[0]
            if not DRY_RUN:
                # Abundances are keyed (sample, taxon). Where a sample holds both
                # the wrong and the correct taxon the values must be SUMMED, not
                # overwritten, or that sample's rank total silently drops. Done in
                # Python because MySQL cannot reference the target table inside
                # ON DUPLICATE KEY UPDATE when it is also the source.
                cur.execute(
                    "SELECT sample_id, relative_abundance FROM sample_taxon_abundances "
                    "WHERE taxon_id=%s", (taxon_id,),
                )
                for sample_id, value in cur.fetchall():
                    cur.execute(
                        "SELECT relative_abundance FROM sample_taxon_abundances "
                        "WHERE sample_id=%s AND taxon_id=%s", (sample_id, canonical),
                    )
                    existing_value = cur.fetchone()
                    if existing_value:
                        cur.execute(
                            "UPDATE sample_taxon_abundances SET relative_abundance=%s "
                            "WHERE sample_id=%s AND taxon_id=%s",
                            (float(existing_value[0]) + float(value), sample_id, canonical),
                        )
                    else:
                        cur.execute(
                            "INSERT INTO sample_taxon_abundances "
                            "(sample_id, taxon_id, relative_abundance) VALUES (%s,%s,%s)",
                            (sample_id, canonical, value),
                        )
                    repointed_ab += 1
                cur.execute("DELETE FROM sample_taxon_abundances WHERE taxon_id=%s", (taxon_id,))
                cur.execute(
                    "UPDATE IGNORE taxon_disease_associations SET taxon_id=%s WHERE taxon_id=%s",
                    (canonical, taxon_id),
                )
                repointed_as += cur.rowcount
                cur.execute("DELETE FROM taxon_disease_associations WHERE taxon_id=%s", (taxon_id,))
                cur.execute("DELETE FROM taxa WHERE id=%s", (taxon_id,))
            merged += 1
        else:
            if not DRY_RUN:
                cur.execute(
                    """UPDATE taxa SET scientific_name=%s, genus=%s, species='',
                           genus_key=%s, species_key='', taxonomic_rank='genus'
                       WHERE id=%s""",
                    (correct_genus, correct_genus, genus_key, taxon_id),
                )
            renamed += 1

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    else:
        conn.commit()

    print(f"  merged into existing taxon : {merged}")
    print(f"  renamed in place           : {renamed}")
    print(f"  abundance rows repointed   : {repointed_ab}")
    print(f"  association rows repointed : {repointed_as}")
    cur.execute("SELECT COUNT(*) FROM taxa")
    print(f"  taxa {before} -> {cur.fetchone()[0]}")
    cur.close()
    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
