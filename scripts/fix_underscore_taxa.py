#!/usr/bin/env python3
"""Repair MGnify-derived taxa whose names kept MGnify's underscore separator.

"Abiotrophia_defectiva" was stored as a single-token genus, which (a) duplicates
the properly-named row already loaded from another source, (b) blocks NCBI
Taxonomy lookups, and (c) splits one organism across two ML features.

Where a properly-named row already exists the two are merged and dependent rows
repointed; otherwise the row is renamed in place with genus/species re-derived.
Runs in one transaction and reports before/after counts.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.transform import parse_scientific_name

DRY_RUN = "--apply" not in sys.argv


def main() -> int:
    settings = Settings.from_env(".env")
    conn = connect(settings)
    cur = conn.cursor()

    cur.execute("SELECT COUNT(*) FROM taxa")
    taxa_before = cur.fetchone()[0]
    cur.execute("SELECT COUNT(*) FROM sample_taxon_abundances")
    abund_before = cur.fetchone()[0]
    cur.execute("SELECT COUNT(*) FROM taxon_disease_associations")
    assoc_before = cur.fetchone()[0]

    cur.execute(
        "SELECT id, scientific_name, source_database FROM taxa "
        "WHERE scientific_name LIKE %s ORDER BY id",
        ("%\\_%",),
    )
    targets = cur.fetchall()
    print(f"underscore taxa found: {len(targets)}")

    merged = renamed = skipped = repointed_abund = repointed_assoc = 0

    for taxon_id, name, _source in targets:
        clean = name.replace("_", " ").strip()
        genus, species = parse_scientific_name(clean)
        if not genus:
            skipped += 1
            continue
        genus_key, species_key = genus.casefold(), species.casefold()

        cur.execute(
            "SELECT id FROM taxa WHERE genus_key=%s AND species_key=%s AND id<>%s",
            (genus_key, species_key, taxon_id),
        )
        existing = cur.fetchone()

        if existing:
            canonical = existing[0]
            if not DRY_RUN:
                cur.execute(
                    "UPDATE IGNORE sample_taxon_abundances SET taxon_id=%s WHERE taxon_id=%s",
                    (canonical, taxon_id),
                )
                repointed_abund += cur.rowcount
                cur.execute(
                    "UPDATE IGNORE taxon_disease_associations SET taxon_id=%s WHERE taxon_id=%s",
                    (canonical, taxon_id),
                )
                repointed_assoc += cur.rowcount
                # anything left is a genuine duplicate row that would violate the
                # unique key; drop it rather than leaving it orphaned
                cur.execute("DELETE FROM sample_taxon_abundances WHERE taxon_id=%s", (taxon_id,))
                cur.execute("DELETE FROM taxon_disease_associations WHERE taxon_id=%s", (taxon_id,))
                cur.execute("DELETE FROM taxa WHERE id=%s", (taxon_id,))
            merged += 1
        else:
            if not DRY_RUN:
                cur.execute(
                    """UPDATE taxa
                       SET scientific_name=%s, genus=%s, species=%s,
                           genus_key=%s, species_key=%s,
                           taxonomic_rank=%s
                       WHERE id=%s""",
                    (clean, genus, species, genus_key, species_key,
                     "species" if species else "genus", taxon_id),
                )
            renamed += 1

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    else:
        conn.commit()

    print(f"  merged into existing taxon : {merged}")
    print(f"  renamed in place           : {renamed}")
    print(f"  skipped (unparseable)      : {skipped}")
    print(f"  abundance rows repointed   : {repointed_abund}")
    print(f"  association rows repointed : {repointed_assoc}")

    cur.execute("SELECT COUNT(*) FROM taxa")
    taxa_after = cur.fetchone()[0]
    cur.execute("SELECT COUNT(*) FROM sample_taxon_abundances")
    abund_after = cur.fetchone()[0]
    cur.execute("SELECT COUNT(*) FROM taxon_disease_associations")
    assoc_after = cur.fetchone()[0]
    print(f"\n  taxa         {taxa_before} -> {taxa_after}")
    print(f"  abundances   {abund_before} -> {abund_after}")
    print(f"  associations {assoc_before} -> {assoc_after}")

    cur.close()
    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
