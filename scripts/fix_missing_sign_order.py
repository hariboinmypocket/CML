#!/usr/bin/env python3
"""Record the sign order wherever a scored association lacks one.

An LDA score's sign is meaningless without its comparison order -- this schema's
central design point -- and phenotype_comparisons stores which phenotype a
POSITIVE score refers to for exactly that reason. One association per store has a
magnitude and no sign order, so its number cannot be read:

    MySQL   Bacteroides / Vitiligo, LDA 4.31, pmid39967732|health|vitiligo|lefse|genus
    mirror  the same finding, under PMID39967732|health|vitiligo|lefse

56 comparisons in MySQL and 45 in the mirror carry no sign order, and that is
fine: the rest hold only unscored associations -- a direction with no magnitude,
where sign order has nothing to disambiguate. The gap is one scored row each.

WHY DERIVING THE ARMS FROM THE DIRECTION IS NOT CIRCULAR. The objection is real:
association_direction() normally derives `direction` FROM the sign order, so
going back the other way would assert nothing. It does not apply here, and the
reason is structural rather than a judgement about this row. That function
returns 'marker' when it has no enriched_group to work from -- so a row with a
score, no sign order, and a direction of 'enriched' or 'depleted' CANNOT have had
that direction derived from the sign order. It must have come from an explicit
field in the source, and for this row it did: data/literature_vitiligo_markers.csv
carries a `direction` column reading `enriched` with both enriched_in columns
empty, and a note saying "Highest LDA score among vitiligo-group markers".

So the rule is: a comparison gets its arms set when it has at least one scored
association, no sign order at all, exactly two phenotypes one of which is Health,
and every scored association agrees about which group a positive score belongs
to. Anything else is left alone and reported.

A MIXED SET OF DIRECTIONS IS NOT AMBIGUITY. LEfSe reports markers of both
groups, so the MySQL comparison holds 5 'enriched' rows and 3 'depleted' ones.
An earlier version of this script refused on that, asking the wrong question:
only SCORED associations constrain the sign order, and the unscored ones take
their direction from the source.

Both stores. Idempotent. Dry run by default; pass --apply to write.
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

# Scored associations with no sign order, and everything needed to judge them.
CANDIDATES = """
SELECT c.id, c.comparison_key, c.phenotype_a_id, c.phenotype_b_id,
       COUNT(*) AS n_scored,
       SUM(CASE WHEN a.effect_size > 0 AND a.direction = 'enriched' THEN 1
                WHEN a.effect_size < 0 AND a.direction = 'depleted' THEN 1
                ELSE 0 END) AS n_positive_means_disease,
       SUM(CASE WHEN a.direction NOT IN ('enriched', 'depleted') THEN 1
                ELSE 0 END) AS n_undirected,
       MIN(a.disease_id) AS disease_id,
       MAX(a.disease_id) AS disease_id_max
FROM taxon_disease_associations a
JOIN phenotype_comparisons c ON c.id = a.comparison_id
WHERE a.effect_size IS NOT NULL AND a.effect_size <> 0
  AND c.positive_score_enriched_in_id IS NULL
  AND c.negative_score_enriched_in_id IS NULL
GROUP BY c.id, c.comparison_key, c.phenotype_a_id, c.phenotype_b_id
"""


def health_id(cur, placeholder: str):
    cur.execute(f"SELECT id FROM diseases WHERE name = {placeholder}", ("Health",))
    row = cur.fetchone()
    return row[0] if row else None


def repair(label: str, cur, placeholder: str) -> tuple[int, int]:
    health = health_id(cur, placeholder)
    if health is None:
        print(f"   {label}: no 'Health' disease row; cannot judge any arms")
        return 0, 0

    cur.execute(CANDIDATES)
    rows = cur.fetchall()
    if not rows:
        print(f"   {label}: no scored association lacks a sign order")
        return 0, 0

    written = refused = 0
    for (cid, key, a_id, b_id, n_scored, n_pos_disease, n_undirected,
         disease_id, disease_id_max) in rows:
        why = None
        if n_undirected:
            why = f"{n_undirected} of {n_scored} scored rows carry no direction"
        elif disease_id != disease_id_max:
            why = "scored rows span more than one disease"
        elif health not in (a_id, b_id):
            why = f"neither arm is Health (a={a_id}, b={b_id})"
        elif disease_id not in (a_id, b_id):
            why = f"the disease {disease_id} is not one of the arms"
        elif n_pos_disease not in (0, n_scored):
            why = (f"{n_pos_disease} of {n_scored} scored rows say a positive "
                   f"score means the disease, the rest say the opposite")
        if why:
            print(f"   {label}: REFUSED {key}: {why}")
            refused += 1
            continue

        # Unanimous. All scored rows agreeing that a positive score means the
        # disease puts the disease on the positive arm; all disagreeing puts it
        # on the negative one.
        if n_pos_disease == n_scored:
            positive, negative = disease_id, health
        else:
            positive, negative = health, disease_id
        print(f"   {label}: {key}")
        print(f"      {n_scored} scored, positive -> {positive}, "
              f"negative -> {negative}")
        written += 1
        if not DRY_RUN:
            cur.execute(f"""UPDATE phenotype_comparisons
                            SET positive_score_enriched_in_id = {placeholder},
                                negative_score_enriched_in_id = {placeholder}
                            WHERE id = {placeholder}""", (positive, negative, cid))
    return written, refused


def main() -> int:
    total_written = total_refused = 0

    print("mysql:")
    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    w, r = repair("mysql", cur, "%s")
    if not DRY_RUN:
        conn.commit()
    conn.close()
    total_written += w
    total_refused += r

    path = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")
    if Path(path).exists():
        print("\nsqlite mirror:")
        sq = sqlite3.connect(path)
        sc = sq.cursor()
        w, r = repair("sqlite", sc, "?")
        if not DRY_RUN:
            sq.commit()
        sq.close()
        total_written += w
        total_refused += r
    else:
        print(f"\nsqlite mirror {path} not present; skipped")

    print(f"\nset: {total_written}   refused: {total_refused}")
    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    return 1 if total_refused else 0


if __name__ == "__main__":
    raise SystemExit(main())
