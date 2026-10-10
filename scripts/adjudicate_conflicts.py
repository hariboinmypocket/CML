#!/usr/bin/env python3
"""Adjudicate taxon-disease pairs whose studies disagree on direction.

2,285 (taxon, disease) pairs are reported by more than one study, and 911 of
them conflict: 485 agree in more than 60% of studies, 39 are near ties, and 387
are exact ties where the association vote cannot decide anything. A model fed
`consensus_direction` takes a coin flip on those as though it were a finding.

The adjudicator brings in evidence the vote does not use: the abundance matrix
in this same database. For each study holding both arms, the taxon's mean
relative abundance in that study's case samples is compared with its mean in
that study's controls, giving one direction per study. Those are then counted
the same way the association votes are.

Four decisions shape it.

Only case/control contrasts vote. 2,340 of the 16,889 associations compare one
disease against another rather than against Health, and 1,035 of the replicated
pairs mix the two kinds. An earlier v_taxon_disease_evidence voted across both as
though they answered the same question, which is how Faecalibacterium prausnitzii
in ulcerative colitis came out 'enriched': four studies report it enriched in UC
*versus Crohn disease* and two report it depleted *versus Health*. That view now
applies the same restriction this script does, so the association vote is
computed here a second time, independently, and the two are compared on every
pair -- a cross-check that is cheap to keep and would catch either one drifting.
Pairs whose associations are all disease-vs-disease keep their evidence recorded
but get no association vote.

Within-study, never pooled. q7 pools case and control samples across studies
and finds only 62% concordance with the curated associations, which is what
pooling across different cohorts, protocols and sequencing runs produces. A
direction computed inside one study compares like with like, and only the
directions are aggregated.

Absence is zero, not missing. A taxon that was not detected in a sample has no
row in sample_taxon_abundances. Averaging only the rows that exist would
compare the samples where a taxon was abundant against the samples where it was
abundant, which is no comparison at all. The denominator is therefore every
sample in the arm that was profiled at that taxon's rank, with absence counted
as zero.

A mean is a sum over a count, and the two factorize. The count -- how many
samples in an arm were profiled at a rank -- does not depend on which taxon is
being asked about, and the sum only needs the rows that exist. So both come
from one aggregate pass each, and the division happens per pair in Python. The
first draft of this script asked the database once per pair instead, with a
correlated EXISTS over the 2.58M-row abundance table, and took eleven minutes;
this is the same arithmetic in a few seconds.

Nothing in taxon_disease_associations is modified. The verdict lands in
taxon_disease_adjudication beside the evidence for it, so it can be recomputed,
audited, or ignored.

Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import csv
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect

DRY_RUN = "--apply" not in sys.argv
AUDIT = Path("data/adjudication_audit.csv")

# A vote this lopsided is treated as settled; anything flatter is a conflict.
CLEAR_MAJORITY = 0.6
# Below this many studies with usable abundance, the abundance vote is too thin
# to overturn or confirm anything.
MIN_ABUNDANCE_STUDIES = 2
# The control arm. Every within-study comparison is against these samples.
HEALTH = "Health"

# Per-study association directions, restricted to contrasts against Health.
# One row per (taxon, disease, study, direction, contrast kind), so a study
# that reports a direction from both a vs-Health and a disease-vs-disease
# comparison is recorded under both and the scope label stays truthful. Only the
# vs-Health rows vote. A study appearing with both directions among its
# vs-Health comparisons is internally split and votes for neither.
ASSOC_VOTE_SQL = """
SELECT a.taxon_id, a.disease_id, cp.study_id, a.direction,
       COALESCE(da.name = %s OR db.name = %s, 0) AS vs_health
FROM taxon_disease_associations a
JOIN phenotype_comparisons cp ON cp.id = a.comparison_id
LEFT JOIN diseases da ON da.id = cp.phenotype_a_id
LEFT JOIN diseases db ON db.id = cp.phenotype_b_id
WHERE a.direction IN ('enriched', 'depleted')
GROUP BY a.taxon_id, a.disease_id, cp.study_id, a.direction, vs_health
"""

# How many samples each (study, disease) arm contributes at each rank. This is
# the denominator of every mean, and it is taxon-independent.
ARM_SIZE_SQL = """
SELECT s.study_id, s.disease_id, tx.taxonomic_rank, COUNT(DISTINCT ab.sample_id)
FROM sample_taxon_abundances ab
JOIN taxa tx ON tx.id = ab.taxon_id
JOIN samples s ON s.id = ab.sample_id
GROUP BY s.study_id, s.disease_id, tx.taxonomic_rank
"""

# The numerator: total abundance of each taxon of interest within each
# (study, disease) arm. Absent samples contribute nothing, which is what
# counting absence as zero means.
TAXON_SUM_SQL = """
SELECT ab.taxon_id, s.study_id, s.disease_id, SUM(ab.relative_abundance)
FROM sample_taxon_abundances ab
JOIN samples s ON s.id = ab.sample_id
WHERE ab.taxon_id IN ({placeholders})
GROUP BY ab.taxon_id, s.study_id, s.disease_id
"""


def association_vote(per_study: dict) -> dict:
    """Count one vote per study, from vs-Health contrasts only.

    `per_study` maps study_id -> set of directions seen in that study's
    case/control comparisons. A study holding both directions is split and
    counted as neither, rather than voting twice as the view lets it.
    """
    enriched = depleted = split = 0
    for directions in per_study.values():
        if len(directions) > 1:
            split += 1
        elif "enriched" in directions:
            enriched += 1
        else:
            depleted += 1
    total = enriched + depleted
    ratio = (max(enriched, depleted) / total) if total else None
    return {
        "assoc_enriched": enriched,
        "assoc_depleted": depleted,
        "assoc_split_studies": split,
        "assoc_consensus": (
            None if not total or enriched == depleted
            else "enriched" if enriched > depleted else "depleted"
        ),
        "assoc_agreement": round(ratio, 3) if ratio is not None else None,
        "assoc_studies": total,
    }


def abundance_vote(taxon_id: int, disease_id: int, rank: str,
                   arm_size: dict, taxon_sum: dict, health_id: int) -> dict:
    """One direction per study that holds both arms, then counted."""
    enriched = depleted = ties = 0
    if disease_id != health_id:
        controls = arm_size.get((health_id, rank), {})
        for study_id, n_case in arm_size.get((disease_id, rank), {}).items():
            n_health = controls.get(study_id)
            if not n_health or not n_case:
                continue
            case = taxon_sum.get((taxon_id, study_id, disease_id), 0.0) / n_case
            health = taxon_sum.get((taxon_id, study_id, health_id), 0.0) / n_health
            if case == health:
                ties += 1
            elif case > health:
                enriched += 1
            else:
                depleted += 1
    total = enriched + depleted
    ratio = (max(enriched, depleted) / total) if total else None
    return {
        "ab_enriched": enriched,
        "ab_depleted": depleted,
        "ab_ties": ties,
        "ab_studies": total,
        "ab_ratio": round(ratio, 3) if ratio is not None else None,
        "ab_direction": (
            None if not total or enriched == depleted
            else "enriched" if enriched > depleted else "depleted"
        ),
    }


def adjudicate(assoc: dict, ab: dict) -> tuple[str | None, str]:
    """Return (verdict, basis). Honest about what cannot be decided."""
    if assoc["assoc_studies"] == 0:
        # Every association for this pair compares one disease against another,
        # so nothing here speaks to the disease-vs-Health direction a model
        # needs. Abundance is the only witness.
        if ab["ab_studies"] < MIN_ABUNDANCE_STUDIES:
            return None, "unresolved_no_case_control_evidence"
        if ab["ab_direction"] is None:
            return None, "unresolved_abundance_also_tied"
        if ab["ab_ratio"] is not None and ab["ab_ratio"] <= CLEAR_MAJORITY:
            return None, "unresolved_abundance_inconclusive"
        return ab["ab_direction"], "abundance_only_no_case_control_association"

    settled = assoc["assoc_agreement"] is not None and assoc["assoc_agreement"] > CLEAR_MAJORITY
    if settled:
        if ab["ab_direction"] is None or ab["ab_studies"] < MIN_ABUNDANCE_STUDIES:
            return assoc["assoc_consensus"], "association_consensus"
        if ab["ab_direction"] == assoc["assoc_consensus"]:
            return assoc["assoc_consensus"], "association_and_abundance_agree"
        # The studies mostly agree and the measurements say otherwise. Not
        # overruled -- flagged, because one of the two is wrong and which is not
        # decidable from here.
        return assoc["assoc_consensus"], "association_consensus_abundance_disagrees"

    # The association vote is a tie or near-tie, so abundance decides if it can.
    if ab["ab_studies"] < MIN_ABUNDANCE_STUDIES:
        return None, "unresolved_no_abundance_evidence"
    if ab["ab_direction"] is None:
        return None, "unresolved_abundance_also_tied"
    if ab["ab_ratio"] is not None and ab["ab_ratio"] <= CLEAR_MAJORITY:
        return None, "unresolved_abundance_inconclusive"
    return ab["ab_direction"], "abundance_majority"


def main() -> int:
    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    work = conn.cursor()

    cur.execute("SELECT id FROM diseases WHERE name = %s", (HEALTH,))
    row = cur.fetchone()
    if not row:
        print(f"no '{HEALTH}' disease row -- nothing to compare against")
        return 1
    health_id = row[0]

    cur.execute(
        """SELECT e.taxon_id, e.scientific_name, t.taxonomic_rank, e.disease_id,
                  e.disease_name, e.n_studies, e.n_studies_enriched,
                  e.n_studies_depleted, e.consensus_direction, e.agreement_ratio
           FROM v_taxon_disease_evidence e
           JOIN taxa t ON t.id = e.taxon_id
           WHERE e.n_studies > 1
           ORDER BY e.n_studies DESC, e.taxon_id""")
    pairs = cur.fetchall()
    print(f"replicated (taxon, disease) pairs: {len(pairs)}")

    started = time.time()
    # (taxon, disease) -> {study_id: {directions}} for vs-Health contrasts, and
    # the set of studies contributing only disease-vs-disease contrasts.
    vs_health: dict[tuple[int, int], dict[int, set]] = defaultdict(lambda: defaultdict(set))
    other_contrast: dict[tuple[int, int], set] = defaultdict(set)
    cur.execute(ASSOC_VOTE_SQL, (HEALTH, HEALTH))
    for taxon_id, disease_id, study_id, direction, is_vs_health in cur.fetchall():
        key = (taxon_id, disease_id)
        if is_vs_health:
            vs_health[key][study_id].add(direction)
        else:
            other_contrast[key].add(study_id)
    print(f"association votes: {len(vs_health)} pairs with case/control contrasts, "
          f"{len(other_contrast)} with disease-vs-disease  ({time.time() - started:.0f}s)")

    # (disease_id, rank) -> {study_id: n_samples_profiled_at_that_rank}
    arm_size: dict[tuple[int, str], dict[int, int]] = defaultdict(dict)
    cur.execute(ARM_SIZE_SQL)
    for study_id, disease_id, rank, n in cur.fetchall():
        if disease_id is not None:
            arm_size[(disease_id, rank)][study_id] = int(n)
    print(f"arm sizes: {sum(len(v) for v in arm_size.values())} "
          f"(study, disease, rank) cells  ({time.time() - started:.0f}s)")

    taxon_ids = sorted({p[0] for p in pairs})
    taxon_sum: dict[tuple[int, int, int], float] = {}
    chunk = 500
    for start in range(0, len(taxon_ids), chunk):
        batch = taxon_ids[start:start + chunk]
        sql = TAXON_SUM_SQL.format(placeholders=",".join(["%s"] * len(batch)))
        cur.execute(sql, batch)
        for taxon_id, study_id, disease_id, total in cur.fetchall():
            if disease_id is not None:
                taxon_sum[(taxon_id, study_id, disease_id)] = float(total or 0.0)
    print(f"taxon sums: {len(taxon_sum)} (taxon, study, disease) cells "
          f"for {len(taxon_ids)} taxa  ({time.time() - started:.0f}s)")

    rows: list[dict] = []
    basis_counts: Counter = Counter()
    scope_counts: Counter = Counter()
    # v_taxon_disease_evidence now applies the same two rules this script does,
    # so its counts and the ones computed here should agree on every pair. They
    # are derived independently -- SQL in the view, Python here -- so a
    # disagreement means one of them has drifted, and saying so is worth more
    # than quietly preferring either.
    view_disagreements: list[str] = []
    for (taxon_id, name, rank, disease_id, disease, n_studies,
         pooled_enr, pooled_dep, pooled_consensus, pooled_agreement) in pairs:
        key = (taxon_id, disease_id)
        case_control = vs_health.get(key, {})
        assoc = association_vote(case_control)
        if case_control and other_contrast.get(key):
            scope = "mixed"
        elif case_control:
            scope = "vs_health"
        elif other_contrast.get(key):
            scope = "disease_vs_disease"
        else:
            scope = "none"
        scope_counts[scope] += 1
        view_cons = None if pooled_consensus == "tied" else pooled_consensus
        if (int(pooled_enr) != assoc["assoc_enriched"]
                or int(pooled_dep) != assoc["assoc_depleted"]
                or view_cons != assoc["assoc_consensus"]):
            view_disagreements.append(
                f"{name} / {disease}: view {pooled_enr}/{pooled_dep} {view_cons}, "
                f"computed {assoc['assoc_enriched']}/{assoc['assoc_depleted']} "
                f"{assoc['assoc_consensus']}")
        ab = abundance_vote(taxon_id, disease_id, rank, arm_size, taxon_sum, health_id)
        verdict, basis = adjudicate(assoc, ab)
        basis_counts[basis] += 1
        rows.append({
            "taxon_id": taxon_id, "scientific_name": name, "taxonomic_rank": rank,
            "disease_id": disease_id, "disease_name": disease,
            "n_studies": n_studies, "assoc_contrast_scope": scope,
            **assoc, **ab, "verdict": verdict, "basis": basis,
            "view_enriched": pooled_enr, "view_depleted": pooled_dep,
            "view_consensus": pooled_consensus,
            "view_agreement": pooled_agreement,
        })
    print(f"adjudicated in {time.time() - started:.0f}s")
    if view_disagreements:
        print(f"\n!! v_taxon_disease_evidence disagrees with this script on "
              f"{len(view_disagreements)} of {len(pairs)} pairs:")
        for line in view_disagreements[:10]:
            print(f"     {line}")
        if len(view_disagreements) > 10:
            print(f"     ... and {len(view_disagreements) - 10} more")
    else:
        print(f"view cross-check: v_taxon_disease_evidence agrees on all "
              f"{len(pairs)} pairs")

    if not DRY_RUN:
        work.execute("DELETE FROM taxon_disease_adjudication")
        work.executemany(
            """INSERT INTO taxon_disease_adjudication
                 (taxon_id, disease_id, n_studies, assoc_contrast_scope,
                  assoc_enriched, assoc_depleted, assoc_split_studies,
                  assoc_consensus, assoc_agreement, abundance_enriched,
                  abundance_depleted, abundance_studies, abundance_agreement,
                  verdict, basis)
               VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)""",
            [(r["taxon_id"], r["disease_id"], r["n_studies"], r["assoc_contrast_scope"],
              r["assoc_enriched"], r["assoc_depleted"], r["assoc_split_studies"],
              r["assoc_consensus"], r["assoc_agreement"],
              r["ab_enriched"], r["ab_depleted"], r["ab_studies"], r["ab_ratio"],
              r["verdict"], r["basis"]) for r in rows])
        conn.commit()
        print(f"wrote {len(rows)} rows to taxon_disease_adjudication")
    else:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")

    print("\nassociation contrast scope:")
    for scope, n in scope_counts.most_common():
        print(f"   {scope:46} {n:>5}")

    print("\nverdict basis:")
    for basis, n in basis_counts.most_common():
        print(f"   {basis:46} {n:>5}")
    resolved = sum(n for b, n in basis_counts.items() if not b.startswith("unresolved"))
    print(f"\nresolved: {resolved}   unresolved: {len(rows) - resolved}")

    if rows:
        AUDIT.parent.mkdir(parents=True, exist_ok=True)
        with AUDIT.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"audit written to {AUDIT} ({len(rows)} rows)")

    work.close()
    cur.close()
    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
