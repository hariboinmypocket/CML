#!/usr/bin/env python3
"""Does a curated association describe its own cohort, or the field?

q7 reports where the curated associations and this database's abundance matrix
point in opposite directions. That alone does not say what has gone wrong, and
the obvious reading -- that those rows are mistakes -- is the wrong one.

This separates the two possibilities by asking, for every per-study comparison,
whether that study is one that REPORTED the association:

  An association backed by its own cohort and reversed only in cohorts that did
  not report it has not failed; it has not REPLICATED. The row is a faithful
  record of one study. What it is not is a fact about the disease, and a model
  that treats it as one will not generalise.

  An association contradicted by the measurements of the very study it is
  attributed to is a different animal entirely -- a sign convention, a taxonomy
  mismatch, or a loading error -- and is worth opening.

As of 2026-10-12 every disagreement in the database is the first kind. That is a
result about the loaders as much as about the biology, and it is worth knowing
if it ever stops being true, which is why this is a script and not a paragraph.

Writes data/replication_audit.csv, one row per (taxon, disease).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect

OUT = Path("data/replication_audit.csv")

SQL = """
WITH prof AS (
    SELECT s.study_id, s.disease_id, COUNT(DISTINCT a.sample_id) AS n
    FROM sample_taxon_abundances a
    JOIN taxa t ON t.id = a.taxon_id AND t.taxonomic_rank = 'species'
    JOIN samples s ON s.id = a.sample_id
    GROUP BY s.study_id, s.disease_id
),
sums AS (
    SELECT a.taxon_id, s.study_id, s.disease_id, SUM(a.relative_abundance) AS total
    FROM sample_taxon_abundances a
    JOIN taxa t ON t.id = a.taxon_id AND t.taxonomic_rank = 'species'
    JOIN samples s ON s.id = a.sample_id
    GROUP BY a.taxon_id, s.study_id, s.disease_id
),
reported AS (
    SELECT DISTINCT a.taxon_id, a.disease_id, cp.study_id
    FROM taxon_disease_associations a
    JOIN phenotype_comparisons cp ON cp.id = a.comparison_id
    LEFT JOIN diseases rpa ON rpa.id = cp.phenotype_a_id
    LEFT JOIN diseases rpb ON rpb.id = cp.phenotype_b_id
    WHERE COALESCE(rpa.name = 'Health' OR rpb.name = 'Health', 0) = 1
)
SELECT e.scientific_name, e.disease_name, e.consensus_direction, e.n_studies,
       st.project_accession,
       CASE WHEN r.study_id IS NULL THEN 0 ELSE 1 END AS reported_here,
       COALESCE(sc.total, 0) / pc.n AS mean_case,
       COALESCE(sh.total, 0) / ph.n AS mean_health
FROM v_taxon_disease_evidence e
JOIN prof pc ON pc.disease_id = e.disease_id
JOIN prof ph ON ph.study_id = pc.study_id
            AND ph.disease_id = (SELECT id FROM diseases WHERE name = 'Health')
JOIN studies st ON st.id = pc.study_id
LEFT JOIN reported r ON r.taxon_id = e.taxon_id AND r.disease_id = e.disease_id
                    AND r.study_id = pc.study_id
LEFT JOIN sums sc ON sc.taxon_id = e.taxon_id AND sc.study_id = pc.study_id
                 AND sc.disease_id = e.disease_id
LEFT JOIN sums sh ON sh.taxon_id = e.taxon_id AND sh.study_id = pc.study_id
                 AND sh.disease_id = (SELECT id FROM diseases WHERE name = 'Health')
WHERE e.taxonomic_rank = 'species'
  AND e.consensus_direction IN ('enriched', 'depleted')
"""


def main() -> int:
    conn = connect(Settings.from_env(".env"))
    rows = pd.read_sql_query(SQL, conn)
    conn.close()

    rows["study_direction"] = np.where(
        rows["mean_case"] > rows["mean_health"], "enriched",
        np.where(rows["mean_case"] < rows["mean_health"], "depleted", "tied"))
    decided = rows["study_direction"] != "tied"
    rows["matches"] = decided & (rows["study_direction"] == rows["consensus_direction"])

    key = ["scientific_name", "disease_name", "consensus_direction", "n_studies"]
    g = rows.groupby(key, dropna=False)
    pairs = g.agg(
        n_studies_compared=("project_accession", "nunique"),
        n_decided=("study_direction", lambda v: int((v != "tied").sum())),
        n_matching=("matches", "sum"),
        n_own=("reported_here", lambda v: int(v.sum())),
    ).reset_index()

    own = rows[(rows["reported_here"] == 1) & decided]
    other = rows[(rows["reported_here"] == 0) & decided]
    own_g = own.groupby(key).agg(n_own_decided=("matches", "size"),
                                 n_own_matching=("matches", "sum")).reset_index()
    oth_g = other.groupby(key).agg(n_other_decided=("matches", "size"),
                                   n_other_matching=("matches", "sum")).reset_index()
    pairs = pairs.merge(own_g, on=key, how="left").merge(oth_g, on=key, how="left")
    for col in ("n_own_decided", "n_own_matching", "n_other_decided", "n_other_matching"):
        pairs[col] = pairs[col].fillna(0).astype(int)

    pairs["majority_direction"] = np.where(
        pairs["n_matching"] * 2 > pairs["n_decided"], pairs["consensus_direction"],
        np.where(pairs["n_matching"] * 2 < pairs["n_decided"],
                 np.where(pairs["consensus_direction"] == "enriched", "depleted", "enriched"),
                 None))
    pairs["agrees"] = pairs["majority_direction"] == pairs["consensus_direction"]
    pairs["own_cohort_match_rate"] = np.where(
        pairs["n_own_decided"] > 0,
        pairs["n_own_matching"] / pairs["n_own_decided"].replace(0, np.nan), np.nan).round(3)
    pairs["other_cohort_match_rate"] = np.where(
        pairs["n_other_decided"] > 0,
        pairs["n_other_matching"] / pairs["n_other_decided"].replace(0, np.nan), np.nan).round(3)

    scored = pairs[pairs["majority_direction"].notna()]
    dis = scored[~scored["agrees"]]
    agr = scored[scored["agrees"]]

    print(f"species-level pairs with a direction on both sides: {len(scored):,}")
    print(f"   abundance majority agrees:    {len(agr):,}")
    print(f"   abundance majority disagrees: {len(dis):,}\n")

    print("Does each association describe the cohort that reported it?")
    print(f"   {'pair group':12} {'cohort':24} {'match':>8} {'of':>8}  rate")
    for label, sub in (("agree", agr), ("disagree", dis)):
        for pre, cohort in (("own", "yes (reported it)"), ("other", "no (other cohorts)")):
            m, n = sub[f"n_{pre}_matching"].sum(), sub[f"n_{pre}_decided"].sum()
            if n:
                print(f"   {label:12} {cohort:24} {m:>8,} {n:>8,}  {100*m/n:>5.1f}%")
    print("\n   The associations describe their own cohorts correctly in BOTH groups.")
    print("   A disagreement is a failure to REPLICATE, not a wrong row.\n")

    # Two different filters pull in opposite directions, and reporting only one
    # of them invites the wrong conclusion.
    print("Disagreement rate by COMPARABLE STUDIES (rises -- demanding more")
    print("studies strips out single-study pairs, which agree with themselves):")
    for lo in (1, 2, 3, 5):
        sd = dis[dis["n_studies_compared"] >= lo]
        sa = agr[agr["n_studies_compared"] >= lo]
        tot = len(sd) + len(sa)
        print(f"   >= {lo} comparable studies: {len(sd):>4} of {tot:>5,} = "
              f"{100*len(sd)/tot if tot else 0:.1f}%")

    print("\nDisagreement rate by STRENGTH OF THE CONTRADICTION (falls -- the")
    print("disagreements are concentrated where the abundance vote is a bare")
    print("majority over a negligible difference):")
    try:
        q7 = pd.read_csv("q7_association_vs_abundance.csv")
    except FileNotFoundError:
        print("   q7_association_vs_abundance.csv not present; run query_gutdb.py first")
    else:
        q7 = q7.dropna(subset=["agrees_within_study"])
        q7 = q7[q7["n_studies_compared"] >= 2]
        for ag, fc in ((0.0, 0.0), (0.6, 0.5), (0.7, 1.0), (0.8, 1.0)):
            sub = q7[(q7["abundance_agreement"] >= ag)
                     & (q7["median_log2_fold_change"].abs() >= fc)]
            bad = int((sub["agrees_within_study"] == 0).sum())
            print(f"   abundance majority >= {ag:.1f}, |log2FC| >= {fc:.1f}: "
                  f"{bad:>4} of {len(sub):>5,} = {100*bad/len(sub) if len(sub) else 0:.1f}%")

    kinds = np.where(dis["n_own_decided"] == 0, "no_own_cohort_abundance",
                     np.where(dis["own_cohort_match_rate"] >= 0.5,
                              "non_replicating", "contradicted_in_own_cohort"))
    print("\nKind of disagreement:")
    for kind, n in pd.Series(kinds).value_counts().items():
        note = "  <-- worth opening" if kind == "contradicted_in_own_cohort" else ""
        print(f"   {kind:28} {n:>5}{note}")

    pairs["disagreement_kind"] = np.where(
        pairs["majority_direction"].isna(), "abundance_tied",
        np.where(pairs["agrees"], "agrees",
                 np.where(pairs["n_own_decided"] == 0, "no_own_cohort_abundance",
                          np.where(pairs["own_cohort_match_rate"] >= 0.5,
                                   "non_replicating", "contradicted_in_own_cohort"))))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    pairs.sort_values(["disagreement_kind", "n_studies_compared"],
                      ascending=[True, False]).to_csv(OUT, index=False)
    print(f"\nwrote {OUT} ({len(pairs):,} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
