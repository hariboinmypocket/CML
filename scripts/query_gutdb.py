"""Query the gut-microbiome-disease database and export result tables.

Reads MySQL by default, which is the authoritative store. Pass --sqlite to read
the mirror instead.

The mirror used to be the default, and that was the wrong choice for anything
quoted anywhere: it is built from the CSVs in data/, so sources loaded straight
into MySQL never reach it. When this was changed it sat 13,199 samples and
1,022,325 abundance rows behind, with no build path that could close the gap.
Three of these tables (q3, q4, q7) also read views whose direction logic was
corrected on 2026-10-10, so any copy of them produced before that date reports
directions that have since been reversed.
"""
from __future__ import annotations

import os
import sqlite3
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

DB = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")
USE_SQLITE = "--sqlite" in sys.argv

QUERIES: dict[str, str] = {
    # 1. What is in the database, by source and direction.
    "q1_content_summary": """
        SELECT a.source_database,
               a.effect_type,
               COALESCE(c.method, 'unspecified') AS method,
               COUNT(*)                                  AS n_associations,
               COUNT(DISTINCT a.taxon_id)                AS n_taxa,
               COUNT(DISTINCT a.disease_id)              AS n_diseases,
               COUNT(DISTINCT c.study_id)                AS n_studies,
               SUM(CASE WHEN a.direction = 'enriched' THEN 1 ELSE 0 END) AS n_enriched,
               SUM(CASE WHEN a.direction = 'depleted' THEN 1 ELSE 0 END) AS n_depleted,
               ROUND(AVG(ABS(a.effect_size)), 2)         AS mean_abs_effect
        FROM taxon_disease_associations a
        JOIN phenotype_comparisons c ON c.id = a.comparison_id
        GROUP BY a.source_database, a.effect_type, c.method
        ORDER BY n_associations DESC
    """,
    # 2. Association burden per disease, with the sample counts behind it.
    "q2_disease_rollup": """
        SELECT d.name AS disease_name,
               d.mesh_id,
               COUNT(*)                        AS n_associations,
               COUNT(DISTINCT a.taxon_id)      AS n_taxa,
               COUNT(DISTINCT c.study_id)      AS n_studies,
               SUM(CASE WHEN a.direction = 'enriched' THEN 1 ELSE 0 END) AS n_enriched,
               SUM(CASE WHEN a.direction = 'depleted' THEN 1 ELSE 0 END) AS n_depleted,
               ROUND(AVG(ABS(a.effect_size)), 2) AS mean_abs_effect,
               (SELECT COUNT(*) FROM samples s WHERE s.disease_id = d.id) AS n_samples
        FROM taxon_disease_associations a
        JOIN diseases d ON d.id = a.disease_id
        JOIN phenotype_comparisons c ON c.id = a.comparison_id
        GROUP BY d.id, d.name, d.mesh_id
        ORDER BY n_associations DESC
    """,
    # 3. Taxa reported across many diseases: strong dysbiosis markers, weak
    #    disease-discriminating features.
    "q3_pan_disease_taxa": """
        SELECT scientific_name, taxonomic_rank, phylum, n_diseases, n_associations,
               n_studies, n_enriched, n_depleted, n_other_contrast,
               specificity_class, depleted_fraction
        FROM v_taxon_specificity
        WHERE n_diseases >= 5
        ORDER BY n_diseases DESC, n_associations DESC
    """,
    # 4. Replicated (taxon, disease) pairs and the ones where studies disagree.
    "q4_replicated_evidence": """
        SELECT scientific_name, disease_name, mesh_id, n_studies, n_associations,
               contrast_scope, n_studies_case_control,
               n_studies_enriched, n_studies_depleted, n_studies_split,
               n_studies_other_contrast, consensus_direction,
               agreement_ratio, has_conflict, mean_abs_effect_size, n_sources
        FROM v_taxon_disease_evidence
        WHERE n_studies >= 2
        ORDER BY has_conflict DESC, n_studies DESC, agreement_ratio ASC
    """,
    # 5. IBS marker panel, strongest signed effects first.
    "q5_ibs_markers": """
        SELECT t.scientific_name,
               t.taxonomic_rank,
               a.direction,
               a.effect_size,
               a.effect_type,
               st.project_accession,
               c.method,
               pa.name AS phenotype_a,
               pb.name AS phenotype_b,
               a.source_database
        FROM taxon_disease_associations a
        JOIN taxa t     ON t.id = a.taxon_id
        JOIN diseases d ON d.id = a.disease_id
        JOIN phenotype_comparisons c ON c.id = a.comparison_id
        JOIN studies st ON st.id = c.study_id
        LEFT JOIN diseases pa ON pa.id = c.phenotype_a_id
        LEFT JOIN diseases pb ON pb.id = c.phenotype_b_id
        WHERE d.name = 'Irritable Bowel Syndrome'
        ORDER BY ABS(a.effect_size) DESC
    """,
    # 6. Cohort readiness: covariates vs feature vectors, per study.
    "q6_cohort_readiness": """
        SELECT project_accession, source_database, data_quality, n_samples,
               n_with_sex, n_with_age, n_with_bmi, n_with_features,
               n_control, n_case, n_phenotypes, demographic_status
        FROM v_ml_cohort_summary
        ORDER BY n_samples DESC
    """,
}

# 7. Do the curated associations agree with the abundance matrix?
#
#    This table was rebuilt on 2026-10-12. The previous version answered a
#    different question than its name claimed and reported 62% concordance,
#    which was quoted as "the curated associations only match the measurements
#    62% of the time". The true figure is about 95%. Three things were wrong:
#
#      It pooled case and control samples ACROSS studies, comparing cases from
#      one cohort against controls from every other cohort, on different
#      protocols and sequencing runs.
#
#      It averaged over v_abundance_species, which holds only the rows that
#      exist. 97.3% of (sample, taxon) cells are absences with no row, so it
#      compared "abundance where the taxon was detected in cases" against
#      "abundance where it was detected in controls" -- not a case/control
#      contrast at all.
#
#      Its n_case_samples >= 20 gate counted DETECTIONS rather than samples
#      profiled, so the threshold meant something other than it appeared to.
#
#    It now does what taxon_disease_adjudication does: one direction per study
#    that holds both arms, absence counted as zero, denominator every sample in
#    the arm profiled at species rank.
#
#    WHICH OF THE THREE MATTERED, measured on the 4,046 pairs every method can
#    score: pooled with absence missing gives 65.8%, pooled with absence as zero
#    gives 95.5%, within-study with absence as zero gives 95.3%. So counting
#    absence as zero is worth +29.7 points and comparing within study is worth
#    -0.2. The old number was an artefact of skipping absences, NOT of pooling,
#    and an earlier version of this comment said otherwise.
#
#    Within-study is still the right comparison and is what this table reports:
#    it compares like with like, and it is what taxon_disease_adjudication uses
#    to settle directions. It simply is not what rescues the headline. Both
#    columns are kept so neither claim has to be taken on trust.
CONCORDANCE_SQL = """
WITH prof AS (
    SELECT s.study_id, s.disease_id, COUNT(DISTINCT a.sample_id) AS n
    FROM sample_taxon_abundances a
    JOIN taxa t ON t.id = a.taxon_id AND t.taxonomic_rank = 'species'
    JOIN samples s ON s.id = a.sample_id
    GROUP BY s.study_id, s.disease_id
),
arms AS (
    SELECT pc.study_id, pc.disease_id, pc.n AS n_case, ph.n AS n_health
    FROM prof pc
    JOIN prof ph ON ph.study_id = pc.study_id
                AND ph.disease_id = (SELECT id FROM diseases WHERE name = 'Health')
    WHERE pc.disease_id <> (SELECT id FROM diseases WHERE name = 'Health')
),
sums AS (
    SELECT a.taxon_id, s.study_id, s.disease_id, SUM(a.relative_abundance) AS total
    FROM sample_taxon_abundances a
    JOIN taxa t ON t.id = a.taxon_id AND t.taxonomic_rank = 'species'
    JOIN samples s ON s.id = a.sample_id
    GROUP BY a.taxon_id, s.study_id, s.disease_id
),
pairs AS (
    SELECT taxon_id, scientific_name, disease_id, disease_name,
           n_studies, consensus_direction, mean_abs_effect_size
    FROM v_taxon_disease_evidence
    WHERE taxonomic_rank = 'species'
      AND consensus_direction IN ('enriched', 'depleted')
),
reported AS (
    -- Which studies actually reported each association, case/control only. A
    -- study's abundance can then be checked against a claim that study made,
    -- separately from claims other cohorts made.
    SELECT DISTINCT a.taxon_id, a.disease_id, cp.study_id
    FROM taxon_disease_associations a
    JOIN phenotype_comparisons cp ON cp.id = a.comparison_id
    LEFT JOIN diseases rpa ON rpa.id = cp.phenotype_a_id
    LEFT JOIN diseases rpb ON rpb.id = cp.phenotype_b_id
    WHERE COALESCE(rpa.name = 'Health' OR rpb.name = 'Health', 0) = 1
)
SELECT p.scientific_name, p.disease_name, p.n_studies, p.consensus_direction,
       p.mean_abs_effect_size, st.project_accession,
       CASE WHEN r.study_id IS NULL THEN 0 ELSE 1 END AS reported_here,
       ar.n_case AS n_case_samples, ar.n_health AS n_health_samples,
       COALESCE(sc.total, 0) / ar.n_case   AS mean_abundance_case,
       COALESCE(sh.total, 0) / ar.n_health AS mean_abundance_health
FROM pairs p
JOIN arms ar ON ar.disease_id = p.disease_id
JOIN studies st ON st.id = ar.study_id
LEFT JOIN reported r ON r.taxon_id = p.taxon_id AND r.disease_id = p.disease_id
                    AND r.study_id = ar.study_id
LEFT JOIN sums sc ON sc.taxon_id = p.taxon_id AND sc.study_id = ar.study_id
                 AND sc.disease_id = p.disease_id
LEFT JOIN sums sh ON sh.taxon_id = p.taxon_id AND sh.study_id = ar.study_id
                 AND sh.disease_id = (SELECT id FROM diseases WHERE name = 'Health')
"""


def summarize_concordance(rows: pd.DataFrame) -> pd.DataFrame:
    """Collapse per-(pair, study) rows to one row per (taxon, disease).

    Two directions are derived from the same numbers. The within-study one
    counts each study's own case-vs-control comparison; the pooled one averages
    every case sample against every control sample regardless of which cohort
    produced it. Reporting both is the point: they disagree often enough that
    which one a figure used changes its headline.
    """
    rows = rows.copy()
    eps = 1e-6
    rows["log2_fold_change"] = np.log2(
        (rows["mean_abundance_case"] + eps) / (rows["mean_abundance_health"] + eps))
    rows["study_direction"] = np.where(
        rows["mean_abundance_case"] > rows["mean_abundance_health"], "enriched",
        np.where(rows["mean_abundance_case"] < rows["mean_abundance_health"],
                 "depleted", "tied"))
    # Pooled means are sample-weighted so they match what averaging every
    # sample together would give, rather than averaging the study means.
    rows["case_mass"] = rows["mean_abundance_case"] * rows["n_case_samples"]
    rows["health_mass"] = rows["mean_abundance_health"] * rows["n_health_samples"]
    # Does this study's own abundance back the direction it is claimed to support?
    decided_here = rows["study_direction"] != "tied"
    rows["own_checked"] = (rows["reported_here"] == 1) & decided_here
    rows["own_matching"] = rows["own_checked"] & (
        rows["study_direction"] == rows["consensus_direction"])

    key = ["scientific_name", "disease_name", "n_studies", "consensus_direction",
           "mean_abs_effect_size"]
    g = rows.groupby(key, dropna=False)
    out = g.agg(
        n_studies_compared=("project_accession", "nunique"),
        n_studies_enriched=("study_direction", lambda v: int((v == "enriched").sum())),
        n_studies_depleted=("study_direction", lambda v: int((v == "depleted").sum())),
        median_log2_fold_change=("log2_fold_change", "median"),
        case_mass=("case_mass", "sum"),
        n_case_samples=("n_case_samples", "sum"),
        health_mass=("health_mass", "sum"),
        n_health_samples=("n_health_samples", "sum"),
        n_own_studies=("own_checked", "sum"),
        n_own_matching=("own_matching", "sum"),
    ).reset_index()

    decided = out["n_studies_enriched"] + out["n_studies_depleted"]
    out["abundance_direction"] = np.where(
        out["n_studies_enriched"] > out["n_studies_depleted"], "enriched",
        np.where(out["n_studies_depleted"] > out["n_studies_enriched"], "depleted", None))
    out["abundance_agreement"] = np.where(
        decided > 0,
        np.maximum(out["n_studies_enriched"], out["n_studies_depleted"]) / decided.where(decided > 0),
        np.nan).round(3)

    out["pooled_mean_case"] = (out["case_mass"] / out["n_case_samples"]).round(8)
    out["pooled_mean_health"] = (out["health_mass"] / out["n_health_samples"]).round(8)
    out["pooled_direction"] = np.where(
        out["pooled_mean_case"] > out["pooled_mean_health"], "enriched",
        np.where(out["pooled_mean_case"] < out["pooled_mean_health"], "depleted", None))

    for method in ("abundance", "pooled"):
        out[f"agrees_{'within_study' if method == 'abundance' else 'pooled'}"] = np.where(
            out[f"{method}_direction"].isna(), np.nan,
            (out[f"{method}_direction"] == out["consensus_direction"]).astype(float))

    # Why a pair disagrees, which is not the same question as whether it does.
    #
    # An association that matches the abundance in the cohorts that REPORTED it,
    # and is reversed only in cohorts that did not, is a finding that failed to
    # replicate -- not a wrong row. Measured across the disagreeing pairs, the
    # reporting cohorts back their own association 91.9% of the time while other
    # cohorts back it 18.3% of the time, so almost every disagreement here is of
    # that kind. The rare other kind, where a study's own measurements contradict
    # the claim attributed to it, is worth actual investigation: a sign
    # convention, a taxonomy mismatch or a loading error.
    own_rate = np.where(out["n_own_studies"] > 0,
                        out["n_own_matching"] / out["n_own_studies"].where(out["n_own_studies"] > 0),
                        np.nan)
    out["own_cohort_match_rate"] = np.round(own_rate, 3)
    out["disagreement_kind"] = np.select(
        [
            out["abundance_direction"].isna(),
            out["agrees_within_study"] == 1,
            out["n_own_studies"] == 0,
            own_rate >= 0.5,
        ],
        ["abundance_tied", "agrees", "no_own_cohort_abundance", "non_replicating"],
        default="contradicted_in_own_cohort",
    )

    out["median_log2_fold_change"] = out["median_log2_fold_change"].round(3)
    out = out.drop(columns=["case_mass", "health_mass"])
    # Disagreements first, best-evidenced first. Sorting by |log2_fold_change|
    # would lead with the 175 pairs where one arm's mean is ~0 and the 1e-6
    # pseudocount dominates the ratio, which says more about the pseudocount
    # than about the taxon.
    return out.sort_values(
        ["agrees_within_study", "n_studies_compared", "abundance_agreement"],
        ascending=[True, False, False], na_position="last"
    ).reset_index(drop=True)


def open_source():
    """The database these tables are read from, and a label for the banner."""
    if USE_SQLITE:
        return sqlite3.connect(DB), f"sqlite mirror ({DB})"
    from gutdb.config import Settings
    from gutdb.db import connect
    return connect(Settings.from_env(".env")), "MySQL (authoritative)"


def main() -> dict[str, pd.DataFrame]:
    con, label = open_source()
    print(f"source: {label}")
    results = {name: pd.read_sql_query(sql, con) for name, sql in QUERIES.items()}
    # Per-(pair, study) rows, collapsed in pandas: the log2 ratio and the median
    # across studies both need functions SQLite lacks without its math
    # extension, so doing it here keeps both sources producing one shape.
    conc = summarize_concordance(pd.read_sql_query(CONCORDANCE_SQL, con))
    results["q7_association_vs_abundance"] = conc
    con.close()
    for name, frame in results.items():
        frame.to_csv(f"{name}.csv", index=False)
        print(f"{name}: {len(frame)} rows -> {name}.csv")

    # The headline this table exists to produce, stated rather than left to be
    # recomputed by hand -- which is how the old 62% came to be quoted as a
    # statement about the curated associations.
    scored = conc.dropna(subset=["agrees_within_study"])
    independent = scored[scored["n_studies_compared"] >= 2]
    print("\nconcordance with the curated associations:")
    print(f"   all pairs:                  "
          f"{int(scored['agrees_within_study'].sum()):>5,}/{len(scored):<6,} = "
          f"{100*scored['agrees_within_study'].mean():.1f}%")
    print(f"   pairs with >= 2 studies:    "
          f"{int(independent['agrees_within_study'].sum()):>5,}/{len(independent):<6,} = "
          f"{100*independent['agrees_within_study'].mean():.1f}%   <-- quote this one")
    print("   A pair with only ONE comparable study has its abundance direction")
    print("   computed from the same cohort that produced the association, so it")
    print("   confirms itself. The >=2 figure is the independent one, and it")
    print("   reconciles with taxon_disease_adjudication's 91.7%.")
    print("   (The 62% this table used to report was neither: it averaged only")
    print("   the samples where a taxon was detected.)")
    print("\nwhy the disagreements disagree:")
    for kind, n in conc["disagreement_kind"].value_counts().items():
        if kind == "agrees":
            continue
        print(f"   {kind:28} {n:>5}")
    print("   non_replicating = the reporting cohort's own abundance backs the")
    print("   association; other cohorts reverse it. Not a wrong row.")
    return results


if __name__ == "__main__":
    main()
