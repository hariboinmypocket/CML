"""Query the SQLite gut-microbiome-disease database and export result tables."""
from __future__ import annotations

import os
import sqlite3
import numpy as np
import pandas as pd

DB = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")

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
               n_studies, n_enriched, n_depleted, specificity_class, depleted_fraction
        FROM v_taxon_specificity
        WHERE n_diseases >= 5
        ORDER BY n_diseases DESC, n_associations DESC
    """,
    # 4. Replicated (taxon, disease) pairs and the ones where studies disagree.
    "q4_replicated_evidence": """
        SELECT scientific_name, disease_name, mesh_id, n_studies, n_associations,
               n_studies_enriched, n_studies_depleted, consensus_direction,
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

# 7. Do the curated associations agree with the abundance matrix? For every
#    species-level (taxon, disease) pair that has both an association and
#    abundance data, compare its mean relative abundance in case samples of that
#    disease against Health samples. Abundances are species-rank only here, so
#    the comparison is rank-clean.
CONCORDANCE_SQL = """
WITH case_ab AS (
    SELECT sp.taxon_id, s.disease_id,
           COUNT(*) AS n_case_samples,
           AVG(sp.relative_abundance) AS mean_case
    FROM v_abundance_species sp
    JOIN samples s ON s.id = sp.sample_id
    JOIN diseases d ON d.id = s.disease_id
    WHERE d.name <> 'Health'
    GROUP BY sp.taxon_id, s.disease_id
),
health_ab AS (
    SELECT sp.taxon_id,
           COUNT(*) AS n_health_samples,
           AVG(sp.relative_abundance) AS mean_health
    FROM v_abundance_species sp
    JOIN samples s ON s.id = sp.sample_id
    JOIN diseases d ON d.id = s.disease_id
    WHERE d.name = 'Health'
    GROUP BY sp.taxon_id
)
SELECT e.scientific_name, e.disease_name, e.n_studies, e.consensus_direction,
       e.mean_abs_effect_size,
       c.n_case_samples, h.n_health_samples,
       ROUND(c.mean_case, 6)   AS mean_abundance_case,
       ROUND(h.mean_health, 6) AS mean_abundance_health,
       CASE WHEN e.consensus_direction = 'enriched' AND c.mean_case > h.mean_health THEN 1
            WHEN e.consensus_direction = 'depleted' AND c.mean_case < h.mean_health THEN 1
            WHEN e.consensus_direction IN ('enriched', 'depleted') THEN 0
       END AS abundance_agrees
FROM v_taxon_disease_evidence e
JOIN case_ab  c ON c.taxon_id = e.taxon_id AND c.disease_id = e.disease_id
JOIN health_ab h ON h.taxon_id = e.taxon_id
WHERE e.taxonomic_rank = 'species'
  AND e.consensus_direction IN ('enriched', 'depleted')
  AND c.n_case_samples >= 20 AND h.n_health_samples >= 20
"""


def main() -> dict[str, pd.DataFrame]:
    con = sqlite3.connect(DB)
    results = {name: pd.read_sql_query(sql, con) for name, sql in QUERIES.items()}
    conc = pd.read_sql_query(CONCORDANCE_SQL, con)
    # log2 fold change computed here rather than in SQL: SQLite's LOG() is only
    # present in builds compiled with the math extension.
    pseudo = 1e-6
    conc["log2_fold_change"] = np.log2(
        (conc["mean_abundance_case"] + pseudo) / (conc["mean_abundance_health"] + pseudo)
    ).round(3)
    results["q7_association_vs_abundance"] = conc.reindex(
        conc["log2_fold_change"].abs().sort_values(ascending=False).index
    ).reset_index(drop=True)
    con.close()
    for name, frame in results.items():
        frame.to_csv(f"{name}.csv", index=False)
        print(f"{name}: {len(frame)} rows -> {name}.csv")
    return results


if __name__ == "__main__":
    main()
