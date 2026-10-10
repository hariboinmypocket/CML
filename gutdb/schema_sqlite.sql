-- SQLite translation of gutdb/schema.sql (MySQL 8).
-- Table and column names, keys and semantics are preserved 1:1; only
-- engine-specific syntax changes (AUTO_INCREMENT -> AUTOINCREMENT, ENUM ->
-- TEXT + CHECK, boolean expressions -> CASE, integer division -> *1.0).
-- One deliberate addition: taxon_disease_associations.p_value / .q_value, which
-- the literature marker CSVs carry and the MySQL schema currently discards.
--
-- Five deliberate omissions, all empty in every row this mirror can load:
--   ingestion_runs.error_message  -- only gutdb/pipeline.py writes it (MySQL side);
--     this mirror's Loader.finish_run takes no error argument, so a failed load
--     leaves its row at status = 'running' rather than recording a message.
--   samples.subject_id, samples.timepoint  -- written only by
--     scripts/load_biomapai_study.py (PRJNA1125469, longitudinal), a study this
--     mirror does not load; still live in gutdb/schema.sql.
--   samples.unclassified_fraction  -- computed by pipeline.sync_gmrepo_abundances
--     and written straight to MySQL; it is not in that function's CSV export
--     fieldnames, so no value ever reaches this build.
--   sample_taxon_abundances.detection_threshold  -- dropped from the MySQL side
--     too: pipeline.py read it from a source field no ingested export carries,
--     so it was NULL in all 2.58M rows there as well.
-- Restore any of them here if the corresponding source becomes loadable.
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS ingestion_runs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    source_name TEXT NOT NULL,
    source_uri TEXT,
    started_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
    completed_at TEXT,
    status TEXT NOT NULL DEFAULT 'running'
        CHECK (status IN ('running', 'completed', 'failed')),
    rows_read INTEGER NOT NULL DEFAULT 0,
    rows_inserted INTEGER NOT NULL DEFAULT 0,
    rows_updated INTEGER NOT NULL DEFAULT 0,
    rows_skipped INTEGER NOT NULL DEFAULT 0
);

CREATE TABLE IF NOT EXISTS taxa (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    superkingdom TEXT,
    phylum TEXT,
    class_name TEXT,
    order_name TEXT,
    family TEXT,
    genus TEXT NOT NULL,
    species TEXT NOT NULL DEFAULT '',
    genus_key TEXT NOT NULL,
    species_key TEXT NOT NULL DEFAULT '',
    scientific_name TEXT NOT NULL,
    taxonomic_rank TEXT NOT NULL DEFAULT 'species'
        CHECK (taxonomic_rank IN ('species', 'genus', 'other')),
    ncbi_tax_id INTEGER,
    gram_stain TEXT,
    oxygen_requirement TEXT,
    ph_preference TEXT,
    sporulation TEXT,
    shape TEXT,
    cell_arrangement TEXT,
    mobility TEXT,
    flagella_presence TEXT,
    number_of_membranes INTEGER,
    biotic_relationship TEXT,
    habitat TEXT,
    temperature_range TEXT,
    optimal_temperature REAL,
    metabolism TEXT,
    energy_source TEXT,
    -- 1 = documented human pathogen. NULL means UNKNOWN, never "not a pathogen".
    human_pathogen INTEGER,
    energy_mode TEXT NOT NULL DEFAULT 'N/A'
        CHECK (energy_mode IN ('fermenter', 'respirator', 'mixed', 'N/A')),
    primary_food_source TEXT NOT NULL DEFAULT 'N/A',
    source_database TEXT NOT NULL DEFAULT 'MiMeDB',
    UNIQUE (genus_key, species_key)
);
CREATE INDEX IF NOT EXISTS ix_taxa_ncbi ON taxa (ncbi_tax_id);
CREATE INDEX IF NOT EXISTS ix_taxa_scientific_name ON taxa (scientific_name);

CREATE TABLE IF NOT EXISTS diseases (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    mesh_id TEXT UNIQUE,
    name TEXT NOT NULL,
    name_key TEXT NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS studies (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    project_accession TEXT NOT NULL UNIQUE,
    title TEXT,
    description TEXT,
    data_type TEXT,
    source_database TEXT NOT NULL DEFAULT 'GMrepo',
    data_quality TEXT NOT NULL DEFAULT 'unknown'
        CHECK (data_quality IN ('curated', 'qualified', 'unknown'))
);

CREATE TABLE IF NOT EXISTS phenotype_comparisons (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    study_id INTEGER NOT NULL REFERENCES studies(id) ON DELETE CASCADE,
    phenotype_a_id INTEGER NOT NULL REFERENCES diseases(id),
    phenotype_b_id INTEGER NOT NULL REFERENCES diseases(id),
    comparison_key TEXT NOT NULL UNIQUE,
    method TEXT,
    taxonomic_level TEXT NOT NULL DEFAULT 'species'
        CHECK (taxonomic_level IN ('species', 'genus', 'mixed')),
    positive_score_enriched_in_id INTEGER REFERENCES diseases(id),
    negative_score_enriched_in_id INTEGER REFERENCES diseases(id),
    notes TEXT
);
CREATE INDEX IF NOT EXISTS ix_comparison_study ON phenotype_comparisons (study_id);

CREATE TABLE IF NOT EXISTS taxon_disease_associations (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    taxon_id INTEGER NOT NULL REFERENCES taxa(id) ON DELETE CASCADE,
    disease_id INTEGER NOT NULL REFERENCES diseases(id) ON DELETE CASCADE,
    comparison_id INTEGER NOT NULL REFERENCES phenotype_comparisons(id) ON DELETE CASCADE,
    direction TEXT NOT NULL
        CHECK (direction IN ('enriched', 'depleted', 'marker', 'no_difference')),
    effect_type TEXT NOT NULL DEFAULT 'LDA',
    effect_size REAL,
    p_value REAL,
    q_value REAL,
    source_database TEXT NOT NULL DEFAULT 'GMrepo',
    UNIQUE (taxon_id, disease_id, comparison_id, effect_type)
);
CREATE INDEX IF NOT EXISTS ix_association_disease_direction
    ON taxon_disease_associations (disease_id, direction);
CREATE INDEX IF NOT EXISTS ix_association_taxon
    ON taxon_disease_associations (taxon_id);

CREATE TABLE IF NOT EXISTS samples (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    study_id INTEGER NOT NULL REFERENCES studies(id) ON DELETE CASCADE,
    disease_id INTEGER REFERENCES diseases(id),
    gmrepo_sample_id TEXT,
    run_accession TEXT NOT NULL UNIQUE,
    sex TEXT,
    age_years REAL,
    bmi REAL,
    country TEXT,
    qc_status TEXT,
    body_site TEXT
);
CREATE INDEX IF NOT EXISTS ix_sample_study ON samples (study_id);
CREATE INDEX IF NOT EXISTS ix_sample_disease ON samples (disease_id);

CREATE TABLE IF NOT EXISTS sample_taxon_abundances (
    sample_id INTEGER NOT NULL REFERENCES samples(id) ON DELETE CASCADE,
    taxon_id INTEGER NOT NULL REFERENCES taxa(id) ON DELETE CASCADE,
    relative_abundance REAL NOT NULL,
    PRIMARY KEY (sample_id, taxon_id)
);
CREATE INDEX IF NOT EXISTS ix_abundance_taxon ON sample_taxon_abundances (taxon_id);

DROP VIEW IF EXISTS v_taxon_disease_summary;
CREATE VIEW v_taxon_disease_summary AS
SELECT
    t.id AS taxon_id,
    t.scientific_name,
    d.id AS disease_id,
    d.name AS disease_name,
    d.mesh_id,
    a.direction,
    a.effect_type,
    a.effect_size,
    s.project_accession
FROM taxon_disease_associations a
JOIN taxa t ON t.id = a.taxon_id
JOIN diseases d ON d.id = a.disease_id
JOIN phenotype_comparisons c ON c.id = a.comparison_id
JOIN studies s ON s.id = c.study_id;

-- Per-sample ML readiness: covariates and feature vectors tracked separately,
-- because demographic completeness and feature availability are independent gaps.
DROP VIEW IF EXISTS v_ml_ready_samples;
CREATE VIEW v_ml_ready_samples AS
SELECT
    s.id AS sample_id,
    s.run_accession,
    st.project_accession,
    st.source_database AS study_source,
    st.data_type,
    d.name AS phenotype,
    d.mesh_id,
    CASE WHEN d.name = 'Health' THEN 1 ELSE 0 END AS is_control,
    s.sex,
    s.age_years,
    s.bmi,
    s.country,
    COALESCE(s.body_site, 'unspecified') AS body_site,
    CASE WHEN s.sex IS NOT NULL THEN 1 ELSE 0 END AS has_sex,
    CASE WHEN s.age_years IS NOT NULL THEN 1 ELSE 0 END AS has_age,
    CASE WHEN s.sex IS NOT NULL AND s.age_years IS NOT NULL THEN 1 ELSE 0 END AS has_demographics,
    COALESCE(ab.n_taxa, 0) AS n_taxa_profiled,
    CASE WHEN COALESCE(ab.n_taxa, 0) > 0 THEN 1 ELSE 0 END AS has_features,
    CASE WHEN s.sex IS NOT NULL AND s.age_years IS NOT NULL
              AND COALESCE(ab.n_taxa, 0) > 0 THEN 1 ELSE 0 END AS ml_ready
FROM samples s
JOIN studies st ON st.id = s.study_id
LEFT JOIN diseases d ON d.id = s.disease_id
LEFT JOIN (
    SELECT sample_id, COUNT(*) AS n_taxa
    FROM sample_taxon_abundances
    GROUP BY sample_id
) ab ON ab.sample_id = s.id;

-- Study-level rollup: demographic missingness here is study-level, not random.
DROP VIEW IF EXISTS v_ml_cohort_summary;
CREATE VIEW v_ml_cohort_summary AS
SELECT
    st.project_accession,
    st.source_database,
    st.data_quality,
    COUNT(*) AS n_samples,
    SUM(CASE WHEN s.sex IS NOT NULL THEN 1 ELSE 0 END) AS n_with_sex,
    SUM(CASE WHEN s.age_years IS NOT NULL THEN 1 ELSE 0 END) AS n_with_age,
    SUM(CASE WHEN s.bmi IS NOT NULL THEN 1 ELSE 0 END) AS n_with_bmi,
    SUM(CASE WHEN ab.sample_id IS NOT NULL THEN 1 ELSE 0 END) AS n_with_features,
    SUM(CASE WHEN d.name = 'Health' THEN 1 ELSE 0 END) AS n_control,
    SUM(CASE WHEN d.name <> 'Health' THEN 1 ELSE 0 END) AS n_case,
    COUNT(DISTINCT s.disease_id) AS n_phenotypes,
    CASE
        WHEN SUM(CASE WHEN s.sex IS NOT NULL THEN 1 ELSE 0 END) = 0 THEN 'no_demographics'
        WHEN SUM(CASE WHEN s.sex IS NOT NULL THEN 1 ELSE 0 END) = COUNT(*) THEN 'full_demographics'
        ELSE 'partial_demographics'
    END AS demographic_status
FROM samples s
JOIN studies st ON st.id = s.study_id
LEFT JOIN diseases d ON d.id = s.disease_id
LEFT JOIN (
    SELECT DISTINCT sample_id FROM sample_taxon_abundances
) ab ON ab.sample_id = s.id
GROUP BY st.project_accession, st.source_database, st.data_quality;

-- Evidence strength per (taxon, disease). Studies, not association rows, are the
-- unit of agreement, and agreement is a ratio rather than a boolean.
-- Evidence strength per (taxon, disease). Mirrors gutdb/schema.sql; read that
-- file's comment for why the direction vote is restricted the way it is.
--
-- In short: only case/control contrasts may vote, because "enriched in UC
-- relative to Crohn disease" does not answer a case/control question, and one
-- study casts one vote, because two independent COUNT(DISTINCT study_id)
-- expressions let a study holding both directions count on both sides. The
-- earlier version of this view did neither and reported 33 directions
-- backwards, Faecalibacterium prausnitzii in ulcerative colitis among them.
DROP VIEW IF EXISTS v_taxon_disease_evidence;
CREATE VIEW v_taxon_disease_evidence AS
WITH assoc AS (
    SELECT a.taxon_id,
           a.disease_id,
           c.study_id,
           a.direction,
           a.effect_size,
           a.source_database,
           COALESCE(pa.name = 'Health' OR pb.name = 'Health', 0) AS vs_health
    FROM taxon_disease_associations a
    JOIN phenotype_comparisons c ON c.id = a.comparison_id
    LEFT JOIN diseases pa ON pa.id = c.phenotype_a_id
    LEFT JOIN diseases pb ON pb.id = c.phenotype_b_id
),
per_study AS (
    SELECT taxon_id,
           disease_id,
           study_id,
           MAX(vs_health) AS has_case_control,
           CASE
               WHEN MAX(vs_health) = 0 THEN 'other_contrast'
               WHEN COUNT(DISTINCT CASE WHEN vs_health = 1 THEN direction END) > 1 THEN 'split'
               ELSE MAX(CASE WHEN vs_health = 1 THEN direction END)
           END AS vote
    FROM assoc
    GROUP BY taxon_id, disease_id, study_id
),
votes AS (
    SELECT taxon_id,
           disease_id,
           SUM(CASE WHEN vote = 'enriched' THEN 1 ELSE 0 END) AS n_studies_enriched,
           SUM(CASE WHEN vote = 'depleted' THEN 1 ELSE 0 END) AS n_studies_depleted,
           SUM(CASE WHEN vote = 'split' THEN 1 ELSE 0 END) AS n_studies_split,
           SUM(CASE WHEN vote = 'other_contrast' THEN 1 ELSE 0 END) AS n_studies_other_contrast,
           SUM(has_case_control) AS n_studies_case_control
    FROM per_study
    GROUP BY taxon_id, disease_id
),
totals AS (
    SELECT taxon_id,
           disease_id,
           COUNT(*) AS n_associations,
           COUNT(DISTINCT study_id) AS n_studies,
           COUNT(DISTINCT source_database) AS n_sources,
           ROUND(AVG(ABS(effect_size)), 3) AS mean_abs_effect_size,
           CASE
               WHEN MIN(vs_health) = 1 THEN 'vs_health'
               WHEN MAX(vs_health) = 0 THEN 'disease_vs_disease'
               ELSE 'mixed'
           END AS contrast_scope
    FROM assoc
    GROUP BY taxon_id, disease_id
)
SELECT t.id AS taxon_id,
       t.scientific_name,
       t.taxonomic_rank,
       t.phylum,
       d.id AS disease_id,
       d.name AS disease_name,
       d.mesh_id,
       o.n_associations,
       o.n_studies,
       o.n_sources,
       o.mean_abs_effect_size,
       o.contrast_scope,
       v.n_studies_case_control,
       v.n_studies_enriched,
       v.n_studies_depleted,
       v.n_studies_split,
       v.n_studies_other_contrast,
       CASE
           WHEN v.n_studies_enriched > v.n_studies_depleted THEN 'enriched'
           WHEN v.n_studies_depleted > v.n_studies_enriched THEN 'depleted'
           WHEN (v.n_studies_enriched + v.n_studies_depleted) = 0 THEN NULL
           ELSE 'tied'
       END AS consensus_direction,
       CASE
           WHEN (v.n_studies_enriched + v.n_studies_depleted) > 0
           THEN ROUND(
               MAX(v.n_studies_enriched, v.n_studies_depleted) * 1.0
               / (v.n_studies_enriched + v.n_studies_depleted), 3)
       END AS agreement_ratio,
       CASE WHEN v.n_studies_enriched > 0 AND v.n_studies_depleted > 0
            THEN 1 ELSE 0 END AS has_conflict
FROM totals o
JOIN votes v ON v.taxon_id = o.taxon_id AND v.disease_id = o.disease_id
JOIN taxa t ON t.id = o.taxon_id
JOIN diseases d ON d.id = o.disease_id;

-- How disease-specific is each taxon? depleted_fraction near 1 across many
-- diseases is the signature of a general dysbiosis marker, not a specific one.
-- It counts case/control associations only, for the same reason as above:
-- pooling both kinds moved this fraction by over 10 percentage points for 85 of
-- the 821 taxa with three or more associations.
DROP VIEW IF EXISTS v_taxon_specificity;
CREATE VIEW v_taxon_specificity AS
SELECT
    s.taxon_id,
    s.scientific_name,
    s.taxonomic_rank,
    s.phylum,
    s.n_diseases,
    s.n_associations,
    s.n_studies,
    s.n_enriched,
    s.n_depleted,
    s.n_other_contrast,
    CASE
        WHEN s.n_diseases = 1 THEN 'single_disease'
        WHEN s.n_diseases <= 5 THEN 'narrow'
        WHEN s.n_diseases <= 20 THEN 'broad'
        ELSE 'pan_disease'
    END AS specificity_class,
    CASE
        WHEN (s.n_enriched + s.n_depleted) > 0
        THEN ROUND(s.n_depleted * 1.0 / (s.n_enriched + s.n_depleted), 3)
    END AS depleted_fraction
FROM (
    SELECT
        a.taxon_id,
        t.scientific_name,
        t.taxonomic_rank,
        t.phylum,
        COUNT(DISTINCT a.disease_id) AS n_diseases,
        COUNT(*) AS n_associations,
        COUNT(DISTINCT c.study_id) AS n_studies,
        SUM(CASE WHEN a.direction = 'enriched'
                  AND COALESCE(pa.name = 'Health' OR pb.name = 'Health', 0) = 1
                 THEN 1 ELSE 0 END) AS n_enriched,
        SUM(CASE WHEN a.direction = 'depleted'
                  AND COALESCE(pa.name = 'Health' OR pb.name = 'Health', 0) = 1
                 THEN 1 ELSE 0 END) AS n_depleted,
        SUM(CASE WHEN COALESCE(pa.name = 'Health' OR pb.name = 'Health', 0) = 0
                 THEN 1 ELSE 0 END) AS n_other_contrast
    FROM taxon_disease_associations a
    JOIN taxa t ON t.id = a.taxon_id
    JOIN phenotype_comparisons c ON c.id = a.comparison_id
    LEFT JOIN diseases pa ON pa.id = c.phenotype_a_id
    LEFT JOIN diseases pb ON pb.id = c.phenotype_b_id
    GROUP BY a.taxon_id, t.scientific_name, t.taxonomic_rank, t.phylum
) s;


-- Rank-safe access to the abundance matrix: abundances sum to 1 WITHIN a rank,
-- not across a sample, so rank-filtered access is the only correct path.
DROP VIEW IF EXISTS v_abundance_genus;
CREATE VIEW v_abundance_genus AS
SELECT
    a.sample_id, s.run_accession, s.study_id, s.disease_id,
    a.taxon_id, t.scientific_name AS genus, t.phylum, t.superkingdom,
    a.relative_abundance
FROM sample_taxon_abundances a
JOIN taxa t ON t.id = a.taxon_id
JOIN samples s ON s.id = a.sample_id
WHERE t.taxonomic_rank = 'genus';

DROP VIEW IF EXISTS v_abundance_species;
CREATE VIEW v_abundance_species AS
SELECT
    a.sample_id, s.run_accession, s.study_id, s.disease_id,
    a.taxon_id, t.scientific_name AS species, t.genus, t.phylum, t.superkingdom,
    a.relative_abundance
FROM sample_taxon_abundances a
JOIN taxa t ON t.id = a.taxon_id
JOIN samples s ON s.id = a.sample_id
WHERE t.taxonomic_rank = 'species';

DROP VIEW IF EXISTS v_abundance_coverage;
CREATE VIEW v_abundance_coverage AS
SELECT
    a.sample_id,
    s.run_accession,
    t.taxonomic_rank,
    COUNT(*) AS n_taxa,
    ROUND(SUM(a.relative_abundance), 6) AS rank_total,
    CASE WHEN ABS(SUM(a.relative_abundance) - 1.0) <= 0.01 THEN 1 ELSE 0 END AS sums_to_one
FROM sample_taxon_abundances a
JOIN taxa t ON t.id = a.taxon_id
JOIN samples s ON s.id = a.sample_id
GROUP BY a.sample_id, s.run_accession, t.taxonomic_rank;
