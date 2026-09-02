CREATE TABLE IF NOT EXISTS ingestion_runs (
    id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    source_name VARCHAR(100) NOT NULL,
    source_uri TEXT NULL,
    started_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    completed_at TIMESTAMP NULL,
    status ENUM('running', 'completed', 'failed') NOT NULL DEFAULT 'running',
    rows_read BIGINT UNSIGNED NOT NULL DEFAULT 0,
    rows_inserted BIGINT UNSIGNED NOT NULL DEFAULT 0,
    rows_updated BIGINT UNSIGNED NOT NULL DEFAULT 0,
    rows_skipped BIGINT UNSIGNED NOT NULL DEFAULT 0,
    error_message TEXT NULL,
    PRIMARY KEY (id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE IF NOT EXISTS taxa (
    id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    -- superkingdom replaces the former `kingdom` column, which mixed two
    -- vocabularies (Eubacteria/Fungi alongside NCBI's newer Bacillati/
    -- Pseudomonadati) and carried no information superkingdom lacked.
    superkingdom VARCHAR(50) NULL,
    phylum VARCHAR(150) NULL,
    class_name VARCHAR(150) NULL,
    order_name VARCHAR(150) NULL,
    family VARCHAR(150) NULL,
    genus VARCHAR(150) NOT NULL,
    species VARCHAR(255) NOT NULL DEFAULT '',
    genus_key VARCHAR(150) NOT NULL,
    species_key VARCHAR(255) NOT NULL DEFAULT '',
    scientific_name VARCHAR(420) NOT NULL,
    taxonomic_rank ENUM('species', 'genus', 'other') NOT NULL DEFAULT 'species',
    ncbi_tax_id BIGINT NULL,
    gram_stain VARCHAR(50) NULL,
    oxygen_requirement VARCHAR(100) NULL,
    ph_preference VARCHAR(100) NULL,
    sporulation VARCHAR(100) NULL,
    shape VARCHAR(100) NULL,
    cell_arrangement VARCHAR(150) NULL,
    mobility VARCHAR(10) NULL,
    flagella_presence VARCHAR(10) NULL,
    number_of_membranes TINYINT UNSIGNED NULL,
    biotic_relationship VARCHAR(100) NULL,
    habitat VARCHAR(150) NULL,
    temperature_range VARCHAR(50) NULL,
    optimal_temperature DOUBLE NULL,
    metabolism VARCHAR(255) NULL,
    energy_source VARCHAR(150) NULL,
    -- 1 = documented human pathogen in MiMeDB. NULL means UNKNOWN, not "not a
    -- pathogen": MiMeDB only ever records the positive flag, so absence carries
    -- no negative evidence and must not be read as 0.
    human_pathogen TINYINT UNSIGNED NULL,
    energy_mode ENUM('fermenter', 'respirator', 'mixed', 'N/A') NOT NULL DEFAULT 'N/A',
    primary_food_source VARCHAR(100) NOT NULL DEFAULT 'N/A',
    source_database VARCHAR(100) NOT NULL DEFAULT 'MiMeDB',
    PRIMARY KEY (id),
    UNIQUE KEY uq_taxa_key (genus_key, species_key),
    KEY ix_taxa_ncbi (ncbi_tax_id),
    KEY ix_taxa_scientific_name (scientific_name)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE IF NOT EXISTS diseases (
    id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    mesh_id VARCHAR(32) NULL,
    name VARCHAR(255) NOT NULL,
    name_key VARCHAR(255) NOT NULL,
    PRIMARY KEY (id),
    UNIQUE KEY uq_disease_name (name_key),
    UNIQUE KEY uq_disease_mesh (mesh_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE IF NOT EXISTS studies (
    id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    project_accession VARCHAR(64) NOT NULL,
    title TEXT NULL,
    description TEXT NULL,
    data_type VARCHAR(100) NULL,
    source_database VARCHAR(100) NOT NULL DEFAULT 'GMrepo',
    data_quality ENUM('curated', 'qualified', 'unknown') NOT NULL DEFAULT 'unknown',
    PRIMARY KEY (id),
    UNIQUE KEY uq_study_accession (project_accession)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE IF NOT EXISTS phenotype_comparisons (
    id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    study_id BIGINT UNSIGNED NOT NULL,
    phenotype_a_id BIGINT UNSIGNED NOT NULL,
    phenotype_b_id BIGINT UNSIGNED NOT NULL,
    comparison_key VARCHAR(600) NOT NULL,
    method VARCHAR(100) NULL,
    taxonomic_level ENUM('species', 'genus', 'mixed') NOT NULL DEFAULT 'species',
    positive_score_enriched_in_id BIGINT UNSIGNED NULL,
    negative_score_enriched_in_id BIGINT UNSIGNED NULL,
    notes TEXT NULL,
    PRIMARY KEY (id),
    UNIQUE KEY uq_comparison (comparison_key),
    KEY ix_comparison_study (study_id),
    CONSTRAINT fk_comparison_study FOREIGN KEY (study_id) REFERENCES studies(id)
        ON UPDATE CASCADE ON DELETE CASCADE,
    CONSTRAINT fk_comparison_a FOREIGN KEY (phenotype_a_id) REFERENCES diseases(id)
        ON UPDATE CASCADE ON DELETE RESTRICT,
    CONSTRAINT fk_comparison_b FOREIGN KEY (phenotype_b_id) REFERENCES diseases(id)
        ON UPDATE CASCADE ON DELETE RESTRICT,
    CONSTRAINT fk_comparison_positive FOREIGN KEY (positive_score_enriched_in_id) REFERENCES diseases(id)
        ON UPDATE CASCADE ON DELETE SET NULL,
    CONSTRAINT fk_comparison_negative FOREIGN KEY (negative_score_enriched_in_id) REFERENCES diseases(id)
        ON UPDATE CASCADE ON DELETE SET NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE IF NOT EXISTS taxon_disease_associations (
    id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    taxon_id BIGINT UNSIGNED NOT NULL,
    disease_id BIGINT UNSIGNED NOT NULL,
    comparison_id BIGINT UNSIGNED NOT NULL,
    direction ENUM('enriched', 'depleted', 'marker', 'no_difference') NOT NULL,
    effect_type VARCHAR(50) NOT NULL DEFAULT 'LDA',
    effect_size DOUBLE NULL,
    source_database VARCHAR(100) NOT NULL DEFAULT 'GMrepo',
    PRIMARY KEY (id),
    UNIQUE KEY uq_taxon_disease_evidence (taxon_id, disease_id, comparison_id, effect_type),
    KEY ix_association_disease_direction (disease_id, direction),
    KEY ix_association_taxon (taxon_id),
    CONSTRAINT fk_association_taxon FOREIGN KEY (taxon_id) REFERENCES taxa(id)
        ON UPDATE CASCADE ON DELETE CASCADE,
    CONSTRAINT fk_association_disease FOREIGN KEY (disease_id) REFERENCES diseases(id)
        ON UPDATE CASCADE ON DELETE CASCADE,
    CONSTRAINT fk_association_comparison FOREIGN KEY (comparison_id) REFERENCES phenotype_comparisons(id)
        ON UPDATE CASCADE ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE IF NOT EXISTS samples (
    id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    study_id BIGINT UNSIGNED NOT NULL,
    disease_id BIGINT UNSIGNED NULL,
    gmrepo_sample_id VARCHAR(128) NULL,
    run_accession VARCHAR(128) NOT NULL,
    -- Repeated measures: several runs can come from the same person. Without a
    -- subject key these look like independent observations, which inflates the
    -- effective sample size and leaks one individual across train/test splits.
    subject_id VARCHAR(128) NULL,
    timepoint VARCHAR(32) NULL,
    sex VARCHAR(50) NULL,
    age_years DOUBLE NULL,
    bmi DOUBLE NULL,
    country VARCHAR(100) NULL,
    qc_status VARCHAR(100) NULL,
    body_site VARCHAR(100) NULL,
    -- share of reads GMrepo could not assign to a taxon. Recorded rather than
    -- discarded: it is a per-sample quality signal, and the classified taxa are
    -- renormalized without it so relative_abundance still sums to 1.
    unclassified_fraction DOUBLE NULL,
    PRIMARY KEY (id),
    UNIQUE KEY uq_sample_run (run_accession),
    KEY ix_sample_study (study_id),
    KEY ix_sample_subject (subject_id),
    KEY ix_sample_disease (disease_id),
    CONSTRAINT fk_sample_study FOREIGN KEY (study_id) REFERENCES studies(id)
        ON UPDATE CASCADE ON DELETE CASCADE,
    CONSTRAINT fk_sample_disease FOREIGN KEY (disease_id) REFERENCES diseases(id)
        ON UPDATE CASCADE ON DELETE SET NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE IF NOT EXISTS sample_taxon_abundances (
    sample_id BIGINT UNSIGNED NOT NULL,
    taxon_id BIGINT UNSIGNED NOT NULL,
    relative_abundance DOUBLE NOT NULL,
    detection_threshold DOUBLE NULL,
    PRIMARY KEY (sample_id, taxon_id),
    KEY ix_abundance_taxon (taxon_id),
    CONSTRAINT fk_abundance_sample FOREIGN KEY (sample_id) REFERENCES samples(id)
        ON UPDATE CASCADE ON DELETE CASCADE,
    CONSTRAINT fk_abundance_taxon FOREIGN KEY (taxon_id) REFERENCES taxa(id)
        ON UPDATE CASCADE ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE OR REPLACE VIEW v_taxon_disease_summary AS
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


-- Per-sample ML readiness.
--
-- A sample is only usable for supervised per-sample modelling if it has BOTH
-- covariates (sex/age) AND a taxonomic feature vector in sample_taxon_abundances.
-- Those two are tracked separately here on purpose: demographic completeness and
-- feature availability are independent gaps in this database, and collapsing them
-- into a single "ready" flag hides which one is actually blocking a given cohort.
CREATE OR REPLACE VIEW v_ml_ready_samples AS
SELECT
    s.id                                   AS sample_id,
    s.run_accession,
    st.project_accession,
    st.source_database                     AS study_source,
    st.data_type,
    d.name                                 AS phenotype,
    d.mesh_id,
    (d.name = 'Health')                    AS is_control,
    s.sex,
    s.age_years,
    s.bmi,
    s.country,
    COALESCE(s.body_site, 'unspecified')   AS body_site,
    (s.sex IS NOT NULL)                    AS has_sex,
    (s.age_years IS NOT NULL)              AS has_age,
    (s.sex IS NOT NULL AND s.age_years IS NOT NULL) AS has_demographics,
    COALESCE(ab.n_taxa, 0)                 AS n_taxa_profiled,
    (COALESCE(ab.n_taxa, 0) > 0)           AS has_features,
    (s.sex IS NOT NULL
     AND s.age_years IS NOT NULL
     AND COALESCE(ab.n_taxa, 0) > 0)       AS ml_ready
FROM samples s
JOIN studies st ON st.id = s.study_id
LEFT JOIN diseases d ON d.id = s.disease_id
LEFT JOIN (
    SELECT sample_id, COUNT(*) AS n_taxa
    FROM sample_taxon_abundances
    GROUP BY sample_id
) ab ON ab.sample_id = s.id;


-- Study-level rollup: where each cohort stands on covariates vs features.
-- Demographic missingness in this database is study-level rather than random
-- (a study reports sex/age for everyone or for nobody), so the study is the
-- correct unit for reasoning about what is recoverable and what must be modelled
-- as a batch effect.
CREATE OR REPLACE VIEW v_ml_cohort_summary AS
SELECT
    st.project_accession,
    st.source_database,
    st.data_quality,
    COUNT(*)                                AS n_samples,
    SUM(s.sex IS NOT NULL)                  AS n_with_sex,
    SUM(s.age_years IS NOT NULL)            AS n_with_age,
    SUM(s.bmi IS NOT NULL)                  AS n_with_bmi,
    SUM(ab.sample_id IS NOT NULL)           AS n_with_features,
    SUM(d.name = 'Health')                  AS n_control,
    SUM(d.name <> 'Health')                 AS n_case,
    COUNT(DISTINCT s.disease_id)            AS n_phenotypes,
    CASE
        WHEN SUM(s.sex IS NOT NULL) = 0 THEN 'no_demographics'
        WHEN SUM(s.sex IS NOT NULL) = COUNT(*) THEN 'full_demographics'
        ELSE 'partial_demographics'
    END                                     AS demographic_status
FROM samples s
JOIN studies st ON st.id = s.study_id
LEFT JOIN diseases d ON d.id = s.disease_id
LEFT JOIN (
    SELECT DISTINCT sample_id FROM sample_taxon_abundances
) ab ON ab.sample_id = s.id
GROUP BY st.project_accession, st.source_database, st.data_quality;


-- Evidence strength per (taxon, disease).
--
-- 75% of pairs in this table rest on a single study, and 913 pairs have studies
-- that disagree on direction, but every association row looks identical in the
-- base schema. This view exposes replication depth and directional agreement so
-- a 15-study unanimous finding can be told apart from a one-off.
--
-- Agreement is a RATIO, not a boolean: most disagreements here are lopsided
-- majorities with a single outlier (16 studies enriched vs 1 depleted is 94%
-- agreement, not a genuine controversy). Callers pick their own threshold.
--
-- Studies, not association rows, are the unit of agreement, so a single study
-- contributing several comparisons cannot outvote several independent cohorts.
CREATE OR REPLACE VIEW v_taxon_disease_evidence AS
SELECT
    e.taxon_id,
    e.scientific_name,
    e.taxonomic_rank,
    e.phylum,
    e.disease_id,
    e.disease_name,
    e.mesh_id,
    e.n_associations,
    e.n_studies,
    e.n_studies_enriched,
    e.n_studies_depleted,
    e.n_sources,
    e.mean_abs_effect_size,
    CASE
        WHEN e.n_studies_enriched > e.n_studies_depleted THEN 'enriched'
        WHEN e.n_studies_depleted > e.n_studies_enriched THEN 'depleted'
        WHEN (e.n_studies_enriched + e.n_studies_depleted) = 0 THEN NULL
        ELSE 'tied'
    END AS consensus_direction,
    CASE
        WHEN (e.n_studies_enriched + e.n_studies_depleted) > 0
        THEN ROUND(
            GREATEST(e.n_studies_enriched, e.n_studies_depleted)
            / (e.n_studies_enriched + e.n_studies_depleted), 3)
    END AS agreement_ratio,
    (e.n_studies_enriched > 0 AND e.n_studies_depleted > 0) AS has_conflict
FROM (
    SELECT
        a.taxon_id,
        t.scientific_name,
        t.taxonomic_rank,
        t.phylum,
        a.disease_id,
        d.name AS disease_name,
        d.mesh_id,
        COUNT(*) AS n_associations,
        COUNT(DISTINCT c.study_id) AS n_studies,
        COUNT(DISTINCT CASE WHEN a.direction = 'enriched' THEN c.study_id END) AS n_studies_enriched,
        COUNT(DISTINCT CASE WHEN a.direction = 'depleted' THEN c.study_id END) AS n_studies_depleted,
        COUNT(DISTINCT a.source_database) AS n_sources,
        ROUND(AVG(ABS(a.effect_size)), 3) AS mean_abs_effect_size
    FROM taxon_disease_associations a
    JOIN taxa t ON t.id = a.taxon_id
    JOIN diseases d ON d.id = a.disease_id
    JOIN phenotype_comparisons c ON c.id = a.comparison_id
    GROUP BY a.taxon_id, t.scientific_name, t.taxonomic_rank, t.phylum,
             a.disease_id, d.name, d.mesh_id
) e;


-- How disease-specific is each taxon?
--
-- Faecalibacterium is depleted across 74 different diseases: a strong
-- disease-vs-health marker and a near-useless disease-vs-disease one. That
-- distinction is invisible in the base tables, and picking such a taxon as a
-- discriminative feature is a silent modelling error rather than a loud one.
--
-- depleted_fraction near 1 across many diseases is the signature of a general
-- dysbiosis marker (loss of a commensal in illness generally) rather than
-- anything specific to one condition.
CREATE OR REPLACE VIEW v_taxon_specificity AS
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
    CASE
        WHEN s.n_diseases = 1 THEN 'single_disease'
        WHEN s.n_diseases <= 5 THEN 'narrow'
        WHEN s.n_diseases <= 20 THEN 'broad'
        ELSE 'pan_disease'
    END AS specificity_class,
    CASE
        WHEN (s.n_enriched + s.n_depleted) > 0
        THEN ROUND(s.n_depleted / (s.n_enriched + s.n_depleted), 3)
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
        SUM(a.direction = 'enriched') AS n_enriched,
        SUM(a.direction = 'depleted') AS n_depleted
    FROM taxon_disease_associations a
    JOIN taxa t ON t.id = a.taxon_id
    JOIN phenotype_comparisons c ON c.id = a.comparison_id
    GROUP BY a.taxon_id, t.scientific_name, t.taxonomic_rank, t.phylum
) s;


-- Rank-safe access to the abundance matrix.
--
-- A sample can carry more than one taxonomic rank: genus for every sample, and
-- species additionally where the sequencing supports it (shotgun only; GMrepo
-- does not report species for 16S runs at all). Abundances therefore sum to 1
-- WITHIN a rank, not across the sample, so an unfiltered
--   SELECT ... FROM sample_taxon_abundances
-- over a dual-rank sample returns a mixture that sums to the number of ranks
-- present and means nothing.
--
-- These views make the rank-filtered path the easy one, so callers do not have
-- to remember the invariant to get a correct feature matrix.

CREATE OR REPLACE VIEW v_abundance_genus AS
SELECT
    a.sample_id,
    s.run_accession,
    s.subject_id,
    s.study_id,
    s.disease_id,
    a.taxon_id,
    t.scientific_name AS genus,
    t.phylum,
    t.superkingdom,
    a.relative_abundance
FROM sample_taxon_abundances a
JOIN taxa t ON t.id = a.taxon_id
JOIN samples s ON s.id = a.sample_id
WHERE t.taxonomic_rank = 'genus';


CREATE OR REPLACE VIEW v_abundance_species AS
SELECT
    a.sample_id,
    s.run_accession,
    s.subject_id,
    s.study_id,
    s.disease_id,
    a.taxon_id,
    t.scientific_name AS species,
    t.genus,
    t.phylum,
    t.superkingdom,
    a.relative_abundance
FROM sample_taxon_abundances a
JOIN taxa t ON t.id = a.taxon_id
JOIN samples s ON s.id = a.sample_id
WHERE t.taxonomic_rank = 'species';


-- Which ranks each sample actually has, and whether each rank sums to 1.
-- Intended as the first thing to check before building a feature matrix:
-- a sample listed here with sums_to_one = 0 should not be fed to a model.
CREATE OR REPLACE VIEW v_abundance_coverage AS
SELECT
    a.sample_id,
    s.run_accession,
    t.taxonomic_rank,
    COUNT(*) AS n_taxa,
    ROUND(SUM(a.relative_abundance), 6) AS rank_total,
    (ABS(SUM(a.relative_abundance) - 1.0) <= 0.01) AS sums_to_one,
    s.unclassified_fraction
FROM sample_taxon_abundances a
JOIN taxa t ON t.id = a.taxon_id
JOIN samples s ON s.id = a.sample_id
GROUP BY a.sample_id, s.run_accession, t.taxonomic_rank, s.unclassified_fraction;
