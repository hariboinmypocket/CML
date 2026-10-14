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
    -- The rank NCBI assigns this tax ID, copied verbatim from nodes.dmp. It sits
    -- BESIDE taxonomic_rank rather than replacing it, because the two answer
    -- different questions: taxonomic_rank is the abundance LAYER the row was
    -- loaded into, and ncbi_rank is what the organism actually is.
    --
    -- 408 rows disagree. 125 sit in the genus layer while NCBI calls them
    -- coarser -- 90 families (Muribaculaceae 2005473, Tannerellaceae 2005525),
    -- 19 orders, 7 classes, 5 phyla -- and 283 sit in the species layer while
    -- NCBI calls them finer, almost all strains.
    --
    -- THE GENUS-LAYER ONES MUST NOT BE "FIXED" BY MOVING THEM. GMrepo's own
    -- export labels them taxonomic_rank = genus, and they are the unclassified
    -- RESIDUAL within that family: their mean share is 0.2-0.6% while the real
    -- genera beside them carry far more, so they are the reads a 16S classifier
    -- placed to family but no further. 3,443 samples hold such a row alongside
    -- genuine genera of the same family; the parts are disjoint, so the
    -- per-(sample, rank) sum to 1 is correct. Taking them out of the genus layer
    -- would break that invariant and discard real mass.
    --
    -- What was missing is only that nothing SAID they were coarser, so a model
    -- treated Muribaculaceae as a peer of Bacteroides. v_abundance_genus exposes
    -- coarser_than_genus for exactly that.
    ncbi_rank VARCHAR(32) NULL,
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
    -- Whether this genus responds to disease NON-specifically, from MicrobiomeHD
    -- (Duvallet et al. 2017, file-S3): a genus qualifies when it is significant
    -- at q < 0.05 in the same direction in at least two DIFFERENT diseases.
    -- 'health' and 'disease' say which direction, 'mixed' means it is enriched
    -- in cases for two diseases and in controls for two others.
    --
    -- This is one meta-analysis's verdict over 28 studies, not an intrinsic
    -- property, and NULL means "not assessed by that analysis" rather than
    -- "specific to one disease". Its value is directional: v_taxon_specificity
    -- can tell that a genus is pan-disease but not whether it marks health or
    -- disease, which is exactly what this adds.
    nonspecific_response ENUM('health', 'disease', 'mixed') NULL,
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
    -- 'doi:<doi>' for each of the 13 studies that carry metabolite levels: the
    -- DOI of the paper describing the protocol, NOT an instrument name.
    --
    -- The distinction is the whole point. The Muller collection redistributes
    -- the matrices, not the protocols, so no platform string exists to store --
    -- and the column spent its life holding the literal 'see mtb.map.tsv', a
    -- file that maps compound names to KEGG and HMDB ids and contains no method
    -- at all. Every study gave the same answer, so the standing rule below was
    -- uncheckable. A DOI makes it mechanical instead: levels are poolable only
    -- where it matches, and no two cohorts here share one.
    metabolomics_method VARCHAR(120) NULL,
    -- 'centred' or 'non_negative', derived from the study's own level range by
    -- derive_metabolite_scale(). Two of the 13 cohorts arrive already log- or
    -- z-transformed and hold 70,067 negative levels; the other 11 are
    -- non-negative, with maxima spanning nine orders of magnitude (461 to
    -- 3.07e10).
    --
    -- THE ZERO IS WHY THIS COLUMN EXISTS. On a non-negative scale a zero means
    -- the compound was not detected -- 78,400 of YACHIDA_CRC_2019's 130,242
    -- levels, 60% of them. On a centred scale a zero is mid-distribution:
    -- SINHA_CRC_2016 holds 23,553 negatives against 22,540 positives with a mean
    -- of -0.04, so its 281 zeros mean "at this cohort's average" -- the opposite
    -- of absence. The same stored number meant opposite things with nothing
    -- recording which, and "normalize within study" does not cover it: a uniform
    -- log1p over this table returns NaN on 70,067 rows and takes the log of a
    -- log on 136,134 more.
    --
    -- 'non_negative' claims no unit, only non-negativity. The source records
    -- units nowhere; metabolomics_method above carries the protocol DOI.
    metabolite_scale VARCHAR(20) NULL,
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
    -- Significance, where the source reports it. These are independent of
    -- effect_size, not derived from it: GMrepo's LEfSe export carries an LDA
    -- score and no p/q value, while MicrobiomeHD's standardized re-analysis
    -- carries an FDR-corrected q value and no effect size. A row may therefore
    -- have one, the other, or both, and a NULL means "the source did not report
    -- it" rather than "not significant".
    p_value DOUBLE NULL,
    q_value DOUBLE NULL,
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
    -- GMrepo's own QC verdict for the run, '1' passed and '0' failed, copied
    -- from its QCStatus field. The vocabulary belongs to that source: a sample
    -- GMrepo never assessed has NULL here, which is why 3,342 rows are empty.
    -- Do not use it to record provenance or curation tier -- studies.data_quality
    -- is for that, and a loader that wrote 'curated' here put a word from a
    -- different axis into a column others filter as a boolean.
    qc_status VARCHAR(100) NULL,
    -- Controlled vocabulary, two values: 'feces' (voided stool, 36,573) and
    -- 'gastric mucosa' (endoscopic biopsy, 120). It is a PER-STUDY CONSTANT
    -- everywhere -- no source in this database reports it per sample -- and for
    -- every source but MGnify it is the source's documented scope rather than a
    -- reading, which data/body_site_audit.csv records in its `basis` column.
    -- It used to be NULL for 36,094 rows in MySQL while the SQLite mirror said
    -- 'gut' for the same samples, so the same query answered differently per
    -- store. 'gut' was also the wrong word: the stomach is part of the gut, so
    -- it could not separate the biopsies from the stool, which is the only
    -- distinction this column is for. scripts/normalize_body_site.py folds both.
    -- The platform that produced THIS sample's profile, which studies.data_type
    -- cannot express. GMrepo records it per run -- 'AMPLICON' or 'Metagenomics'
    -- in data/gmrepo_project_samples.csv -- and 10 projects contain both, so the
    -- loader collapsed them to data_type = 'mixed' and 3,786 samples lost a
    -- platform the source had for every one of them (2,683 amplicon, 1,103
    -- shotgun).
    --
    -- That loss is not cosmetic. Platform sets the taxonomic floor: amplicon
    -- reads one marker gene and cannot separate species within a genus, so
    -- inside a 'mixed' study some samples resolve species and some stop at genus
    -- with nothing saying which. It is also the strongest batch covariate here
    -- -- scripts/loso_audit.py measures platform mixing at +0.175 optimism --
    -- and both that audit and export_ml_matrix.py read the STUDY's value, so for
    -- these samples the covariate they carried was the word 'mixed'.
    sequencing_platform VARCHAR(40) NULL,
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
    -- detection_threshold was dropped: no source this pipeline ingests has
    -- ever supplied the column, so it was NULL in all 2.58M rows.
    PRIMARY KEY (sample_id, taxon_id),
    KEY ix_abundance_taxon (taxon_id),
    CONSTRAINT fk_abundance_sample FOREIGN KEY (sample_id) REFERENCES samples(id)
        ON UPDATE CASCADE ON DELETE CASCADE,
    CONSTRAINT fk_abundance_taxon FOREIGN KEY (taxon_id) REFERENCES taxa(id)
        ON UPDATE CASCADE ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

-- Sparse per-sample metadata that does not deserve a column.
--
-- Cohorts carry traits that matter for confound control but that almost no
-- other cohort reports: smoking, alcohol, blood pressure, comorbidities,
-- surgery type. Across the Muller 2022 collection alone that is 43,700
-- populated values spread over 366 distinct column names, and ERAWIJANTARI
-- contributes 41 clinical columns by itself. Widening `samples` to hold them
-- would mean 366 columns at 95-99% NULL each -- smoking alone would be
-- populated for 1,011 of 36,693 samples, under 3%.
--
-- So absence is the absence of a row, exactly as in sample_taxon_abundances.
-- There are no empty cells here by construction, and the whole table is ~44k
-- rows against that one's 3.07M.
--
-- VALUES ARE NOT COMPARABLE ACROSS STUDIES, and often not even commensurable.
-- Smoking arrives as `Brinkman Index` (a pack-year product), `SmokingStatus`
-- (categorical), `Tobacco_Use` and `smoking status` -- four names, four scales,
-- one concept. Nothing is harmonized on the way in, because harmonizing would
-- mean discarding the original scale, and the use case is per-cohort anyway:
-- confound control belongs inside a study, which is also where
-- scripts/adjudicate_conflicts.py does its comparisons. Filter by study before
-- trusting an attribute name to mean one thing.
CREATE TABLE IF NOT EXISTS sample_attributes (
    sample_id BIGINT UNSIGNED NOT NULL,
    -- The source's own column name, kept verbatim rather than cleaned. Some are
    -- R make.names artifacts (`Weight..kg.`, `Gout...22` and `Gout...46`, which
    -- are two different gout columns the source named identically); renaming
    -- them would break the link back to the file a value came from.
    attribute VARCHAR(80) NOT NULL,
    -- 512 rather than 255: the longest value in the Muller collection is a
    -- 289-character clinical course note. Loaders report truncation instead of
    -- silently cutting.
    value VARCHAR(512) NOT NULL,
    -- A number you may do arithmetic on -- which is NARROWER than "value parsed
    -- as a number", and the earlier wording here claimed the looser thing while
    -- the data broke even that. Two classes were cleared by
    -- scripts/fix_attribute_numerics.py:
    --
    --   434 rows under `unique_id` held '10007539_0' as 100075390, the
    --   underscore dropped and the digits run together. Not a parse of the
    --   value, and it collides: '1000753_90' gives the same number.
    --
    --   687 rows under two keys whose names spell out their own code mapping
    --   ('..._1_no_2_yes') stored 1 and 2 as labels. Averaging labels is
    --   meaningless, and the encoding INVERTS against its neighbour -- 1 is no
    --   there, 1 is yes under `Antibiotics` -- so pooling them flips 343 rows.
    --
    -- 13 keys still MIX numbers with words and are left that way, because the
    -- number is correct wherever there is one; what is missing is the
    -- denominator. AVG(value_numeric) over `Stage` is the mean of 27 zeros out
    -- of 150 rows whose real values run I to IVB, and under
    -- `Alcohol.drink.week.` the non-numeric values are RANGES, so the numeric
    -- subset is precisely the respondents who gave a point answer rather than a
    -- random sample. Read any mean against v_sample_attribute_keys.
    value_numeric DOUBLE NULL,
    PRIMARY KEY (sample_id, attribute),
    KEY ix_sample_attribute (attribute),
    KEY ix_sample_attribute_numeric (attribute, value_numeric),
    CONSTRAINT fk_sample_attribute_sample FOREIGN KEY (sample_id) REFERENCES samples(id)
        ON UPDATE CASCADE ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

-- Per-key coverage for sample_attributes, so a mean always has its denominator.
--
-- 362 keys, and 13 of them mix numbers with words. Without this view
-- `SELECT AVG(value_numeric) ... WHERE attribute = 'Stage'` returns 0 and looks
-- like an answer: it is the mean of the 27 rows that happen to hold a numeric
-- stage out of 150 whose values run I to IVB. numeric_fraction is what tells you
-- not to trust it. Same shape as absence in the abundance matrix -- the number
-- was never wrong, the denominator was missing.
--
-- `mixes_numbers_and_words` is the flag to filter on. `looks_like_an_interval`
-- narrows it further to the keys whose non-numeric values are RANGES ('0-1',
-- '2-4', '3+'), where the numeric subset is biased by construction.
--
-- The activity keys show how bad that gets. Under
-- `PhysicalExercise_LastWeek.h.` the ONLY point answer is 0, given by 4 of 49
-- respondents; the other 45 answered '1-', '1-3' or '3+'. So
-- AVG(value_numeric) is exactly 0.0 -- not merely incomplete but anti-correlated
-- with the truth, because the rows that carry a number are precisely the people
-- who did nothing. Cycling_LastWeek.h. and Housework_Childcare_LastWeek.h.
-- are the same shape.
CREATE OR REPLACE VIEW v_sample_attribute_keys AS
SELECT
    a.attribute,
    COUNT(*)                                   AS n_rows,
    COUNT(DISTINCT a.sample_id)                AS n_samples,
    COUNT(DISTINCT a.value)                    AS n_distinct_values,
    SUM(a.value_numeric IS NOT NULL)           AS n_numeric,
    ROUND(SUM(a.value_numeric IS NOT NULL) / COUNT(*), 3) AS numeric_fraction,
    (SUM(a.value_numeric IS NOT NULL) > 0
     AND SUM(a.value_numeric IS NOT NULL) < COUNT(*)) AS mixes_numbers_and_words,
    SUM(a.value REGEXP '^[0-9]+(\\.[0-9]+)?[-+]') > 0 AS looks_like_an_interval,
    MIN(a.value_numeric)                       AS numeric_min,
    MAX(a.value_numeric)                       AS numeric_max,
    AVG(a.value_numeric)                       AS numeric_mean
FROM sample_attributes a
GROUP BY a.attribute;


-- The metabolome axis.
--
-- Added for the Muller et al. 2022 collection (doi:10.1038/s41522-022-00345-5),
-- which pairs faecal microbiome profiles with metabolite levels for 14 cohorts.
-- Metabolites are a genuinely different feature class from taxa, not more rows
-- of the same, and they obey none of the invariants sample_taxon_abundances
-- does.
--
-- LEVELS ARE NOT COMPARABLE ACROSS STUDIES. Instruments, targeted versus
-- untargeted designs, extraction protocols and units all differ between
-- cohorts, and nothing sums to 1. A model may rank or correlate levels WITHIN a
-- study, or compare a study's case arm against its own control arm, but pooling
-- raw levels across studies is meaningless in a way that pooling relative
-- abundances is not. studies.metabolomics_method holds the protocol's DOI so
-- that constraint is checkable rather than merely stated: two studies' levels
-- may be pooled only where that value matches, and in this database no two
-- cohorts share one.
CREATE TABLE IF NOT EXISTS metabolites (
    id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    -- Normalized to the 7-digit form: upstream carries both HMDB00002 and
    -- HMDB0000002 for the same compound. Either identifier may be absent, and
    -- MySQL permits repeated NULLs under a UNIQUE key, so a compound known only
    -- to KEGG and one known only to HMDB both store cleanly.
    hmdb_id VARCHAR(16) NULL,
    kegg_id VARCHAR(10) NULL,
    name VARCHAR(255) NOT NULL,
    PRIMARY KEY (id),
    UNIQUE KEY uq_metabolite_hmdb (hmdb_id),
    UNIQUE KEY uq_metabolite_kegg (kegg_id),
    KEY ix_metabolite_name (name)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE IF NOT EXISTS sample_metabolite_levels (
    sample_id BIGINT UNSIGNED NOT NULL,
    metabolite_id BIGINT UNSIGNED NOT NULL,
    level DOUBLE NOT NULL,
    -- Upstream flags annotations the original authors considered uncertain
    -- (High.Confidence.Annotation=FALSE: 643 of 98,933 features). Carried
    -- rather than filtered, because which confidence floor is acceptable is the
    -- caller's decision, not the loader's.
    high_confidence_annotation BOOLEAN NOT NULL DEFAULT TRUE,
    PRIMARY KEY (sample_id, metabolite_id),
    KEY ix_metabolite_level_metabolite (metabolite_id),
    CONSTRAINT fk_metabolite_level_sample FOREIGN KEY (sample_id) REFERENCES samples(id)
        ON UPDATE CASCADE ON DELETE CASCADE,
    CONSTRAINT fk_metabolite_level_metabolite FOREIGN KEY (metabolite_id) REFERENCES metabolites(id)
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
--
-- has_features asks only whether n_taxa_profiled is above zero, which is a lower
-- bar than it sounds. 613 (sample, rank) profiles hold two taxa or fewer: 269
-- with a single genus, 316 with two genera, 28 with two species. All of them are
-- faithful to GMrepo, whose own numbers for them sum to exactly 1.0, so these
-- are reported upstream as complete profiles and there is nothing here that was
-- lost at load. A one-part composition has a CLR of identically zero, so such a
-- row is not a feature vector in any usable sense, and gutdb/compositional.py
-- drops anything under two parts because the transform is undefined there. The
-- median sample holds 62 genera and the 5th percentile 9.
--
-- FILTER ON n_genera_profiled OR n_species_profiled, NEVER ON n_taxa_profiled.
-- n_taxa_profiled counts rows across every rank at once, so for the 12,566
-- samples profiled at both genus and species it is a sum of two richnesses and
-- not the richness of anything. 44 samples are thin at one rank while the sum
-- looks healthy, and a filter written against the sum keeps them.
--
-- Depth is also mostly a STUDY-level property rather than a sample-level one:
-- six studies holding 516 samples have no sample at all above 15 genera, so
-- inside them a thin profile is the cohort and not an outlier. That matters
-- because depth can predict the label on its own -- Neuromyelitis Optica's
-- entire case arm is one such study at a median of 3 genera against controls at
-- 45, and richness used alone as the only predictor scores AUC 0.982 there.
-- v_ml_cohort_summary carries the per-study depth columns and
-- scripts/richness_audit.py measures the confound per cohort.
-- Those two are tracked separately here on purpose: demographic completeness and
-- feature availability are independent gaps in this database, and collapsing them
-- into a single "ready" flag hides which one is actually blocking a given cohort.
CREATE OR REPLACE VIEW v_ml_ready_samples AS
SELECT
    s.id                                   AS sample_id,
    s.run_accession,
    st.project_accession,
    st.source_database                     AS study_source,
    -- The STUDY's summary. For the 10 projects holding both platforms it reads
    -- 'mixed', which is not a batch and not a resolution; use
    -- sequencing_platform for anything per-sample.
    st.data_type,
    s.sequencing_platform,
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
    COALESCE(ab.n_genera, 0)               AS n_genera_profiled,
    COALESCE(ab.n_species, 0)              AS n_species_profiled,
    -- Richness at the sample's deepest profiled rank, which is the number a
    -- compositional transform is actually applied over. GREATEST rather than a
    -- sum: the ranks are separate compositions, each summing to 1 on its own.
    GREATEST(COALESCE(ab.n_genera, 0),
             COALESCE(ab.n_species, 0))    AS profile_richness,
    (COALESCE(ab.n_taxa, 0) > 0)           AS has_features,
    -- Two parts is the mathematical floor, not a quality bar: below it the CLR
    -- is identically zero. Anything stricter is a modelling choice and is left
    -- to whoever is modelling.
    (GREATEST(COALESCE(ab.n_genera, 0),
              COALESCE(ab.n_species, 0)) >= 2) AS has_transformable_features,
    -- The biological sample, not the sequencing run. A study may sequence one
    -- sample several times, so `samples` can hold two to four rows that are one
    -- piece of biological material: PRJEB28543 profiles each of its samples four
    -- times, and 681 sample groups here cover 1,575 rows. Their profiles differ
    -- genuinely (independent runs, independently profiled), so none is a
    -- duplicate to delete -- but treating them as independent observations is
    -- pseudo-replication, inflating the effective sample size and putting the
    -- same material on both sides of a train/test split.
    s.gmrepo_sample_id                     AS biological_sample_id,
    COALESCE(rep.n_runs, 1)                AS runs_for_this_sample,
    (COALESCE(rep.n_runs, 1) > 1
     AND s.id <> rep.first_id)             AS is_replicate_run,
    (s.sex IS NOT NULL
     AND s.age_years IS NOT NULL
     AND COALESCE(ab.n_taxa, 0) > 0)       AS ml_ready,
    -- ml_ready says the row has covariates and features; this additionally
    -- keeps one run per biological sample.
    --
    -- IT DOES NOT KEEP ONE SAMPLE PER PERSON, and an earlier version of this
    -- comment claimed it was "what a split should use", which overstated it. A
    -- longitudinal study sampling one subject at seven timepoints produces seven
    -- distinct biological samples, all of which survive this flag: 389 subjects
    -- contribute 1,340 of its rows, so 971 are repeat visits by someone already
    -- in the cohort. Irritable Bowel Syndrome is 305 rows from 51 people and
    -- Diabetes Mellitus, Type 1 is 66 rows from 11 -- a sixfold overstatement of
    -- the effective sample size, and a random split puts the same person on both
    -- sides of it. Use ml_ready_one_per_subject for a split.
    (s.sex IS NOT NULL
     AND s.age_years IS NOT NULL
     AND COALESCE(ab.n_taxa, 0) > 0
     AND (COALESCE(rep.n_runs, 1) = 1 OR s.id = rep.first_id)) AS ml_ready_deduplicated,
    s.subject_id,
    COALESCE(subj.n_visits, 1)             AS samples_for_this_subject,
    -- One row per person where the source names the person, falling back to the
    -- sample-level dedup where it does not.
    --
    -- THE FALLBACK IS MOST OF THE COHORT, so read this flag knowing what it can
    -- and cannot see. subject_id is NULL for all 33,630 GMrepo samples and all
    -- 120 from MGnify, which is 10,961 of the 13,328 deduplicated rows. A
    -- longitudinal GMrepo study is therefore invisible to this flag, and it only
    -- genuinely constrains the 2,367 MullerMM and PubMed rows. It is an
    -- improvement on ml_ready_deduplicated, not a solution to
    -- pseudo-replication.
    (s.sex IS NOT NULL
     AND s.age_years IS NOT NULL
     AND COALESCE(ab.n_taxa, 0) > 0
     AND (COALESCE(rep.n_runs, 1) = 1 OR s.id = rep.first_id)
     AND (s.subject_id IS NULL OR s.id = subj.first_visit_id)) AS ml_ready_one_per_subject
FROM samples s
JOIN studies st ON st.id = s.study_id
LEFT JOIN diseases d ON d.id = s.disease_id
LEFT JOIN (
    SELECT a.sample_id,
           COUNT(*)                                        AS n_taxa,
           SUM(t.taxonomic_rank = 'genus')                 AS n_genera,
           SUM(t.taxonomic_rank = 'species')               AS n_species
    FROM sample_taxon_abundances a
    JOIN taxa t ON t.id = a.taxon_id
    GROUP BY a.sample_id
) ab ON ab.sample_id = s.id
LEFT JOIN (
    SELECT study_id, gmrepo_sample_id, COUNT(*) AS n_runs, MIN(id) AS first_id
    FROM samples
    WHERE gmrepo_sample_id IS NOT NULL
    GROUP BY study_id, gmrepo_sample_id
) rep ON rep.study_id = s.study_id
     AND rep.gmrepo_sample_id = s.gmrepo_sample_id
LEFT JOIN (
    SELECT study_id, subject_id, COUNT(*) AS n_visits,
           MIN(id) AS first_visit_id
    FROM samples
    WHERE subject_id IS NOT NULL
    GROUP BY study_id, subject_id
) subj ON subj.study_id = s.study_id
      AND subj.subject_id = s.subject_id;


-- Study-level rollup: where each cohort stands on covariates vs features.
-- Demographic missingness in this database is study-level rather than random
-- (a study reports sex/age for everyone or for nobody), so the study is the
-- correct unit for reasoning about what is recoverable and what must be modelled
-- as a batch effect.
--
-- Taxonomic DEPTH is study-level for the same reason, and the columns for it are
-- here rather than only per sample because the per-sample view cannot express
-- the problem. Six studies -- PRJNA1031953, PRJNA832282, PRJDB6133, PRJNA422961,
-- PRJNA923180, PRJNA889687, 516 samples between them -- contain no sample above
-- 15 genera, so a per-sample richness filter either keeps all of them or drops
-- every one, and there is no outlier inside them to catch. feature_depth_status
-- names that situation: `genus_depth_degenerate` means nothing in the study can
-- support a genus-level model, whatever the per-sample numbers say.
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
    MAX(COALESCE(ab.n_genera, 0))           AS max_genera,
    ROUND(AVG(COALESCE(ab.n_genera, 0)), 1) AS mean_genera,
    SUM(COALESCE(ab.n_genera, 0) BETWEEN 1 AND 2) AS n_thin_genus_profiles,
    CASE
        WHEN SUM(ab.sample_id IS NOT NULL) = 0 THEN 'no_features'
        WHEN MAX(COALESCE(ab.n_genera, 0)) = 0 THEN 'no_genus_layer'
        WHEN MAX(COALESCE(ab.n_genera, 0)) < 15 THEN 'genus_depth_degenerate'
        ELSE 'ok'
    END                                     AS feature_depth_status,
    CASE
        WHEN SUM(s.sex IS NOT NULL) = 0 THEN 'no_demographics'
        WHEN SUM(s.sex IS NOT NULL) = COUNT(*) THEN 'full_demographics'
        ELSE 'partial_demographics'
    END                                     AS demographic_status
FROM samples s
JOIN studies st ON st.id = s.study_id
LEFT JOIN diseases d ON d.id = s.disease_id
LEFT JOIN (
    SELECT a.sample_id,
           SUM(t.taxonomic_rank = 'genus') AS n_genera
    FROM sample_taxon_abundances a
    JOIN taxa t ON t.id = a.taxon_id
    GROUP BY a.sample_id
) ab ON ab.sample_id = s.id
GROUP BY st.project_accession, st.source_database, st.data_quality;


-- Evidence strength per (taxon, disease).
--
-- 80% of pairs in this table rest on a single study, and the studies behind the
-- rest do not always agree, but every association row looks identical in the
-- base schema. This view exposes replication depth and directional agreement so
-- a 15-study unanimous finding can be told apart from a one-off.
--
-- TWO RULES DECIDE WHETHER A DIRECTION HERE MEANS ANYTHING, and an earlier
-- version of this view broke both. It reported Faecalibacterium prausnitzii as
-- ENRICHED in ulcerative colitis, reversing one of the most replicated findings
-- in the field, and it was not alone: enforcing the first rule below reversed 33
-- directions, 27 of them in ulcerative colitis and 6 in Crohn disease, almost
-- all of them butyrate-producing commensals whose depletion in IBD is textbook
-- (Faecalibacterium, Lachnospira, Agathobacter rectalis, Dorea, Ruminococcus).
-- The abundance matrix in this same database independently backs the corrected
-- direction over the old one 24 times to 7.
--
--   Only case/control contrasts may vote. 2,340 of the 16,889 associations
--   compare one disease against another rather than against Health. "Enriched
--   in UC relative to Crohn disease" and "enriched in UC relative to health"
--   are different claims, and a model trained on case/control data needs the
--   second. Pairs whose associations are ALL disease-vs-disease therefore get
--   consensus_direction = NULL rather than a direction built from the wrong
--   question; contrast_scope says which situation a pair is in.
--
--   One study casts one vote. n_studies_enriched and n_studies_depleted used to
--   be two independent COUNT(DISTINCT study_id) expressions, so a study holding
--   both an enriched and a depleted association was counted in BOTH sides at
--   once. A study that disagrees with itself now lands in n_studies_split and
--   votes for neither.
--
--   Measured honestly, this second rule currently changes NOTHING: no study in
--   the database reports both directions among its own vs-Health associations,
--   so n_studies_split is 0 for all 11,472 pairs and every one of the 226
--   behaviour changes below traces to the contrast restriction instead. It is
--   kept as a guard rather than a fix -- two comparison methods inside one
--   study could disagree tomorrow, and the old expression would have silently
--   counted that study twice.
--
-- What the contrast restriction actually changed, over the 2,285 replicated
-- pairs: 118 that read 'tied' became decided, 41 that read as decided turned
-- out to be genuine ties, 33 reversed outright, and 34 lost their direction
-- entirely because every association behind them was disease-vs-disease. Every
-- one of those 226 pairs is 'mixed' or 'disease_vs_disease' scope; no pure
-- vs_health pair moved.
--
-- Replication depth (n_associations, n_studies, n_sources) still counts every
-- contrast, because a disease-vs-disease result is a real result; it just
-- cannot answer a case/control question.
--
-- Agreement is a RATIO, not a boolean, so callers pick their own threshold. Of
-- the 2,285 replicated pairs, 387 were exact ties under the old pooled vote and
-- 265 remain exact ties under this one, where consensus_direction is 'tied' and
-- genuinely undecided. taxon_disease_adjudication carries those further by
-- bringing this database's own abundance measurements to bear; this view stays
-- purely a summary of what the curated associations say.
CREATE OR REPLACE VIEW v_taxon_disease_evidence AS
WITH assoc AS (
    -- Every association, tagged with whether its comparison is case/control.
    -- COALESCE because a comparison may carry no phenotype ids at all, and
    -- NULL OR NULL is NULL rather than false.
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
    -- One row per (taxon, disease, study), carrying that study's single vote.
    -- A study whose own case/control associations disagree with itself is
    -- 'split' and votes for neither side, instead of being counted on both.
    SELECT taxon_id,
           disease_id,
           study_id,
           MAX(vs_health) AS has_case_control,
           CASE
               WHEN MAX(vs_health) = 0 THEN 'other_contrast'
               WHEN COUNT(DISTINCT CASE WHEN vs_health THEN direction END) > 1 THEN 'split'
               ELSE MAX(CASE WHEN vs_health THEN direction END)
           END AS vote
    FROM assoc
    GROUP BY taxon_id, disease_id, study_id
),
votes AS (
    SELECT taxon_id,
           disease_id,
           SUM(vote = 'enriched') AS n_studies_enriched,
           SUM(vote = 'depleted') AS n_studies_depleted,
           SUM(vote = 'split') AS n_studies_split,
           SUM(vote = 'other_contrast') AS n_studies_other_contrast,
           SUM(has_case_control = 1) AS n_studies_case_control
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
       -- Replication depth counts EVERY contrast: a disease-vs-disease result
       -- is still a result, it just cannot vote on a case/control direction.
       o.n_associations,
       o.n_studies,
       o.n_sources,
       o.mean_abs_effect_size,
       o.contrast_scope,
       -- Direction counts case/control studies ONLY, one vote each.
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
               GREATEST(v.n_studies_enriched, v.n_studies_depleted)
               / (v.n_studies_enriched + v.n_studies_depleted), 3)
       END AS agreement_ratio,
       (v.n_studies_enriched > 0 AND v.n_studies_depleted > 0) AS has_conflict
FROM totals o
JOIN votes v ON v.taxon_id = o.taxon_id AND v.disease_id = o.disease_id
JOIN taxa t ON t.id = o.taxon_id
JOIN diseases d ON d.id = o.disease_id;


-- How disease-specific is each taxon?
--
-- Faecalibacterium is depleted across dozens of different diseases: a strong
-- disease-vs-health marker and a near-useless disease-vs-disease one. That
-- distinction is invisible in the base tables, and picking such a taxon as a
-- discriminative feature is a silent modelling error rather than a loud one.
--
-- depleted_fraction near 1 across many diseases is the signature of a general
-- dysbiosis marker (loss of a commensal in illness generally) rather than
-- anything specific to one condition.
--
-- depleted_fraction counts CASE/CONTROL associations only, for the same reason
-- v_taxon_disease_evidence does: "depleted in UC relative to Crohn disease"
-- says nothing about depletion in illness. Pooling both kinds moved this
-- fraction by more than 10 percentage points for 85 of the 821 taxa with at
-- least three associations, and by more than 25 points for 21 of them, with a
-- worst case of 0.40 -- enough to turn a dysbiosis marker into an apparently
-- specific one or the reverse. It is NULL where a taxon has no case/control
-- association at all.
--
-- The breadth columns (n_diseases, n_associations, n_studies) still count every
-- contrast, and these counts are association-level rather than one-vote-per-
-- study on purpose: specificity is a question about how widely a taxon has been
-- reported, not about how cohorts voted on one disease.
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
    s.n_other_contrast,
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
        SUM(a.direction = 'enriched'
            AND COALESCE(pa.name = 'Health' OR pb.name = 'Health', 0)) AS n_enriched,
        SUM(a.direction = 'depleted'
            AND COALESCE(pa.name = 'Health' OR pb.name = 'Health', 0)) AS n_depleted,
        SUM(COALESCE(pa.name = 'Health' OR pb.name = 'Health', 0) = 0) AS n_other_contrast
    FROM taxon_disease_associations a
    JOIN taxa t ON t.id = a.taxon_id
    JOIN phenotype_comparisons c ON c.id = a.comparison_id
    LEFT JOIN diseases pa ON pa.id = c.phenotype_a_id
    LEFT JOIN diseases pb ON pb.id = c.phenotype_b_id
    GROUP BY a.taxon_id, t.scientific_name, t.taxonomic_rank, t.phylum
) s;


-- Where studies disagree on a taxon's direction, and what the abundance data
-- says about it. Derived, not authoritative: taxon_disease_associations is left
-- exactly as loaded, and this table records a verdict beside the evidence for it
-- so it can be recomputed or ignored. Populated by scripts/adjudicate_conflicts.py.
CREATE TABLE IF NOT EXISTS taxon_disease_adjudication (
    taxon_id BIGINT UNSIGNED NOT NULL,
    disease_id BIGINT UNSIGNED NOT NULL,
    n_studies INT NOT NULL,
    -- Which contrasts the curated associations for this pair actually use.
    -- 2,340 of the 16,889 associations compare one disease against another
    -- rather than against Health, and they answer a different question than a
    -- case/control model asks, so only the vs-Health ones are allowed to vote.
    -- 'mixed' means both kinds are present; 'disease_vs_disease' means no
    -- case/control association exists and the verdict rests on abundance alone.
    assoc_contrast_scope VARCHAR(20) NULL,
    -- the vote taken from the curated associations, vs-Health contrasts only,
    -- one vote per study
    assoc_enriched INT NOT NULL,
    assoc_depleted INT NOT NULL,
    -- studies whose own vs-Health associations disagree with themselves; these
    -- vote for neither side instead of voting twice
    assoc_split_studies INT NOT NULL DEFAULT 0,
    assoc_consensus VARCHAR(20) NULL,
    assoc_agreement DECIMAL(4,3) NULL,
    -- the vote taken from this database's own abundance measurements, one
    -- direction per study that holds both arms, absence counted as zero
    abundance_enriched INT NOT NULL,
    abundance_depleted INT NOT NULL,
    abundance_studies INT NOT NULL,
    abundance_agreement DECIMAL(4,3) NULL,
    -- NULL verdict means the evidence does not decide, which is a result rather
    -- than a gap: some pairs are genuinely heterogeneous across cohorts.
    verdict ENUM('enriched', 'depleted') NULL,
    basis VARCHAR(60) NOT NULL,
    PRIMARY KEY (taxon_id, disease_id),
    KEY ix_adjudication_verdict (verdict),
    KEY ix_adjudication_basis (basis),
    KEY ix_adjudication_scope (assoc_contrast_scope),
    CONSTRAINT fk_adjudication_taxon FOREIGN KEY (taxon_id) REFERENCES taxa(id)
        ON UPDATE CASCADE ON DELETE CASCADE,
    CONSTRAINT fk_adjudication_disease FOREIGN KEY (disease_id) REFERENCES diseases(id)
        ON UPDATE CASCADE ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;


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
    -- 125 of these features are not genera. See taxa.ncbi_rank: they are the
    -- family- or order-level residual a 16S classifier could not place further,
    -- faithfully reported by GMrepo at genus rank, and disjoint from the real
    -- genera beside them. Filter on this when a model needs features at one
    -- resolution; do not drop the rows, which would break the sum to 1.
    t.ncbi_rank,
    (t.ncbi_rank IS NOT NULL
     AND t.ncbi_rank NOT IN ('genus', 'subgenus', 'species', 'subspecies',
                             'strain', 'no rank', 'clade')) AS coarser_than_genus,
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
