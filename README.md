# Gut microbiome database population pipeline

This project turns the earlier notebook/Workbench workflow into a repeatable MySQL 8 pipeline. It loads the existing gut microbiome CSV, enriches taxa from MiMeDB, and stores GMrepo diseases, studies, comparisons, marker associations, clinical samples, and sample-level relative abundances in normalized tables.

## Why this schema

A bacterium can be related to many diseases, and a disease can be related to many bacteria. The relationship therefore belongs in `taxon_disease_associations`, which contains foreign keys to both `taxa` and `diseases`. `phenotype_comparisons` preserves the study, comparison groups, method, and the meaning of positive/negative LDA scores. This avoids losing the evidence context behind an association.

The core relationship is:

```text
taxa 1---* taxon_disease_associations *---1 diseases
                     |
                     *---1 phenotype_comparisons *---1 studies
```

Raw-source provenance is retained in `taxon_source_records`, while each load is audited in `ingestion_runs`. Failed loads roll back their data changes and retain the error in the audit row. All loaders are idempotent and can be rerun.

## 1. Install and configure

Use Python 3.10 or newer and MySQL 8:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
cp .env.example .env
```

Edit `.env` with the MySQL credentials you use in MySQL Workbench. The MySQL user needs permission to create the configured database and tables.

## 2. Populate the database

Load the two local data sources and the IBS marker set recovered from the supplied workflow:

```bash
python populate_database.py all \
  --associations data/gmrepo_PRJNA705217_ibs_markers.csv
```

This command:

1. Creates `gut_microbiome` and all tables/views.
2. Loads `gut_microbiome_new_full.csv`.
3. Collapses MiMeDB strain rows to species and derives `energy_mode` and `primary_food_source`.
4. Fills missing taxon attributes without overwriting curated values.
5. Loads the PRJNA705217 Health-vs-IBS LDA associations.

To intentionally replace existing taxon attributes/statistics, add `--overwrite`.

## 3. Add more GMrepo projects

GMrepo's documented API exposes phenotype, run, and abundance endpoints, while project marker tables can be downloaded as TSV from the project/comparison page. Normalize a marker export to the headers in [data/gmrepo_associations_template.csv](data/gmrepo_associations_template.csv), then run:

```bash
python populate_database.py load-associations path/to/markers.tsv
```

The required association fields are:

- `project_id`
- `scientific_name` (or separate `genus` and `species`)
- `disease_name`
- `phenotype_a` and `phenotype_b`
- `lda_score` or an explicit `direction`

Always populate `negative_enriched_in` and `positive_enriched_in` when importing signed effect sizes. A negative LDA score is not intrinsically “protective”; its meaning depends on the comparison order.

The parser accepts CSV or TSV and recognizes common alternate headers such as `marker_taxon`, `effect_size`, `pvalue`, `fdr`, and `experiment_type`.

### Clinical sample metadata

Normalize GMrepo project/run metadata to [data/gmrepo_samples_template.csv](data/gmrepo_samples_template.csv):

```bash
python populate_database.py load-samples path/to/project_samples.tsv
```

Required fields are `project_id` and `run_id`. Disease, MeSH ID, age, sex, country, and QC status are optional.

### Sample-level relative abundance

Convert abundance data to long format using [data/gmrepo_abundances_template.csv](data/gmrepo_abundances_template.csv), load samples first, then run:

```bash
python populate_database.py load-abundances path/to/abundances.tsv
```

Required fields are `run_id`, `scientific_name`, and `relative_abundance`.

### Sync GMrepo phenotype names and MeSH IDs

```bash
python populate_database.py sync-phenotypes
```

This uses GMrepo's documented `get_all_phenotypes` endpoint. Network access is only required for this command.

### Synchronize every curated GMrepo disease comparison

```bash
python populate_database.py sync-all-comparisons
```

This traverses every comparison and every project returned by GMrepo's current catalog,
stores project-level evidence in `gut_microbiome`, and mirrors the results into the legacy
`microbiome_dataset` schema. New taxa are added to
`microbiome_dataset.gut_microbiome_new`; the updatable
`microbiome_dataset.gut_microbiome_new_full` view exposes the same canonical table under
the requested full-table name. Aggregated disease links are written to
`microbiome_dataset.microbe_disease`.

The command is authoritative and idempotent: rows absent from the current GMrepo export
are pruned, while manually curated values in the taxon table are preserved.

### Backfill missing sample sex/age/BMI from NCBI

GMrepo's own export leaves `sex`/`age_years` blank for roughly 60% of samples because many
submitters never populated those BioSample fields. Where NCBI's SRA/BioSample records do carry
`host_sex`, `host_age`, or a BMI attribute, they can be pulled in directly:

```bash
python populate_database.py enrich-demographics
```

This only fills currently-NULL `sex`/`age_years`/`bmi` fields (`COALESCE`-style) and never
overwrites a curated GMrepo value. `age_years` parsing is unit-aware (days/weeks/months/years,
plus `90+`-style censored values) and rejects anything that converts to an implausible age rather
than trusting a mislabeled attribute. Pass `--limit N` to test against a small batch first. Every
run writes a full audit CSV (default `data/ncbi_biosample_demographics.csv`) recording the
resolved value and parsing note for every sample it touched, and logs to `ingestion_runs` like
every other loader. Not every gap is fillable — many source projects simply never submitted this
metadata to NCBI either, so `sex` and `age_years` will still have real gaps after running it.

## 4. Validate and query

```bash
python populate_database.py validate
python -m unittest discover -s tests -v
```

Example query:

```sql
SELECT scientific_name, disease_name, direction, effect_size, project_accession
FROM v_taxon_disease_summary
WHERE disease_name = 'Irritable Bowel Syndrome'
ORDER BY ABS(effect_size) DESC;
```

To list bacteria linked to more than one disease:

```sql
SELECT t.scientific_name, COUNT(DISTINCT a.disease_id) AS disease_count
FROM taxa t
JOIN taxon_disease_associations a ON a.taxon_id = t.id
GROUP BY t.id, t.scientific_name
HAVING COUNT(DISTINCT a.disease_id) > 1
ORDER BY disease_count DESC, t.scientific_name;
```

## Data notes

- MiMeDB taxonomic and metabolic text is normalized before matching. Species are matched case-insensitively by `(genus, species)`.
- The derivation rules follow the supplied workflow: facultative anaerobes and conflicting pathways become `mixed`; anaerobes default to `fermenter`; aerobes/microaerophiles default to `respirator`; high-specificity food-source terms take precedence.
- The included IBS marker CSV was transcribed from the supplied PDF. Confirm spellings and values against the live GMrepo project export before treating it as publication-ready evidence.
- GMrepo currently documents downloadable project/run data and REST endpoints for phenotype and abundance access. Marker imports are intentionally file-based so the pipeline does not depend on an undocumented web endpoint.
