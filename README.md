# Gut microbiome database population pipeline

A curated and normalized gut microbiome dataset integrating microbial taxonomy, metabolic traits, oxygen requirements, measured metabolite levels with their taxon associations, and disease associations, with a focus on IBS-related microbiome analysis.

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

### The two schemas: working and frozen

Two MySQL schemas exist in parallel on this server.

| schema | role |
|---|---|
| `gut_microbiome` | the **working** database. Everything points here: `.env`, the pre-push hook, the SQLite mirror build, every script. |
| `gut_microbiome_snapshot_20261018` | a **frozen** byte-for-byte copy taken 2026-10-18. Nothing in this repository references it. |

**It is a snapshot "as of" a date, not an uncurated baseline.** It contains every
repair made up to that point, so it is the fixed *before* for curation done
**after** 2026-10-18 and nothing earlier. A genuinely raw baseline is no longer
reconstructible by re-running the loaders, because several repairs now live
inside them -- `normalize_phylum` in the Muller loader,
`numeric_value_is_trustworthy` in `attribute_rows`, `resolve_effect_type`, the
derived `studies.data_type`, the protocol DOI in `metabolomics_method`. A fresh
load would produce something partly curated, which is worse than either end.

**For curation done before that date, the evidence is the audit CSVs, not a
schema diff.** They record per-row before and after values and name the organism
and the number: `data/bracket_repoint_unwind_audit.csv` carries 1,215 rows with
`stored_before`, and `taxa_merge_mirror_audit.csv`, `ncbi_rank_audit.csv`,
`trait_vocabulary_audit.csv`, `body_site_audit.csv` and `richness_audit.csv`
cover the rest.

**Why work continues in `gut_microbiome` rather than in the copy.** Repointing
`.env` and twenty-odd scripts at a new name is twenty chances to miss one, and it
would leave the frozen schema as the one every tool resolves to by default -- so
a single forgotten `--apply` would write to the reference. Freezing the copy
instead means no code path can reach it.

This inverts the original instruction, which was to do the cleanups *in* the
copy, so it is worth being explicit that the inversion is deliberate and was
confirmed on 2026-10-25 rather than drifted into. **Settled; do not re-plumb the
scripts at the snapshot.** The guarantee the instruction was after -- an
untouched reference to compare against -- is what the frozen copy provides, and
it provides it more safely this way round, because the reference is the schema
nothing addresses rather than the schema everything defaults to.

**That protection is naming, not enforcement.** The usual connection is `root`
with full grants, and MySQL has no per-schema read-only mode for a superuser, so
the snapshot is safe because nothing addresses it, not because the server refuses
writes. Treat `USE gut_microbiome_snapshot_20261018` as a read-only context by
convention.

Query it directly, or compare the two:

```sql
SELECT COUNT(*) FROM gut_microbiome_snapshot_20261018.sample_taxon_abundances;

SELECT 'working' AS schema_name, COUNT(*) AS taxa FROM gut_microbiome.taxa
UNION ALL
SELECT 'snapshot', COUNT(*) FROM gut_microbiome_snapshot_20261018.taxa;
```

`validate` runs against whatever `.env` names. To check the snapshot, point a
connection at it and call `validate_database` directly; 27 of the counters apply
there (the three that need `sample_attributes` aliases and the mirror-prefixed
ones do not), and all 27 read zero as of the copy.

To take another snapshot, dump and restore with the server's own `mysqldump` --
match its version to the server, since a newer client can emit syntax an older
server rejects. Put the credentials in a `0600` defaults file rather than on the
command line, where they would be visible in the process list:

```bash
install -m 600 /dev/null /tmp/snap.cnf
printf '[client]\nhost=127.0.0.1\nuser=root\npassword=%s\n' "$GUTDB_PASSWORD" > /tmp/snap.cnf
SNAP=gut_microbiome_snapshot_$(date +%Y%m%d)
/usr/local/mysql/bin/mysql --defaults-extra-file=/tmp/snap.cnf \
  -e "CREATE DATABASE \`$SNAP\` CHARACTER SET utf8mb4 COLLATE utf8mb4_0900_ai_ci;"
/usr/local/mysql/bin/mysqldump --defaults-extra-file=/tmp/snap.cnf \
  --single-transaction --routines --triggers --events --set-gtid-purged=OFF \
  gut_microbiome | /usr/local/mysql/bin/mysql --defaults-extra-file=/tmp/snap.cnf "$SNAP"
rm -f /tmp/snap.cnf
```

`--single-transaction` takes a consistent read without locking the working
database, which matters because more than one agent session writes to it. Verify
a new snapshot by object count, per-table row counts and the generated columns
before trusting it -- a copy with intact tables and broken views passes a row
count and fails everything else.

The SQLite mirror is **not** forked. It is a build artifact rebuilt from the CSVs
in `data/`, and the most frequent defect in this project's history is a repair
that reached only one store; two schemas plus two mirrors would double that
surface for no gain.

## 2. Populate the database

Load the two local taxon sources and the consolidated GMrepo association export
(including all four IBS projects):

```bash
python populate_database.py all \
  --associations data/gmrepo_all_disease_associations.csv
```

This command:

1. Creates `gut_microbiome` and all tables/views.
2. Loads `gut_microbiome_new_full.csv`.
3. Collapses MiMeDB strain rows to species and derives `energy_mode` and `primary_food_source`.
4. Fills missing taxon attributes without overwriting curated values.
5. Loads the consolidated GMrepo disease associations, including PRJNA705217
   and the other IBS projects.

To intentionally replace existing taxon attributes/statistics, add `--overwrite`.

## 3. Add more GMrepo projects

GMrepo's documented API exposes phenotype, run, and abundance endpoints, while project marker tables can be downloaded as TSV from the project/comparison page. Normalize a marker export to the columns below, then run:

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

`p_value` and `q_value` are read where a source reports them, under those names
or `pvalue`/`p` and `qvalue`/`fdr`/`q`. They are independent of `effect_size`
rather than derived from it: GMrepo's export carries an LDA score and no
significance, while a hand-curated marker file often carries the reverse.

The parser accepts CSV or TSV and recognizes common alternate headers such as
`marker_taxon`, `effect_size`, `experiment_type`, and `direction`.

### One GMrepo association input

Use `data/gmrepo_all_disease_associations.csv` for GMrepo association imports.
`project_description` is omitted from this file and future GMrepo association
exports; project-level descriptions can still be held once in `studies`.
The PRJNA705217 transcription and the additional IBS-project export were
consolidated on 2026-10-16. The main file now contains 255 IBS association
records across four projects: PRJEB37924 (153), PRJNA682378 (34), PRJNA637763
(13), and PRJNA705217 (55). Genus and species evidence remain separate.

The older CSVs are preserved under `data/archive/gmrepo_ibs_consolidation_2026-10-16/`
and are no longer loaded separately. Its README records the overlap, score
disagreements, and two unmatched transcribed taxon labels retained pending
taxonomy review. MySQL already held all of these records in
`taxon_disease_associations`; this file consolidation does not reload or alter
the database. `q5_ibs_markers.csv` is a derived query report, not another input.

### Normalized LDA scores

`taxon_disease_associations.lda_score_normalized` is a generated column using
the fixed **softsign** transformation:

```text
lda_score_normalized = ROUND(effect_size / (1 + ABS(effect_size)), 3)
```

It is populated only for `effect_type = 'LDA'`; missing scores and other effect
types return NULL. Existing `effect_size` values are retained as stored. New
inserts and score corrections automatically refresh the generated value. The
summary view also exposes `lda_score_normalized`.

Three decimal places retain more resolution on the smaller scale: -4 maps to
-0.800, 0 to 0.000, and +4 to +0.800. The unrounded transformation is strictly
increasing and bounded between -1 and +1; rounding may create ties and, for
extreme future scores, endpoints. SQLite and some clients may omit trailing
zeros when displaying numbers; this does not change their numeric value.

This is a nonlinear rescaling, not min-max or z-score standardization. It needs
no fitted dataset parameters, so adding studies does not change previous values
and the transform does not use held-out data. It compresses magnitude differences;
retain `effect_size` for original thresholds and analyses requiring that scale.
It is not a probability or a correction for study/method differences. Its sign
follows the source comparison order; use `direction` and the comparison groups
to interpret disease enrichment/depletion. Disease-labelled association scores
must not be treated as independently measured patient features for predicting
those same disease labels.

Apply to an existing MySQL database and, separately, its SQLite mirror:

```bash
python3 scripts/normalize_lda_scores.py --apply
python3 scripts/normalize_lda_scores.py --sqlite gut_microbiome.sqlite --apply
```

Omit `--apply` for a preview. The migration verifies every derived value and
checks that all original association fields are unchanged. New databases get
the generated column from their schema. This does not change the legacy
`microbiome_dataset` tables, which do not store project-level LDA scores.

For a model needing a learned scaler instead, fit it inside each training fold
and apply that fitted scaler to validation/test data, as described in
[scikit-learn's preprocessing guidance](https://scikit-learn.org/stable/modules/preprocessing.html).

### Clinical sample metadata

Normalize GMrepo project/run metadata to these columns:

```bash
python populate_database.py load-samples path/to/project_samples.tsv
```

Required: `project_id`, `run_id`. Optional: `sample_id`, `disease_name`, `mesh_id`,
`sex`, `age_years`, `bmi`, `country`, `qc_status`, `body_site`. Each also accepts
aliases — `project_accession` for `project_id`, `run_accession` for `run_id`,
`gender` for `sex`, `age` for `age_years`, `disease`/`phenotype` for
`disease_name`.

The GMrepo sample export uses the readable `disease_name` and `mesh_id` fields
instead of raw `phenotypes` JSON. `nr_phenotypes`, `longitude`, `latitude`, and
`data_quality` are also omitted. This affects the CSV, not the database schema.
The sample loader preserves a study's existing quality classification when the
input omits `data_quality`.

To refresh the sample CSV from the demographics and curated disease labels
already recovered in MySQL:

```bash
python3 scripts/clean_gmrepo_samples.py --apply
```

Omit `--apply` to preview. Exact run/project/sample identifiers are checked
before merging. The script keeps an original backup (including raw phenotype
JSON) and a field-level change audit in `data/archive/`. Missing measurements
remain blank when no source value exists; BMI 0 and invalid age sentinels are
removed. It does not modify MySQL. Existing disease labels are reconciled with
the database's current curated disease/MeSH pair.

Future `sync-gmrepo-samples` exports use this same reconciliation step, retain
enriched demographics, and write an `_enrichment_audit.csv` alongside the
export. This prevents a raw API refresh from losing the recovered metadata.
New ambiguous phenotype lists stay unlabelled until an actual sample label
is available; choosing a disease from a cohort's list is not sufficient.

`age_years` and `bmi` are checked on the way in: a BMI of 0 or above 200 and an
age below 0 or above 110 are stored as NULL, because GMrepo serves those as
missing-value sentinels. Age 0 is kept — it is a real birth-day sample in at
least one infant cohort.

### Sample-level relative abundance

Convert abundance data to the long format below, load samples first, then run:

```bash
python populate_database.py load-abundances path/to/abundances.tsv
```

Required: `run_id`, `scientific_name`, `relative_abundance`. Optional:
`taxonomic_rank`, `ncbi_tax_id`, `source_database`. A declared
`taxonomic_rank` of `genus` is honoured even when the name reads as a binomial,
which is how GMrepo labels entries like `[Bacteroides] pectinophilus`.

Genus and species are independent layers: each must sum to 1 within its own
rank for a given sample, and `validate` checks that.

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

The CKD supplement audit, `data/ckd_study_demographics_audit.csv`, keeps
`run_accession`, `supplement_id`, `library_name`, `join_basis`, `age_years`,
`sex`, and `bmi`. Duplicate `filled_*` values are omitted; the cohort series
is already encoded in `supplement_id`.
`join_basis` must be retained to describe the source identifier match. Following
source verification, only 57 HzCKD mappings remain active. The 273 HC and 159 CKD
unverified mappings were archived and their imported demographics withdrawn from
MySQL and the sample CSV. SQLite now contains only the verified demographics for
these runs. Samples and abundance data were preserved. The importer no longer
guesses numeric-prefix joins; see [CKD cleanup details](docs/ckd_demographic_cleanup.md).

Clearing those cells did not make the withdrawal durable. `sync_gmrepo_samples`
re-fetches the same demographics from the GMrepo API, and
`prepare_gmrepo_sample_export` fell back to the fetched value whenever the
database held nothing -- which is right for a cell that was never measured and
exactly wrong for one emptied on purpose. The two are indistinguishable in the
data, so the next sync of PRJNA562327 would have refilled all 432 withdrawn
sets, and `load_samples` would have written them back under
`COALESCE(VALUES(sex), sex)`. The export audit could not have caught it either:
the refilled value matches the source, so nothing reads as a change. Verified by
simulation before the fix -- a withdrawn run fed a fetched `female/54/23.4` came
back carrying it, with only the `54` to `54.0` reformat logged.

`demographic_withdrawals` is what makes the decision stick. It lists
(`run_accession`, `field`) pairs whose cells stay blank whatever a source
serves, is seeded by `scripts/seed_demographic_withdrawals.py` by re-deriving
the verified set from `results/ckd_link_verification/verification.json` rather
than transcribing a list, and lives in both stores from
`data/ckd_demographic_withdrawals.csv`. It covers 1,920 cells across 640 runs:
the 432 audited-but-unverified mappings plus the 208 project samples that had no
clinical mapping at all, which the original repair left open. The 57 verified
runs are untouched, and a blocked refill is written to the export audit so the
ledger never acts silently. `withdrawn_demographic_cell_refilled` counts a
stored value for a withdrawn cell in both stores, so a revert by any path is a
validate failure rather than a discovery; it was confirmed to report 1 against a
simulated refill, not merely 0 against clean data.

### "Succinate production" was in the description and in no column

This line has described the dataset as carrying **succinate production** since
the first commit. No column ever did: `git log -S succinate -- gutdb/schema.sql`
returns nothing, and succinate appeared in exactly one tracked file, that
sentence. Only 1 of 5,838 taxa mentions succinate in its MiMeDB `metabolism`
text, so there was nothing to read it from either.

There is **no source in this database for whether an organism produces
succinate**. MiMeDB carries organism-to-metabolite upstream, but the loaded
`gut_microbiome_new_full.csv` is twelve columns and does not. Asserting a
`produces_succinate` boolean would have been the same inference-as-measurement
that the CKD demographics were cleared for.

What the data does support is a measured association, and succinate is well
covered: `C00042_Succinate` has levels in 2,083 samples, every one of which also
carries an abundance profile, 440 of them the MARS IBS cohort.
`taxon_metabolite_associations` holds one Spearman correlation per (taxon,
metabolite, study), and `v_taxon_metabolite_evidence` combines them.

Three things decided the method, each measured rather than assumed:

- **Per study, never pooled.** Succinate levels are not comparable between
  cohorts -- the per-study means run from 0.06 to 3.0e+08, and two cohorts are
  stored centred so they carry negative values. Spearman inside one study is
  invariant to monotone rescaling, so the scale difference stops mattering once
  nothing crosses a study boundary.
- **CLR, not raw proportions.** This was the correction that mattered. Relative
  abundances sum to 1, so when a few genera rise with succinate the rest must
  fall. Raw proportions gave **463 negative against 126 positive** below FDR
  0.05 -- 79% negative, which is not a finding about succinate but about the
  simplex. The same test on CLR values gives **269 negative against 319
  positive**, 46%. The first implementation here used raw proportions and the
  skew is what exposed it.
- **FDR within each study**, over the taxa actually tested there, so the 1,992
  pairs across 11 cohorts do not hand back a page of discoveries that are the
  multiplicity.

The result is coherent: Enterobacterales up with succinate, Oscillospiraceae
butyrate producers down, which is the usual dysbiosis axis and matches succinate
being a dysbiosis marker. **Reading it as ~115 independent genus results would
be wrong.** Closely related genera co-vary, and of the genera agreeing in at
least two studies with none opposing, the positives are 22 Lachnospiraceae, 13
Enterobacteriaceae and 6 Peptostreptococcaceae while the negatives are 25
Oscillospiraceae -- a few family-level observations wearing genus-level counts.
A significant genus was tested in 5.1 of the 11 cohorts on average, because of
the prevalence floor, and 8 genera are significant in **both** directions; the
view reports those as contested rather than averaging them to nothing.

The whole of `FRANZOSA_IBD_2019` -- all 220 samples, the IBD cohort, where
succinate is most discussed -- is flagged low-confidence for this annotation
upstream. It is kept and flagged through
`all_studies_high_confidence`, not dropped, because which confidence floor is
acceptable is the caller's decision.

```bash
python3 scripts/derive_taxon_metabolite_associations.py            # preview
python3 scripts/derive_taxon_metabolite_associations.py --apply
python3 scripts/derive_taxon_metabolite_associations.py --kegg-id C00246 --apply
```

This layer is **MySQL only, and that is not a parity gap**: the mirror carries
no metabolite tables at all, so there is nothing for it to drift from -- unlike
the defects this schema keeps producing, where both stores held the data and one
was repaired. Its four counters are in `_MYSQL_ONLY_CHECKS` for that reason.
Half-adding it would leave `metabolite_id` dangling against a mirror with no
`metabolites`.

### The schema could not be built on a clean machine

`gutdb-cli init` failed on any database that did not already exist.
`v_attribute_concept_coverage` was defined around line 514 and selects from
`v_ml_ready_samples`, which was defined around line 821, so MySQL refused it
with `Table '<db>.v_ml_ready_samples' doesn't exist` and the build stopped
there. It reproduces on every earlier commit, so it had been true for a long
time.

**Nothing caught it because nothing could.** Every integrity counter in this
repository runs against a database that is already built, where all twelve views
exist; a view defined before the view it selects from passes all 57 of them. The
defect only exists during construction, which is the one moment nothing was
looking. Existing databases were unaffected, which is why it survived.

The fix is pure reordering -- the view and its comment block moved below
`v_ml_ready_samples`, verified by sorting both files and diffing, so the content
is identical and only the order changed. A fresh database now builds 16 tables
and 12 views, matching the live schema exactly, and a fresh SQLite mirror builds
10 tables and 8 views with every view resolving when queried.

`tests/test_schema_buildable.py` keeps it closed without needing a server: the
MySQL schema is checked by parsing it for views that select from views defined
later, and the SQLite schema is genuinely built in a temporary file and every
view queried, because SQLite accepts a view over a missing table at creation
time and only fails on use. The detector is itself pinned against a known-bad
and known-good pair, so a regex that silently stopped matching would fail rather
than pass, and it was confirmed to flag the pre-fix schema.

### What is not droppable, measured rather than assumed

A sweep for tables worth dropping found none. All 16 carry rows and all have at
least two readers in tracked code. A column-level sweep for zero-information
columns returned four, and three of them are deliberate:

- `demographic_withdrawals.reason` and `taxon_metabolite_associations.metabolite_id`
  hold one value only because there has been one withdrawal campaign and one
  metabolite so far. Both are keys; a second of either makes them vary.
- `taxa.human_pathogen` is `NULL` in 5,775 rows and `1` in 63, never `0`. That
  is not a dead column, it is the correct encoding: MiMeDB flags known pathogens
  -- *Clostridium perfringens*, *Staphylococcus aureus*, *Bacteroides fragilis* --
  and does not assert that anything is harmless, so absence is unknown rather
  than negative. Dropping it would discard 63 measured facts.
- `taxon_disease_adjudication.assoc_split_studies` is `0` in all 2,284 rows, and
  that zero is the finding. `scripts/adjudicate_conflicts.py` writes it and
  `tests/test_adjudication.py` asserts it reaching 1, so it is live and tested;
  no taxon-disease pair currently has a study that disagrees with itself.

### How far the succinate associations actually reach

394 of 5,838 taxa carry one, which is 6.7%, and the shortfall is data rather
than filtering: **every taxon that appears in the eleven metabolomics cohorts at
all is covered**, none were lost to the prevalence or support floors. The other
5,444 taxa have no metabolomics sample whatsoever. The ceiling is lower still
than it looks, because those cohorts carry only genus-level profiles, so of
2,148 genus-layer taxa only the 394 present in them can ever be reached this
way, and no species row can be. Extending coverage needs new metabolomics, not
new code; `--kegg-id` runs any of the other 1,247 metabolites over the same 394.

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

### Correcting the evidence view

`v_taxon_disease_evidence` reported **wrong directions** until 2026-10-10, and
anything that read `consensus_direction` from it inherited them. Two rules were
broken:

- **It pooled unlike contrasts.** Disease-vs-disease comparisons voted alongside
  disease-vs-Health ones, as though "enriched in UC relative to Crohn disease"
  answered a case/control question.
- **One study could vote twice.** `n_studies_enriched` and `n_studies_depleted`
  were two independent `COUNT(DISTINCT study_id)` expressions, so a study holding
  both an enriched and a depleted association counted on *both* sides.

The view now restricts the direction vote to case/control contrasts at one vote
per study, exposes `contrast_scope`, `n_studies_case_control`, `n_studies_split`
and `n_studies_other_contrast`, and returns `consensus_direction = NULL` for the
753 pairs whose associations are *all* disease-vs-disease. Replication depth
(`n_associations`, `n_studies`, `n_sources`) still counts every contrast, because
a disease-vs-disease result is a real result -- it just cannot answer a
case/control question.

What changed, over the 2,285 replicated pairs:

| change | pairs |
| --- | --- |
| `tied` → decided | 118 |
| decided → `tied` (a genuine tie, previously hidden) | 41 |
| **direction reversed outright** | **33** |
| → `NULL` (no case/control evidence exists) | 34 |
| unchanged | 2,059 |

All 226 changed pairs are `mixed` or `disease_vs_disease` scope; no pure
`vs_health` pair moved. The 33 reversals are **27 in ulcerative colitis and 6 in
Crohn disease** -- the diseases with the most disease-vs-disease comparisons --
and they are mostly butyrate-producing commensals whose depletion in IBD is
textbook: *Faecalibacterium*, *Lachnospira*, *Agathobacter rectalis*, *Dorea*,
*Ruminococcus*. The abundance matrix in this same database independently backs
the corrected direction over the old one **24 times to 7**, with 2 inconclusive.

Measured honestly, the second rule currently changes nothing: no study reports
both directions among its own vs-Health associations, so `n_studies_split` is 0
for all 11,472 pairs and every behaviour change traces to the contrast
restriction. It is kept as a guard, not a fix.

`v_taxon_specificity.depleted_fraction` had the same contrast-mixing defect and
is corrected the same way. Pooling both kinds moved it by more than 10
percentage points for 85 of the 821 taxa with three or more associations, and by
more than 25 points for 21 of them (worst case 0.40) -- enough to turn a general
dysbiosis marker into an apparently specific one. Its breadth columns stay
association-level on purpose: specificity asks how widely a taxon has been
reported, not how cohorts voted on one disease.

### Replication, not error: reading q7's disagreements

```bash
python scripts/replication_audit.py     # writes data/replication_audit.csv
```

`q7` flags 191 of 4,392 species-level pairs where the curated associations and
this database's abundance matrix point opposite ways. The obvious reading -- that
those rows are mistakes -- is wrong, and this script exists to say why. For every
per-study comparison it asks whether that study is one that **reported** the
association:

| pair group | cohorts that reported it | other cohorts |
| --- | --- | --- |
| abundance **agrees** (4,201 pairs) | 98.1% match | 75.1% match |
| abundance **disagrees** (191 pairs) | **91.9% match** | **18.3% match** |

The associations describe their own cohorts correctly in *both* groups. What
makes a pair a "disagreement" is that no other cohort reproduces it. So
`disagreement_kind` in q7 separates the two cases, and as of 2026-10-12 the
split is **190 `non_replicating`, 1 `no_own_cohort_abundance`, and zero
`contradicted_in_own_cohort`** -- not one association is contradicted by the
measurements of the study it is attributed to, which is a statement about the
loaders as much as about the biology.

*Clostridium perfringens* in Crohn disease is the pattern in miniature. It is
flagged depleted while the abundance says enriched 6-1 at a 15-fold difference,
which looks like a sign error until you check the source: the association comes
from PRJNA993675, and that study's own abundance reads 0.00165 in Crohn against
0.00461 in health. It is depleted *there*. Seven other cohorts measure the
opposite.

Most of the 191 are not even substantive. When the abundance disagrees its own
vote is typically a bare majority (median agreement 0.667) whereas agreeing
pairs are usually unanimous (3,459 of 4,201 at 0.9 or above), and the median
absolute log2 fold change is 0.297 against 1.593 across all pairs, on taxa
sitting near 0.04% relative abundance where a handful of reads flips the sign.
Two filters pull in opposite directions and both are reported, because quoting
either alone misleads: demanding more comparable studies **raises** the
disagreement rate (4.3% to 15.5%, since single-study pairs confirm themselves),
while demanding a real contradiction **lowers** it (8.8% to 1.0%). Ten pairs
survive an abundance majority of 0.8 with a two-fold difference over at least
two studies, and nine of those ten rest on a single association study.

**For modelling, this is a generalisation warning list rather than a cleanup
list.** A feature whose direction holds in one cohort and reverses in six others
is exactly what fails on held-out data.

### Adjudicating directional conflicts

2,285 (taxon, disease) pairs are reported by more than one study and the studies
do not always agree: 387 are exact ties, where `v_taxon_disease_evidence`'s
`consensus_direction` is a coin flip presented as a finding. `scripts/adjudicate_conflicts.py`
settles what it can against evidence that vote never uses -- this database's own
abundance matrix -- and records the rest as undecided.

```bash
python scripts/adjudicate_conflicts.py            # dry run, writes data/adjudication_audit.csv
python scripts/adjudicate_conflicts.py --apply     # writes taxon_disease_adjudication
```

For each study holding both arms, the taxon's mean relative abundance in that
study's cases is compared with its mean in that study's controls, giving one
direction per study. Directions are then counted the same way association votes
are. Two details decide whether the number means anything:

- **Absence is zero, not missing.** An undetected taxon has no row in
  `sample_taxon_abundances`, and 97.3% of (sample, taxon) cells at species rank
  are such absences. Averaging only the rows that exist compares the samples
  where a taxon was abundant against the samples where it was abundant. The
  denominator is every sample in the arm profiled at that taxon's rank. **This
  is the rule that matters**: on the 4,046 pairs every method can score,
  counting absence as zero moves concordance with the curated associations from
  65.8% to 95.5%, while the choice of pooling moves it by -0.2.
- **Within-study, never pooled.** A direction computed inside one study compares
  like with like -- same protocol, same sequencing run, same recruitment -- and
  only directions are aggregated, so one large cohort cannot outweigh several
  smaller ones. It is the correct comparison, but it is worth only -0.2 points
  on the concordance figure (95.5% pooled against 95.3% within study), and an
  earlier version of this README claimed the 62% was what pooling costs. It was
  not; it was what skipping absences costs.

**Only case/control contrasts vote.** 2,340 of the 16,889 associations compare
one disease against another rather than against Health. "Enriched in UC relative
to Crohn disease" and "enriched in UC relative to health" are different claims,
and a case/control model needs the second. The association vote is therefore
built from vs-Health comparisons only, one vote per study, and
`assoc_contrast_scope` records which contrasts a pair actually had.

This was originally a workaround for a defect in `v_taxon_disease_evidence`,
which pooled both kinds and let a study vote on both sides at once -- it reported
*Faecalibacterium prausnitzii* as enriched in ulcerative colitis, reversing one
of the field's most replicated findings. **That view is now fixed** (see
[Correcting the evidence view](#correcting-the-evidence-view)), so the rule here
is no longer a workaround but a second, independent implementation: the
adjudicator recomputes the vote in Python and compares it against the view's SQL
on every pair, printing `view cross-check: ... agrees on all 2285 pairs` when
they match. Either one drifting shows up immediately.

Current outcome over the 2,285 pairs:

| basis | pairs |
| --- | --- |
| `association_and_abundance_agree` | 1,620 |
| `association_consensus` (abundance too thin to speak) | 186 |
| `abundance_majority` (association tied, abundance decided) | 145 |
| `association_consensus_abundance_disagrees` | 126 |
| `abundance_only_no_case_control_association` | 22 |
| `unresolved_*` | 186 |

187 of the 387 dead-tied pairs gained a direction. A NULL `verdict` is a result,
not a gap: *Prevotella* in ulcerative colitis stays 3-3 after abundance
adjudication, which is real cohort heterogeneity and should not be handed to a
model as a fact. Nothing in `taxon_disease_associations` is modified, so the
verdict layer can be recomputed, audited, or ignored.

### Muller 2022 microbiome-metabolome collection

[Muller, Algavi & Borenstein 2022](https://doi.org/10.1038/s41522-022-00345-5)
curated 14 human faecal cohorts -- 2,900 samples from 1,849 subjects -- with
paired microbiome **and metabolome** profiles. It matters here for three reasons:
it adds a feature class this database has none of, it is a second *sample-level*
source (98.3% of samples come from GMrepo, and MicrobiomeHD contributes 20
studies but no samples), and it arrives as plain TSV rather than through
Bioconductor the way curatedMetagenomicData's profiles do.

```bash
python scripts/fetch_muller2022.py      # ~140MB into data/muller2022/ (gitignored)
python scripts/map_gtdb_genera.py       # writes data/gtdb_ncbi_genus_map.csv
python scripts/load_muller2022.py       # dry run
python scripts/load_muller2022.py --apply
python scripts/adjudicate_conflicts.py --apply   # new abundance evidence
```

Loaded: **2,464 samples, 491,507 abundance rows, 514,359 metabolite levels over
1,222 compounds**, across 13 cohorts. GMrepo's share of samples falls from 98.3%
to **91.7%**, the deduplicated ML-ready cohort rises from 11,440 to **13,328**,
and `samples.subject_id` goes from 479 populated rows to 2,943. Re-adjudicating
afterwards gave **401 pairs more abundance evidence and flipped no direction at
all** -- 13 cohorts processed by a different pipeline (kraken2/GTDB rather than
GMrepo) reversed nothing, which is the strongest check on the adjudicator so
far.

Sized honestly against the goal of filtering down to ML-ready samples: 2,900
samples is +8.5% on 34,229, but after excluding `iHMP_IBDMDB_2019` (already here
as PRJNA398089) and the two cohorts with no disease contrast, 2,077 case/control
samples remain, and those collapse to **1,581 independent subjects** because 36%
of samples are repeats. Against the deduplicated ML-ready cohort of 11,440 that
is +13.8%. `Subject` is present in all 14 cohorts, which is better than the
GMrepo bulk manages -- `samples.subject_id` is populated for only 479 of 34,229
rows today.

**Taxonomy is GTDB, and that is the whole difficulty.** Columns are full lineage
strings, and the genus field is frequently a genome-bin accession, so the 12,263
distinct labels are nowhere near 12,263 genera. `scripts/map_gtdb_genera.py`
decides each label against the local NCBI taxdump and writes a reviewable row
rather than resolving anything, because handing these strings to `resolve_name`
is how `Ruminococcus2` and `Escherichia/shigella` were minted. By abundance mass:

| decision | labels | mass |
| --- | --- | --- |
| `map` (NCBI genus already in `taxa`) | 2,075 | 88.84% |
| `map_needs_taxon_row` (NCBI knows it, this DB does not) | 1,970 | 2.95% |
| `reject_mag_bin` (`g__UBA1775`, `g__JABGPL01`) | 8,022 | 4.89% |
| `reject_unclassified` (lineage stops above genus) | 126 | 3.24% |
| `reject_gtdb_coined` (NCBI has never heard of it) | 70 | 0.08% |

Four consequences for any loader, each measured rather than assumed:

- **It must SUM, not insert per label.** 4,045 mapped labels collapse onto 3,333
  NCBI genera; 339 targets receive more than one label (`Clostridium` receives
  40, `Enterococcus` 11, `Ruminococcus` 7, `Bacteroides` 6) because GTDB splits
  polyphyletic genera and marks fragments with letter suffixes. **52.19% of all
  abundance mass sits in multi-label targets**, so last-write-wins would silently
  discard half the data. Summing the fragments is also what makes the result
  comparable with GMrepo's NCBI-profiled genera in the first place.
- **Renormalize after dropping.** 8.21% of mass is rejected, so each sample needs
  renormalizing to restore the per-rank sum-to-1 invariant -- the same treatment
  `unclassified_fraction` already receives.
- **Apply a mass floor before creating taxa rows.** Mapping everything would add
  1,853 genus rows, and they start `Abyssibacter`, `Acaryochloris`,
  `Acetohalobium` -- marine and environmental kraken2 noise, not gut flora. A
  floor of 0.001% of total mass keeps 395 genera and **99.81% of mappable mass**
  for only **126 new rows**. Cohort prevalence adds nothing; mass is the
  discriminating filter.
- **`Study.Group` needs a hand-written MeSH mapping per cohort**, in
  `PHENOTYPES`. `1`=CRC in SINHA_CRC_2016; `D`/`H`/`C`=IBS-D/Healthy/IBS-C in
  MARS_IBS_2020, both IBS arms collapsing to one MeSH disease because the schema
  has no subtype field; `MP`=polyps in YACHIDA_CRC_2019. `Adenoma` (D000236) is
  created on load.

Four kinds of double-counting are refused, three of which only the data
revealed:

- **`iHMP_IBDMDB_2019` is skipped entirely** -- already here as PRJNA398089 from
  GMrepo, so its 382 samples are not new.
- **53 samples appear in both YACHIDA_CRC_2019 and ERAWIJANTARI_GASTRIC_CANCER_2020**,
  flagged by the collection itself. All 53 are Healthy controls. They are dropped
  from YACHIDA, which keeps 74 and a working contrast; dropping them from
  ERAWIJANTARI would leave it with one. Their subject ids differ only by a
  suffix (`10025` against `10025.Healthy`), **so a dedup keyed on `subject_id`
  would not have caught them** -- the suffix is stripped on load.
- **Explicit zeros are dropped, not stored.** The GTDB tables are dense and write
  `0.0` for every absent genus, which would have added 233,965 rows meaning
  "measured as absent" to a table where absence is the absence of a row. That
  distinction is load bearing: the adjudicator counts a missing row as zero, so
  stored zeros would leave means unchanged but make every prevalence count wrong.
- **Metabolite columns sharing one identifier within a cohort are dropped and
  reported** (221 columns). Two causes are indistinguishable from the data: the
  same analyte measured twice (`HILIC_NEG_alanine` beside `HILIC_POS_alanine`),
  and an upstream annotation error -- ERAWIJANTARI gives both `C00197_3PG` and
  `C00661_G3P` the identifier HMDB0000807, though 3-phosphoglycerate and
  glycerol-3-phosphate are different molecules. Averaging would invent a value
  across two compounds; last-write-wins would discard a measurement silently.
  The check runs on the resolved `metabolites.id`, not on the identifier string,
  because a column carrying only an HMDB id and one carrying only a KEGG id can
  resolve to the same row.

ERAWIJANTARI's `Gastrectomy` arm (42 samples) gets its own phenotype,
`Gastrectomy` (MeSH D005743), rather than Stomach Neoplasms.
[Erawijantari et al. 2020](https://doi.org/10.1136/gutjnl-2019-319188) studied
"participants with a history of gastrectomy for gastric cancer" and asked what
the *surgery* did to the microbiome; `Surgery_Type` confirms subtotal and total
gastrectomies, so the tumour is resected and calling these samples Stomach
Neoplasms would assert an active cancer that is not there. The paper is indexed
under both MeSH terms and the exposure under study is the procedure. `diseases`
already functions as a phenotype table -- `Health` is not a disease either.

YACHIDA's six group labels collapse the paper's **nine** clinical categories:

| paper's group | collection label | loaded as |
| --- | --- | --- |
| (1) normal, no remarkable findings | `Healthy` | Health |
| (2) a few polyps, up to two small (<5 mm) | `Healthy` | Health |
| (3) MP: >=3 (mostly >=5) polypoid adenomas, low-grade dysplasia, conventional type only, **not** serrated | `MP` | Adenoma |
| (4) intramucosal carcinoma: polypoid adenoma(s) with **high-grade** dysplasia, stage 0/pTis | `Stage_0` | Colorectal Neoplasms |
| (5)-(6) stage I, stage II CRC | `Stage_I_II` | Colorectal Neoplasms |
| (7)-(8) stage III, stage IV CRC (UICC 8th ed.) | `Stage_III_IV` | Colorectal Neoplasms |
| (9) normal with a history of colorectal surgery | `HS` | Colectomy (D003082) |

`Stage_0` is a carcinoma in situ rather than an adenoma, so it joins the neoplasm
arm while `MP`, which is low-grade-dysplasia adenoma, does not. `HS` is group
(9), and it is **not** the 28 stage I-III patients the paper sampled before and
after surgery -- those are a separate subset, and no subject in this collection
has more than one sample. It maps to `Colectomy`, the closest operative MeSH
term and the same `E04.210` branch as `Gastrectomy`; the source says "colorectal
surgery" without naming the procedure, so rectal resections fall under this
label too.

**Caveat on the controls, which outlives this cohort.** The paper defines groups
(1) *and* (2) as the healthy controls, so `Healthy` includes subjects carrying up
to two polyps under 5 mm. Nothing in the collection's metadata separates them --
`Stage` and `Tumor location` are `-` across the whole arm -- so it cannot be
undone at load time and is a property of every contrast drawn against these
controls. It reaches further than YACHIDA: 53 of those 127 are also the **entire**
control arm of ERAWIJANTARI_GASTRIC_CANCER_2020.

**Neither post-surgical arm becomes `Health`.** A resected gut is precisely the
confounder [Erawijantari 2020](https://doi.org/10.1136/gutjnl-2019-319188) was
written to document, so putting these 72 samples in the control pool would bias
every contrast drawn against it. They sit beside their cohort's controls instead:
YACHIDA holds 74 Health, 150 Colorectal Neoplasms, 40 Adenoma and 30 Colectomy;
ERAWIJANTARI holds 54 Health and 42 Gastrectomy. No sample in the database now
lacks a phenotype.

### Making the abundances modellable: the CLR transform

```bash
python scripts/export_ml_matrix.py "Crohn Disease"
python scripts/export_ml_matrix.py "Colorectal Neoplasms" --rank=species --prevalence=0.2
```

Relative abundances are compositional. They carry only ratios, they sum to 1
within a rank, and two taxa cannot both rise -- so an increase in one forces an
apparent fall in the others, and a model trained on raw proportions learns that
artefact. `gutdb/compositional.py` provides the centred log-ratio transform that
removes it, and `scripts/export_ml_matrix.py` writes a matrix for one disease
against Health.

**Nothing is stored in the database, deliberately.** Two measurements decide
that:

- **CLR is not invariant to the part set.** A three-part composition transforms
  to `[0.476, -0.035, -0.441]`; add a fourth part at 0.1% and the same three
  become `[1.910, 1.400, 0.994]`. A stored column would be wrong the moment
  anyone filtered the feature set, which is the first thing a modeller does.
- **The matrix is sparse enough that the transform would be mostly pseudocount.**
  A genus profile holds roughly 50 of 2,076 taxa, so 99% of the log-sum would
  come from replaced zeros and the geometric mean would be a function of the
  replacement constant rather than of the sample. Filtering to prevalent taxa
  first is not tidying; it is what makes the number mean anything.

So the order is load-bearing: **select samples, filter by prevalence,
renormalize, replace zeros, transform.** Renormalizing after filtering is what
keeps each row a composition. Zeros are replaced multiplicatively
([Martín-Fernández et al. 2003](https://doi.org/10.1023/A:1023866030544)) rather
than by adding a constant, because adding one changes the ratios between the
parts that were actually measured -- the only information a composition carries.
`delta` defaults to half the smallest observed value **in that sample**, since a
zero means "below this sample's detection limit" and the limit is not the same
for every sample. Eleven tests assert the properties rather than remembered
numbers: rows sum to zero, the transform is invariant to rescaling a row,
observed ratios survive zero replacement, and a naive pseudocount demonstrably
breaks them.

A Crohn export gives 6,968 samples by 121 genus features, with the prevalence
floor retaining 98.7% of the median sample's mass.

**The export prints a warning it is important not to skip.** The rows pool many
studies, and study is confounded with the label: of the 105 studies behind that
Crohn matrix, **96 contribute only one arm**. A model given a random train/test
split will partly learn which study a sample came from and report a flattering
score. `project_accession` is written into every matrix so cross-validation can
be grouped on it.

### What the cohort can actually support: leave-one-study-out

For a paired no-BMI baseline versus adding BMI **without requiring complete
sex/age metadata**, use `python3 scripts/compare_bmi_models.py`. The new workflow
fits prevalence filtering and demographic preprocessing inside each training
fold, reads either store without changing it, and records snapshot-specific
results. See [BMI comparison](docs/bmi_model_comparison.md). The older audit
below uses the complete-demographics cohort and a different experiment design.

```bash
python scripts/loso_audit.py      # writes data/loso_audit.csv
```

A gut-microbiome classifier is usually scored with a random train/test split,
and in this database that is close to meaningless. Study membership is
confounded with the label -- of the 105 studies behind the Crohn cohort, 96
contribute a single arm -- so a random split puts samples from one study on both
sides of it and the model is rewarded for recognising the cohort. Leave-one-
study-out holds out a whole study, which is what using the model on a new patient
would actually ask of it.

**Only 12 of 65 diseases can be scored at all**, needing three studies with at
least five cases and five controls each. 20 diseases have two such studies, and
none has ten. Plain L2 logistic regression on CLR genus features:

| disease | n | folds | random split | LOSO | LOSO range | optimism |
| --- | --- | --- | --- | --- | --- | --- |
| COVID-19 | 155 | 3 | 0.899 | **0.506** | 0.36–0.71 | +0.392 |
| Irritable Bowel Syndrome | 643 | 3 | 0.926 | **0.577** | 0.52–0.61 | +0.349 |
| Colitis, Ulcerative | 1,339 | 8 | 0.809 | 0.585 | 0.40–0.79 | +0.224 |
| Diabetes Mellitus, Type 1 | 278 | 4 | 0.735 | 0.517 | 0.43–0.58 | +0.218 |
| Multiple Sclerosis | 362 | 4 | 0.729 | 0.555 | 0.49–0.60 | +0.174 |
| Crohn Disease | 1,983 | 9 | 0.825 | **0.687** | 0.28–0.90 | +0.138 |
| Diabetes Mellitus, Type 2 | 226 | 5 | 0.636 | 0.520 | 0.43–0.73 | +0.115 |
| Colorectal Neoplasms | 712 | 7 | 0.731 | 0.647 | 0.42–0.85 | +0.083 |
| Autism Spectrum Disorder | 393 | 5 | 0.595 | 0.515 | 0.45–0.61 | +0.081 |
| Parkinson Disease | 358 | 4 | 0.533 | 0.463 | 0.31–0.52 | +0.070 |
| Depressive Disorder, Major | 246 | 3 | 0.554 | 0.534 | 0.52–0.56 | +0.021 |
| Fatigue Syndrome, Chronic | 1,014 | 4 | 0.748 | **0.743** | 0.69–0.80 | +0.005 |

Mean random-split AUC **0.727**, mean LOSO **0.571**, mean optimism **+0.156**.
IBS would be reported at 0.93 and COVID-19 at 0.90 on a random split; both are
chance-level across cohorts. Only two diseases clear a LOSO AUC of 0.65 --
Crohn disease and chronic fatigue -- and one of those needs a caveat: all four
chronic-fatigue studies are shotgun, so its folds never cross a sequencing
platform, while Crohn's nine mix 16S and shotgun and every fold is a platform
transfer as well. That is most of why Crohn's range is 0.28 to 0.90.

Read the ranges before the means. With three to nine folds these estimates are
noisy, and nine of the twelve diseases have at least one fold at or below chance.
Diseases scoring across platforms lose more on average (+0.175 against +0.060),
but there are only two single-platform diseases to compare against, and the
all-16S one still loses 0.115 -- the study is a confounder within a platform, not
only between them.

**The practical consequence:** the modelling-ready cohort is 13,328 samples, but
the number of *independent cohorts* is what bounds an honest estimate, and for
most diseases that is three to five. Reporting a random-split AUC from this
database overstates what the model would do on a new cohort by about 0.16 on
average, and by more than 0.3 for the two diseases that look most impressive.

### The 269 one-genus profiles are faithful, and still unusable

A sweep for degenerate data finds 269 (sample, rank) groups holding exactly one
taxon at relative abundance 1.0, all genus-rank and all from GMrepo. The obvious
reading is a loading bug that dropped the rest of the profile. **It is not.**
Checked against `data/gmrepo_all_abundances.csv`, the source export gives
**exactly one row for all 210 of them it contains** -- the database is faithful,
and the degeneracy is upstream.

What they are: samples where classification largely failed. Their mean
`unclassified_fraction` is **0.64** against **0.20** database-wide, 190 of the
269 are above 0.5 and 63 above 0.9. GMrepo assigned the one genus it could and
the renormalization to sum-to-1 made it 100%. All 269 carry `qc_status = '1'`,
so GMrepo's own QC does not flag them. The lone survivor is often *Bivalvulida*
(82 of 269), the myxozoan parasite, then *Bifidobacterium* (41) and
*Methanobrevibacter* (24).

So nothing is deleted -- the rows record what the source says. The problem is
that they reach a model. `ml_ready` asks only whether `n_taxa_profiled` exceeds
zero, so **57 one-genus samples sit inside the deduplicated ML cohort**, and a
one-part composition has a **CLR of identically zero**: the part *is* the
geometric mean, so the value carries no information at all regardless of the
abundance.

`clr_matrix` now takes `min_taxa_present`, defaulting to **2** -- a mathematical
floor rather than a quality judgement, since the transform is undefined below it.
The Crohn export drops 50 rows to it and records the count in its provenance
JSON. Raising it further is a modelling choice and is left to the caller; for
scale, the median sample holds 62 genera and the 5th percentile 9.

### The 344 two-taxon profiles are a different problem wearing the same clothes

Widening the sweep to profiles holding **two taxa or fewer** gives 613: the 269
single-genus ones, 316 with two genera, 28 with two species. Every one of the
612 that appears in the GMrepo CSVs matches the source **row for row**, and
GMrepo's own abundances for them sum to exactly **1.0** -- so these are reported
upstream as complete profiles, not truncated ones. (Three apparent mismatches in
PRJNA769052 are the documented genus-layer repair: `[Clostridium] methylpentosum`
and `[Bacteroides] pectinophilus` are binomials GMrepo files under genus, and
`scripts/fix_genus_layer_purity.py` repointed them to the parent genus the sample
already held.)

But the explanation that fits the 269 does **not** fit these. The two-taxon
cohort's mean `unclassified_fraction` is **0.27** against 0.20 database-wide --
close to normal. Classification did not fail. What is actually going on is that
**depth is a property of the study, not of the sample**:

| study | samples | median genera | max genera |
|---|---|---|---|
| PRJNA1031953 | 101 | 1 | 4 |
| PRJNA832282 | 82 | 2 | 6 |
| PRJDB6133 | 58 | 3 | 7 |
| PRJNA422961 | 118 | 3 | 9 |
| PRJNA923180 | 27 | 3 | 7 |
| PRJNA889687 | 130 | 3 | 12 |

Six studies, **516 profiled samples, in which no sample reaches 15 genera**. A
per-sample filter cannot express this: inside these studies a thin profile is not
an outlier, it is the cohort.

**The harm is a confound that reads as a result.** When a disease is carried by
one thin study while its controls are pooled from a hundred rich ones, profile
depth predicts the label by itself. Neuromyelitis Optica is the clean case --
its whole case arm is PRJNA422961 at a median of 3 genera against controls at 45:

| predictor | AUC |
|---|---|
| 130 CLR genus features, 5-fold random split | **0.995** |
| profile richness alone, one integer, pooled | **0.982** |
| profile richness alone, **within study** | **0.497** |

The classifier is a richness detector. And the within-study figure is the one
that names the cause: PRJNA422961's own controls are just as thin as its cases,
so the gap is **pooling**, not a thinner gut. Eight cohorts score 0.75 or better
on richness alone, and **six of those have a single-study case arm** -- which is
exactly the shape `loso_audit.py` declines to score, because leave-one-study-out
needs a held-out study carrying both arms. Their random-split numbers were never
challenged by anything.

Two cases read differently, and the audit does not collapse them:

- **PRJNA788674** (Caliciviridae) holds 31 cases at 15 genera against its own 25
  controls at 30, and the within-study AUC stays at **0.873**. A collapse of gut
  diversity in acute norovirus is documented, so a within-study gap is evidence.
- **PRJNA769052** (COVID-19) holds 37 cases at a median of **2** genera against
  its own **50 controls at 43**, same platform, same project. No virus does that.
  This is the one confound shape nothing removes -- a per-sample filter cannot
  see it, and holding the study out takes both arms with it.

What changed, since again nothing is deleted:

- **`v_ml_ready_samples`** gains `n_genera_profiled`, `n_species_profiled`,
  `profile_richness` and `has_transformable_features`. The old advice to "filter
  on `n_taxa_profiled`" was **broken**: that column counts rows across every rank
  at once, so for the 12,566 dual-rank samples it is a sum of two richnesses and
  not the richness of either, and **44 samples are thin at one rank while the sum
  looks healthy**.
- **`v_ml_cohort_summary`** gains `max_genera`, `mean_genera`,
  `n_thin_genus_profiles` and `feature_depth_status`, whose
  `genus_depth_degenerate` value names the six studies above.
- **`gutdb/diagnostics.py`** computes the pooled and within-study richness AUC
  and the per-study depth ratio.
- **`scripts/export_ml_matrix.py`** reports all of it before writing, and records
  it in the provenance JSON. The Crohn export is unaffected in substance (0.570
  pooled, 0.628 within study) and prints no warning.
- **`scripts/richness_audit.py`** writes `data/richness_audit.csv` and
  `data/richness_arm_imbalance.csv`.

### A KEGG id is a lipid class, and it was being used as an identity

One metabolite held no levels at all —
`1-palmitoyl-2-linoleoyl-GPC (16:0/18:2)`, recorded but unmeasurable. Tracing it
found a correctness bug and a constraint that enforced it.

`ensure_metabolite_rows` resolved a compound by HMDB id and, **if that missed,
fell back to KEGG**. But a KEGG compound id is often a *lipid class* rather than
a molecule:

| KEGG id | what it is | distinct HMDB species in the maps |
|---|---|---|
| `C00422` | triacylglycerol | 28 — **but see the caveat below** |
| `C00157` | phosphatidylcholine | 16 |
| `C04230` | lysophosphatidylcholine | 15 |
| `C00165` | diacylglycerol | 10 |
| `C02530` | cholesterol ester | 8 |

So the second phosphatidylcholine to arrive missed on its own HMDB id, fell back
to `C00157`, found the **first one's row**, and silently reused it — two
different molecules on one row, with the loser's levels colliding on the
`(sample_id, metabolite_id)` key or dropped. `UNIQUE (kegg_id)` then made the
alternative impossible: only one member of each family could exist at all.

Two changes: the KEGG lookup now runs **only** for a compound with no HMDB id of
its own, and `kegg_id` is a non-unique annotation (`ix_metabolite_kegg`) rather
than an identity, with `hmdb_id` carrying identity as before.

**What the reload actually recovered**, measured against a snapshot pinned
beforehand:

| | before | after |
|---|---|---|
| metabolites | 1,222 | **1,248** |
| metabolite levels | 514,359 | **522,707** |
| metabolites with no levels | 1 | **0** |
| samples / abundance rows / attributes | — | **unchanged** |

**50 rows now sit in a shared-KEGG group**, across 24 KEGG ids — molecules that
could not all have existed under the old constraint. `C00157` holds all three of
its species (`HMDB0007972`, `HMDB0007973`, `HMDB0008038`) with 720 levels
between them, and `1-palmitoyl-2-linoleoyl-GPC` is no longer a compound with no
measurements.

**Two numbers above were overstated before the reload, and the reload is what
showed it.** `C00422`'s 28 triacylglycerol species live *only* in
`iHMP_IBDMDB_2019`, which is in `SKIP_COHORTS` — so that family was never in the
database and the reload could not recover it. Using it as the headline example
was wrong. And the earlier claim that "at least 109" compounds were blocked by
the constraint was an upper bound computed across *all* map files, including that
skipped cohort; the realized figure is 26 new rows. The 202-of-1,376 absence
figure stands, but most of it is the skipped cohort and the 216 within-cohort
identity collisions, not the KEGG constraint.

One thing this is *not*: the loader's "216 columns dropped for sharing an
identifier within a cohort" is a separate and deliberate behaviour, keyed on
`identity` (which prefers HMDB), covering two LC modes of one analyte and
genuine upstream annotation errors — ERAWIJANTARI labels both
3-phosphoglycerate and glycerol-3-phosphate `HMDB0000807`. Distinct HMDB species
sharing a KEGG id have distinct identities and never collided there. That
comment's count was stale at 124 and is now 216.

### Two store divergences left after everything else

Re-auditing every column this sweep touched, as value *sets* in both stores,
separated the real problems from the noise. Most differences are **membership** —
the mirror holds a different subset of taxa, so `Caryophyllales` and `Amelidae`
appear only there and `Viruses` only in MySQL. Two were not.

**`studies.data_type` in the mirror still held the raw source vocabulary for
every study with no samples**: `AMPLICON` (59 studies), `16S rRNA` (28),
`shotgun metagenomics` (9). The derivation pass only touched studies that *have*
samples, since a sample-less study has nothing to derive from — but normalizing
a stored value is not a derivation, and skipping it meant
`WHERE data_type = '16S'` matched **122 studies in the mirror against 236 in
MySQL**. Those 96 now fold through `normalize_platform`, which leaves `mixed` and
`meta-analysis` alone so nothing gains a platform it has no evidence for. The
mirror reads `16S` (209), `mNGS` (101), `meta-analysis` (2) — one vocabulary, and
no `mixed` because none of *its* studies span both platforms.

**The `Candidatus` prefix differed across stores for the same taxon.** Taxid
456826 carried `Candidatus Cloacimonadia` / `Cloacimonadales` /
`Cloacimonadaceae` in MySQL and the bare forms in the mirror — 5 values across
class, order and family. Neither store held both forms *internally*, so this was
invisible to any single-store check. `Candidatus` is how NCBI writes the
designation of an uncultivated taxon, so the fold keeps it, matching the
`Saccharimonadales` and `Obscuribacteraceae` entries already in the map.

Prefix splits remaining: **0** across all three ranks. What is left between the
stores is membership only.

### The class, order and family layers had the Firmicutes/Bacillota split too

`PHYLUM_SYNONYMS` solved this for one rank. The other three never got it:
`backfill_lineage.py` passes **only `phylum`** through `normalize_phylum` and
writes class, order and family verbatim from the dump. So a row filled from NCBI
carries the current scientific name while one loaded earlier keeps the old one,
and they sit side by side. Measured with the **tax ID identity test** — two
stored names resolving to the *same* NCBI node are synonyms of one taxon, never
a reclassification:

| rank | nodes split | rows |
|---|---|---|
| `class_name` | 1 | 409 |
| `order_name` | **10** | **1,500** |
| `family` | **13** | 557 |

`GROUP BY order_name` reported Clostridiales as two orders — 413 rows under
`Eubacteriales`, 127 under `Clostridiales`. Others: `Bacillales`/`Caryophanales`
(122/167), `Rhizobiales`/`Hyphomicrobiales` (36/113),
`Ruminococcaceae`/`Oscillospiraceae` (32/161).

**The target is the majority stored name, not NCBI's scientific name**, and that
choice mattered. Folding onto the scientific name would have:

- **re-introduced the nomenclatural code this schema strips.**
  `Saccharimonadales (SeqCode)`, `Saccharimonadaceae (SeqCode)` and
  `Obscuribacteraceae (SeqCode)` are NCBI's scientific names for three of these
  nodes, and writing them back **would have failed the pre-push check** I fixed
  yesterday;
- moved 107 `Streptomycetales` rows onto `Kitasatosporales`, which NCBI prefers
  and almost nobody uses — the same reason `PHYLUM_SYNONYMS` keeps *Firmicutes*
  over *Bacillota*.

The map is **explicit rather than computed at run time**, and that is load-bearing
for a reason this sweep keeps running into: the mirror's majorities differ from
MySQL's (`Bacillales` 183 vs `Caryophanales` 87 there), so a per-store majority
rule would have picked different winners in each store and **re-created the
two-store divergence** the rest of this work removed. One fixed map, both stores,
25 folds, 0 un-folded synonyms left in either.

### Filling the lineage from the row's own genus

143 taxa had no phylum. `backfill_lineage.py` resolves from the NCBI dump, which
cannot help most of them — 91 have no tax ID at all. But **62 did not need NCBI**,
because the answer was already in the database: rows sharing their `genus_key`
carried a phylum and they did not. An internal inconsistency, not missing
knowledge, and not marginal — `Waltera` has 7,944 abundance rows,
*Merdimmobilis hominis* 3,195, *Faecalibacterium butyricigenerans* 3,169 beside
nine `faecalibacterium` peers all saying Firmicutes.

**The guard is the point.** A `genus_key` is derived from a name's first token, so
it can group organisms that are not relatives. Two keys are ambiguous, and
checking them was worth it:

- **`bacteroides`** carries both Bacteroidetes (77 rows) and Firmicutes (1) — and
  **the Firmicutes one is right.** NCBI writes it `[Bacteroides] pectinophilus`,
  brackets meaning the name is misapplied; the organism is a *Clostridium*
  relative. One name whose genus token contradicts its lineage.
- **`unclassified`** carries Proteobacteria and Firmicutes across three MiMeDB
  rows named "Unclassified *&lt;family&gt;* bacterium …", where the first token is
  the word *Unclassified*. Both lineages are correct for their families; the key
  is simply meaningless there.

Neither has a phylum-less row waiting on it, so nothing was skipped — the guard
exists so that when one arrives it is skipped rather than guessed. Phylum-less
taxa: **143 → 81** in MySQL, 40 → 36 in the mirror. The remaining 81 have no
peer carrying a phylum either, so nothing in the database can answer for them.

### The backfill surfaced a cross-kingdom homonym the empty cell had hidden

Filling those 62 phyla made the pre-push check fail with one disagreement:
`Firmicutes -> Arthropoda`. The fill was **correct**; what it did was expose a
collision that an empty cell had kept invisible, because a NULL phylum has
nothing to disagree with.

`Moorella` is two genera. NCBI holds taxid 44260 `Moorella <firmicutes>` as a
**synonym** and 216400 `Moorella <wasps, ants & bees>` as the **scientific
name** — and `normalize_reclassified_phyla.py` matches a row with no tax ID by
name, preferring scientific names. So it resolved our Firmicute to the wasp and
reported an animal phylum. *Bosea* (plant and bacterium), *Bogoriella* and
*Copromorpha* are the same shape, already known here.

Adding the pair to `APPROVED` would have been the wrong fix: that set authorizes
a *reclassification*, so the entry would have permitted moving **every**
Firmicutes row to Arthropoda. Instead the collision is discarded, by a rule that
generalises — when a row matched **by name** has a stored superkingdom that
disagrees with the resolved node's domain, the match is a homonym and is thrown
away. Rows carrying their own tax ID are exempt, since they are not homonym
victims and a genuine cross-domain move would be real.

The counter is reported as `homonym_discarded` rather than left silent, which is
the same lesson as the rest of this sweep: a guard nobody can see is a guard
nobody can check. (The dump also ranks Bacteria and Archaea as `domain` now
rather than `superkingdom`; both spellings are accepted so the guard does not go
inert when the dump is refreshed.)

### Two association-layer "defects" that turned out to be deliberate

Both looked like the vocabulary fragmentation this sweep had already fixed four
times. Neither was, and measuring is the only reason nothing broke.

**`phenotype_comparisons.method`** holds `LEfSe` (468) beside
`as reported (lefse)` (1), plus `as reported (wilcoxon)`, `(t_test)`, `(anova)`
and `(anova_tukey)`. That reads as one method under a wrapper phrase. It is not:
**three of those five sit beside MicrobiomeHD's standardized re-analysis of the
very same study, arms and rank.**

| study | the paper's own result | MicrobiomeHD's re-analysis |
|---|---|---|
| pmid25432777 | `as reported (wilcoxon)` — 13 assoc | `Kruskal-Wallis + BH FDR` — 11 |
| pmid24034618 | `as reported (anova)` — 10 | `Kruskal-Wallis + BH FDR` — 16 |
| pmid25417156 | `as reported (t_test)` — 9 | `Kruskal-Wallis + BH FDR` — 11 |

The prefix is what keeps publication apart from re-analysis, and `method` is part
of `comparison_key` precisely so both can coexist —
`resolve_comparison_method()` in `pipeline.py` already says so. Folding
`as reported (lefse)` onto the 468 GMrepo-sourced `LEfSe` comparisons would also
mix publication-reported findings with repository-sourced ones. **Normalizing
this column is not a cleanup; it is a loss.**

**`taxon_disease_associations.effect_type`** holds three labels for rows with no
`effect_size` — `reported` (359), `q_value` (59), `consensus` (71) — which looks
like three ways of saying the same nothing. But the column is part of
`uq_taxon_disease_evidence`, and `load_microbiomehd.py` uses it as a
**supersession marker**: it deletes rows for the same `(taxon, disease,
comparison)` whose `effect_type` differs, so a finding that arrives as `q_value`
and later gains an effect from the S5 table *flips* to `log2_fold_change` rather
than inserting a second row. A fold would break that, lose the distinction
between "we hold a q-value" and "we hold only a direction", and **be undone on
the next load**.

Neither column had a schema comment, which is why both read as defects an hour
ago. Both now carry one, and `EffectTypeSentinelTests` pins the sentinels so a
future tidy-up fails in the test suite rather than in the data.

### 408 rows disagree with NCBI about their own rank

`taxonomic_rank` in this database is the abundance **layer** a row was loaded
into -- GMrepo publishes a genus array and a species array, and the
per-`(sample, rank)` sum to 1 is defined over that split. It is not a claim about
what the organism is, and for 408 rows the two differ:

| stored | NCBI says | rows | examples |
|---|---|---|---|
| genus | **family** | 90 | *Muribaculaceae* (2005473), *Tannerellaceae* |
| genus | **order** | 19 | *Candidatus Gastranaerophilales*, *Bivalvulida* |
| genus | **class** | 7 | *Candidatus Sericytochromatia*, **Mammalia** |
| genus | **phylum** | 5 | *Candidatus Latescibacterota* |
| species | **strain** | 277 | *Bifidobacterium longum* (1322347) |
| species | subspecies | 4 | *Streptomyces lavendulae* (58340) |

**The 123 coarse genus-layer rows are not an error to move**, and checking that
was the point. GMrepo's own export labels them `taxonomic_rank = genus`, so the
database is faithful. More importantly they are the unclassified **residual**
within a family rather than a duplicate of it: mean shares run 0.2% to 1.6% while
the real genera beside them carry far more, which is exactly what a 16S
classifier produces when it places a read to family but no further. 3,443 samples
hold such a row alongside genuine genera of the same family, the parts are
**disjoint**, and the sum to 1 is correct. Taking them out of the genus layer
would break that invariant and discard real mass -- 16,098 abundance rows and
about 1% of the layer's total.

So nothing moved. What was missing was never the rows; it was that nothing *said*
they were coarser, so a model treated *Muribaculaceae* as a peer of
*Bacteroides*. `taxa.ncbi_rank` records the fact and `v_abundance_genus` exposes
`coarser_than_genus`.

One of them ties back to an earlier finding: **`Bivalvulida` is an order**, and it
was the lone survivor in 82 of the 269 one-genus profiles. A degenerate profile
whose single feature is a myxozoan *order* is now visibly that.

**The 283 species-layer rows are left alone deliberately.** A species-named row
carrying a strain tax ID is the `retaxid` case in `merge_duplicate_taxa.py`'s
`OVERRIDES`, and two concurrent branches are working through taxon identity right
now -- the concatenated-genus repair and the gram-stain corrections both came out
of it. Repointing tax IDs from here would collide with that, so they are recorded
in `data/ncbi_rank_audit.csv` and reported, not touched.

### The same metabolite zero meant absence in 11 cohorts and "average" in 2

`sample_metabolite_levels` mixes scales that differ **in kind**, and nothing
recorded which:

| cohort | levels | negative | zero | max |
|---|---|---|---|---|
| YACHIDA_CRC_2019 | 130,242 | 0 | 78,400 | 665,180 |
| KIM_ADENOMAS_2020 | 89,760 | **46,514** | 513 | 3.3 |
| SINHA_CRC_2016 | 46,374 | **23,553** | 281 | 7.6 |
| KOSTIC_INFANTS_DIABETES_2015 | 13,471 | 0 | 13 | 3.07×10¹⁰ |
| KANG_AUTISM_2017 | 2,596 | 0 | 407 | 461 |

Two cohorts arrive **already log- or z-transformed** -- SINHA holds 23,553
negative levels against 22,540 positive with a mean of −0.04 -- while eleven are
non-negative with maxima spanning **nine orders of magnitude**.

**The hazard is the zero.** On a non-negative scale it means the compound was
not detected: 78,400 of YACHIDA's 130,242 levels, **60% of them**. On a centred
scale a zero sits mid-distribution and means *at this cohort's average* -- the
opposite. Across the table that is **127,836 zeros meaning absence and 794
meaning average**, with nothing to tell them apart.

The standing rule to "normalize within study" does not cover this. A uniform
`log1p` over the table returns **NaN on 70,067 rows** and takes the log of a log
on 136,134 more.

`studies.metabolite_scale` now records it, derived from each cohort's own range
rather than asserted: a negative level cannot be a concentration, so it is
positive evidence of an upstream transform. The absence of negatives is *not*
evidence of a unit, only of non-negativity -- which is exactly what the other
label says, because the source records units nowhere.

### `ml_ready_deduplicated` keeps one run per sample, not one sample per person

The comment on that flag said it "keeps one run per biological sample, which is
what a split should use". The second clause overstated it. A longitudinal study
sampling one subject at seven timepoints produces seven *distinct biological
samples*, and all seven survive the flag:

| phenotype | rows | people | overstatement |
|---|---|---|---|
| Irritable Bowel Syndrome | 305 | **51** | ×6.0 |
| Diabetes Mellitus, Type 1 | 66 | **11** | ×6.0 |
| Fatigue Syndrome, Chronic | 319 | 147 | ×2.2 |
| Health | 1,103 | 613 | ×1.8 |

389 subjects contribute 1,340 of the cohort's rows, so **971 are repeat visits by
someone already in it** -- and a random split puts the same person on both sides.
`ml_ready_one_per_subject` keeps one row per person where the source names the
person, taking the cohort from 13,328 to 12,357.

**Read that flag knowing what it cannot see.** `subject_id` is NULL for all
33,630 GMrepo samples and all 120 from MGnify -- 10,961 of the 13,328
deduplicated rows. A longitudinal *GMrepo* study is invisible to it, and it only
genuinely constrains the 2,367 MullerMM and PubMed rows. It is an improvement on
`ml_ready_deduplicated`, not a solution to pseudo-replication. (The SQLite mirror
has no `subject_id` column at all, so the flag is MySQL-only.)

### A key-alias table, so a question can find all of its data

`sample_attributes` keeps 362 keys exactly as the source wrote them, which is
right — a value can always be traced back to the file and column it came from —
and useless for asking a question. **Antibiotic exposure is spread over seven
yes/no keys**, so `WHERE attribute = 'Antibiotics'` finds 518 of the 1,484 rows
that answer it.

`sample_attribute_aliases` maps each key to a **concept**, and
`sample_attribute_value_aliases` carries the **per-key decoding**.
`v_sample_attributes_canonical` joins both beside the raw columns.

| concept | keys | rows |
|---|---|---|
| `antibiotic_exposure` | 7 | 1,484 |
| `smoking_status` | 3 | 499 |
| `alcohol_amount` | 2 | 455 |
| `immunosuppressant_use` | 2 | 360 |
| `antibiotic_agent` | 5 | 309 |
| `infant_feeding` | 3 | 226 |
| `delivery_mode` | 2 | 173 |

**The value map is the point, not the name map.** `1` means *yes* under
`Antibiotics` and *no* under `Antibiotic_in_the_past_month._1_no_2_yes`. A
name-only alias would let those pool, and the damage is measurable:

| read | "yes" count |
|---|---|
| naive — treat `1`/`yes`/`TRUE` as yes | **602** |
| decoded per key | **309** |
| rows the naive read inverts | **293** |

Same family of error as pooling an LDA score without its comparison order, which
this schema already guards against.

**`concept` is the question, not the topic.** Antibiotics carry three —
`antibiotic_exposure`, `antibiotic_agent`, `antibiotic_start_date` — and one
`antibiotic` concept would put drug names and booleans in the same column.

**Three false clusters were produced by name heuristics and rejected by reading
the values**, and they are in the tests so they stay rejected:

- `steroids` and `mesalamine` are **distinct drugs**, not members of an
  immunosuppressant concept. Only `immunosuppressant` and
  `Immunosuppressive_medication` answer one question.
- **`Breast cancer` is a diagnosis** — a `%breast%` match pulls it in beside the
  infant-feeding keys.
- `GADA_Level`, `IA2A_Level`, `IAA_Level`, `AAB_Level`, `ICA_Level` are **five
  different autoantibodies**, not one "level".

**331 of the 362 keys are unmapped, and that is the normal case** rather than a
backlog: most are a single cohort's own question and need no alias. The view
returns them with `concept` NULL rather than hiding them. `boolean_undecoded`
flags any boolean value the map does not cover — currently **0** — so an
unanticipated answer surfaces instead of being guessed at.

### The flag only covered one direction, and my own fix exposed the other

`coarser_than_genus` flags a genus-layer feature that NCBI ranks as a family,
order, class or phylum — 123 taxa, 16,098 abundance rows. It says nothing about
the opposite case, and **the previous fix created four of those**: the three
bracket organisms NCBI leaves without a genus, plus *Phenoliferia psychrophila*,
carrying 1,233 rows. A modeller saw them as ordinary genera.

Worse, the species layer had the same hole and no flag at all — not even
`ncbi_rank`:

| layer | features | abundance rows | share of layer |
|---|---|---|---|
| genus, coarser | 123 | 16,098 | 0.9% |
| genus, finer | 5 | 1,237 | 0.07% |
| **species, rank differs** | **183** | **189,397** | **14.8%** |

Both views now expose the rank and the flags, with `not_a_genus` as the single
predicate for "features that are genuinely genera".

**The species-layer flag is carefully worded, because the obvious reading is
wrong.** Those 181 are *not* strain-level features: the name is a species and the
stored tax ID points below it — *Bifidobacterium longum* carries 1322347, a
strain, where the species is 216816. That is the `retaxid` class other sessions
are working through, so the flag means "this identifier needs moving up", not
"this is a strain".

A counter keeps the three genus flags consistent — they are three hand-written
rank lists over one fact, so an edit to one that misses another would silently
disagree. It immediately reported **`FAILED_TO_RUN` against the mirror**, which
had no `ncbi_rank` column at all because that backfill was MySQL-only. That is
the one-store gap again, this time in a column rather than a repair, and the
"a check that cannot run is a failure, not a skip" rule is what surfaced it
rather than letting the mirror pass quietly. The column, the backfill and the
view flags now cover both stores: **51 counters, all zero, none failing to run.**

**One thing I did badly:** my first attempt to extend the backfill cut `main()`
in half, leaving a function referencing variables from the one it was split out
of. `git checkout` and a surgical append instead of a clever refactor.

### Unwinding the misattribution: 5.29 of mass taken back out of three genera

The repoint left no audit, so the amount could not be looked up — it had to be
**measured**, and that is what made the reversal safe. For every candidate
sample the script requires MySQL's stored genus value to equal
`source(parent) + source(bracket)` to within 1e-6 before touching it. **All 1,215
reconcile; none was skipped.**

The source files disagree for 682 of them — `gmrepo_reload.csv` holds *exactly
twice* the value in `gmrepo_run_abundances.csv`, a post-duplication artifact — so
rather than pick a file, the reconciliation tries every candidate and accepts the
one that balances against what the database actually stores. **The database
decides, not the CSV.**

| parent genus | samples | what happened |
|---|---|---|
| *Clostridium* | 953 | 1,006 split, **209 whole-row moves** |
| *Bacteroides* | 218 | |
| *Eubacterium* | 44 | |

**209 of those samples had no real *Bacteroides*/*Clostridium*/*Eubacterium* at
all** — the genus row existed *only* because of the repoint. `SRR10377248` is the
clearest: 11.6% of its genus composition was filed under *Bacteroides* and none
of it was *Bacteroides*.

**5.2950** of summed relative abundance moved out of the three genera and onto
three new genus-rank features — `[Bacteroides] pectinophilus` and friends, with
Firmicutes and the NCBI family, matching what the mirror now holds. Mass is
conserved exactly, so sum-to-1 never moved off 0 and that was the end-to-end
check.

**Two of my own bugs surfaced doing it.**

The `taxon_counted_twice_in_a_sample` counter I added last round grouped by
`(sample, tax ID)` **across ranks**. The genus and species layers are separate
compositions, so one organism in both is normal — that is how every organism
works. It read 0 by luck and would have fired the moment a genus-layer feature
took the tax ID its species-layer twin already had, which is precisely what this
fix creates. It now groups per rank.

And a clean re-run **overwrote the recovery record with 0 rows**. This audit
records what was *done*, so it must persist — the opposite of
`phylum_reclassification_audit.csv`, which records what is *outstanding* and is
correctly truncated when a run comes back clean. I drew that distinction after
wiping 1,215 rows, and rebuilt the file from the database: the bracket row holds
the value that moved, the parent holds the remainder, and their sum is what was
stored before. The reconstruction matched the apply exactly — 1,006/209,
218/953/44, 5.2950 — which independently confirms the write did what it planned.

### A repair script whose premise had expired, and would have destroyed 1.27M rows

Extending `fix_genus_layer_purity.py` to the mirror meant reading it first, and
it is **no longer safe to run at all**.

It selects strays with `WHERE t.taxonomic_rank <> 'genus'`, which was exact when
every abundance in the database was a genus-level layer. The **species layer was
loaded afterwards** and now holds 1,277,809 rows — every one of which that
selector matches. Its dry run today reports:

```
abundance rows at a non-genus rank: 1277809
  rows moved to genus : 1277809
```

One `--apply` from collapsing the entire species layer into the genus layer and
adding species mass to genus values. It now **refuses to run** when the species
layer is present, rather than relying on whoever reads the output to notice the
count.

### And the repair it performed is itself questionable

NCBI's square brackets mean the name is **misapplied**, which is exactly the
warning against repointing to the leading token. For all three organisms NCBI
assigns **no genus at all** — they are *incertae sedis*, and all three are
Bacillota:

| tax ID | name | NCBI family | NCBI genus |
|---|---|---|---|
| 84026 | `[Clostridium] methylpentosum` | Oscillospiraceae | **none** |
| 384638 | `[Bacteroides] pectinophilus` | none | **none** |
| 51123 | `[Eubacterium] saphenum` | Anaerovoracaceae | **none** |

So repointing put a **Firmicute's mass under the Bacteroidetes genus
*Bacteroides***. Measured in the mirror, which still holds the un-repointed rows:
**142 rows averaging 1.5% of a sample's genus composition and reaching 11.6%** —
in the most-measured genus in this database. Not a rounding error.

**So I did not replicate it.** `fix_mirror_bracket_names.py` restores the lost
opening bracket — `Clostridium] methylpentosum` → `[Clostridium] methylpentosum`
— and **leaves the abundance exactly where it is**. Mass untouched, so no rank
group changes and no renormalisation is needed; sum-to-1 stayed at 0 violations
through the whole operation. The row remains a genus-layer feature that is now
honestly named, and `taxa.ncbi_rank` already describes it as something NCBI ranks
as a species with no genus — which `v_abundance_genus` surfaces as
`coarser_than_genus`.

Those three stay listed as `spans_two_ranks` in the merge audit, which is
correct: one is genus-layer mass and the other a species-layer taxon, and merging
across ranks is what broke 663 rank groups earlier.

**Outstanding, and yours to call:** MySQL's genus layer still carries that
misattributed mass, and the repoint left no audit — but the mirror holds the
original un-repointed values, so it **is** recoverable by reading across. Undoing
it means taking mass back out of *Bacteroides* in ~142 samples and renormalising,
which is a real change to the most-measured genus here.

### The mirror was double-counting an organism in 582 cells, invisibly

`merge_duplicate_taxa.py` fixed this for MySQL and is MySQL-only, so the mirror
held the pre-merge state: **104 tax IDs over 208 rows**. They are the
species-level form of the Firmicutes/Bacillota split already repaired at four
higher ranks — same tax ID, two names, one organism:

| tax ID | the two names | abundances split |
|---|---|---|
| 33038 | *Ruminococcus gnavus* / *Mediterraneibacter gnavus* | 8,363 |
| 40520 | *Blautia obeum* / *Ruminococcus obeum* | 8,029 |
| 165179 | *Prevotella copri* / *Segatella copri* | 4,688 |

And **582 `(sample, tax ID)` cells counted one organism twice**, across 12,392 of
the mirror's 23,494 samples. The per-rank sums were still exactly 1.0, because
renormalisation happened *after* the duplication — so the sum-to-1 counter could
not see it and **none of validate's 47 counters did**. The duplicated organism
held double share while every other taxon in those samples was squeezed.

The decision logic is **imported** from the MySQL script, not reimplemented: a
second canonicalisation rule would be a second answer to "which name is this
organism", and the stores would disagree about it the way they disagreed about
everything else here. NCBI's current names come from the local taxdump instead of
the API, so it needs no network. 99 merges applied, every one keeping NCBI's
current name — *Klebsiella aerogenes* over *Enterobacter aerogenes*,
*Lancefieldella parvula* over *Atopobium parvulum*.

**I got it wrong twice, and the backup is why that was cheap.**

First run: **663 rank groups stopped summing to 1**, some down to 1e-04. Three
duplicate groups span *two ranks* — the bracket-notation artifacts,
`Eubacterium saphenum` beside `Eubacterium] saphenum` with the opening bracket
lost — and the mangled row sits at genus rank carrying genus-layer abundance.
Repointing it onto a species-rank canonical moved mass from one composition into
another. The abundance layer is defined per rank, so a group spanning two ranks
is never mergeable whatever the tax ID says. MySQL never hit this because
`fix_genus_layer_purity.py` had already repointed its genus-layer binomials; the
script that fixed MySQL never had to refuse them.

Second run: **49 groups still off**, at about 0.99. That one is correct
arithmetic rather than damage — removing a double-counted row *should* drop the
group below 1.0, because the values were normalised while the duplicate was
present. The MySQL script rescales the touched groups afterwards and I had
omitted the step. 557 samples lost a copy; **891 rank groups rescaled**.

Five duplicate groups remain and are listed in the audit: one whose `OVERRIDES`
action names MySQL row ids, one where ***Vibrio fischeri* and *V. furnissii*
share a tax ID** — a mis-assignment, not a duplicate, which the automatic rule
correctly refused — and the three cross-rank artifacts. The new counter therefore
gates on **the harm**, `taxon_counted_twice_in_a_sample`, which is 0 in both
stores, rather than on duplicate tax IDs with five documented exceptions.

### The LDA sign machinery had no invariant either, and two rows needed it

An LDA score's sign is meaningless without its comparison order — this schema's
central design point. Four counters now enforce it, and finding them took one
wrong turn worth recording.

**My first check asked the wrong question.** "Does the arm named by the score's
sign equal the disease?" is `TRUE` for an *enriched* association and `FALSE` for
a *depleted* one, so reading a 0 as failure made **9,183 perfectly correct rows
look like disagreements**. The real invariant is an XOR — the signed arm equals
the disease *iff* the direction is `enriched` — and it holds at **0 violations**
across all 16,399 scored associations. The counter is written that way with the
trap noted, since it is easy to make twice.

**One scored association per store had no sign order at all**: *Bacteroides* /
Vitiligo at LDA 4.31, so the 4.31 could not be read. 56 comparisons in MySQL and
45 in the mirror carry no sign order and that is fine — the rest hold only
unscored associations, where sign order has nothing to disambiguate.

Deriving the arms back from `direction` looks circular, and checking that it
isn't was the real work. `association_direction()` returns `'marker'` when it has
no enriched group to work from — so a row with a score, no sign order, and a
direction of `enriched` **cannot** have had that direction derived from the sign
order. It came from an explicit source field, and for this row
`literature_vitiligo_markers.csv` confirms it: a `direction` column reading
`enriched`, both `enriched_in` columns empty, and a note saying *"Highest LDA
score among vitiligo-group markers"*.

**A second wrong turn:** my guard refused because 3 of the comparison's 8
associations read `depleted`. But LEfSe always reports markers of *both* groups,
so a mixed set is normal — only the *scored* rows constrain the sign order, and
there is one. The guard fired on my own bad premise, which is what a guard is
for; the repair is now a rule rather than a hardcoded row, so it covered the
mirror's copy too under a different `comparison_key`.

### The mirror kept a method name where the number's kind belongs

Found in passing while fixing the above: that row's `effect_type` is `LDA` in
MySQL and **`LEfSe`** in the mirror. `effect_type` records what the stored number
*is*, not how it was produced — and `resolve_effect_type()` exists to convert.
MySQL goes through it; the mirror's build never did, so **123 rows** kept the raw
method.

| | mirror before | after |
|---|---|---|
| `LDA` | 16,072 | **16,073** |
| `reported` | 195 | **317** |
| `LEfSe` | **123** | — |

122 had no score → `reported`; one carried 4.31 → `LDA`. Those are the resolver's
own answers, and the mirror's `LDA` count now matches MySQL's exactly. A counter
rejects any `effect_type` that names a method, so neither store can drift back.

### The view that was once catastrophically wrong had no invariant

`v_taxon_disease_evidence` reported *Faecalibacterium prausnitzii* as **enriched**
in ulcerative colitis — the reverse of one of the field's most replicated
findings — and 32 others with it. It is the most consequential view here, the
only one known to have been badly wrong, and until now **nothing checked it**.

Nine invariants, all reading 0 in both stores:

| counter | what it protects |
|---|---|
| `evidence_direction_on_disease_vs_disease` | the actual historical bug |
| `evidence_study_votes_not_partitioning` | one study, one vote |
| `evidence_votes_exceed_studies` | every study falls in exactly one bucket |
| `evidence_direction_without_a_case_control_vote` | no direction from nothing |
| `evidence_missing_direction_despite_votes` | and none lost |
| `evidence_direction_contradicts_its_counts` | a reversed `CASE` arm |
| `evidence_agreement_ratio_out_of_range` | a majority is in [0.5, 1] |
| `evidence_row_drift` | one row per (taxon, disease) pair |
| `specificity_depleted_fraction_out_of_range` | the same defect's other victim |

**Only the first reproduces that bug, and saying otherwise would overstate
them.** The old view pooled disease-vs-disease with case/control, so all **753**
pure disease-vs-disease pairs carried a direction built from the wrong question —
that counter fires on every one. The second rule it broke, two independent
`COUNT(DISTINCT study_id)` letting one study vote on both sides, **cannot be
caught on today's data**: no `(taxon, disease, study)` group has vs-Health
associations that disagree with each other, so `n_studies_split` is 0 everywhere
and the old form would double-count nothing. That counter is a guard against a
future regression, which is worth keeping on exactly those terms.

**One query, eight counters.** A tuple key returns a single row of counts, so the
CTE-heavy view is materialised once per store while each invariant keeps its
name — a combined total would say a direction was wrong without saying which rule
broke.

**And I had to correct my own justification for that.** I wrote that nine
separate checks "took validate from 32s to 67s", then measured: eight queries
cost **0.89s** against MySQL and the combined one **0.16s**. The collapse saves
0.73s — real, worth having, and *not* why validate takes a minute. That is the
mirror pass, at **36.8s of 65.2s**, mostly the sum-to-1 check over the mirror's
2.05M abundance rows. Validate roughly doubled when the second store started
being checked at all, which is the cost of checking it, not of any one counter.

### Nothing validated the mirror, and it was holding a repaired defect

`validate_database` is typed to MySQL, so the mirror's 8 views and 8 tables had
**no invariants at all** while MySQL's had twelve. That is this project's
signature defect — a guard that reached one store only — showing up in the
*checking* rather than in the data, after the `(SeqCode)` annotations,
`body_site`, the trait vocabularies and `studies.data_type` had each been it in
the data.

`validate_sqlite_mirror()` runs the **same query objects**, not rewritten copies,
so the two stores cannot be checked against different definitions. 17 of the 20
checks are portable; three depend on `sample_attributes`, which the mirror has
no table for. Labels are prefixed `sqlite_`, which means the inverted gate from
the previous round picks them up with **no change on the hook's side** — the
fail-safe default paying for itself one round later.

It found two things immediately.

**`gram_stain_against_phylum` could not run at all**, because it used
`SUBSTRING_INDEX` and `GROUP_CONCAT` to pick a per-phylum majority. Rewritten
with CTEs it runs on both — and the mirror reported **6**. A concurrent session
had corrected those seven gram stains in MySQL; the mirror still carried six of
them, including *Parabacteroides distasonis* with **9,291 abundance rows there**.
The check reports a query that cannot run as a failure rather than skipping it,
because a silently absent counter is how the mirror went unchecked in the first
place.

**The repair script keyed on MySQL `id`.** Extending it to the mirror keyed on
`scientific_name` instead, because the mirror has its own id space: *Vibrio
fischeri* is 1907 in MySQL and **1904** there, so reusing those ids would have
corrected whatever organisms happened to hold them. The phylum and
wrong-value guards still have to pass, and a row with no gram stain recorded is
reported as nothing-to-do rather than counted as a skip, so the exit code stays
reserved for a row that needs a human.

**And the gate needed one fix of mine.** The mirror repeats every label with the
prefix, census figures included, so `sqlite_taxa: 5,308` read as a defect and
would have aborted every push — which I caught by running the gate against real
`validate` output before installing it. It strips the prefix before the census
comparison rather than listing all eight twice, since a new census figure added
in one place and forgotten in the other would block every push.

Both stores now read 0 on all 21 defect counters.

### Seven counters that gated nothing, and a hook list that had drifted twice

Nothing in the test suite reaches a SQL view. The 133 tests are pure functions so
they need no credentials — which is the right trade — but the views carry some of
the most load-bearing logic here: `v_taxon_disease_evidence` reversed 33
directions when it was wrong, `v_ml_ready_samples` overclaimed what its flag
guaranteed, and the two defects in `v_attribute_concept_coverage` above were
caught by reading its output, not by any check.

`validate` now carries seven **view-cardinality** counters, because fan-out is
the failure a view hides best: a `LEFT JOIN` on a non-unique key turns one base
row into several, and every count, mean and fraction built on it is silently
wrong while each individual row still looks right.
`v_sample_attributes_canonical` joins two alias tables and would be first to go.
All seven read 0, and `validate` takes 32s.

**Then the counters turned out to gate nothing.** The pre-push hook checked a
hand-written list of four, and that list had drifted **twice**:
`gram_stain_against_phylum` (added by a concurrent session) and all seven new
ones were printed and ignored, so a defect they caught would have been pushed
anyway. A counter nothing gates on is decorative — the same lesson as the
`homonym_discarded` guard that was silent until it was printed.

So the list is **inverted**: every counter must read zero *except* the census
figures (`taxa`, `samples`, `abundance_rows`, …) named explicitly. That makes the
default fail-safe — a new defect counter is enforced the moment it exists, and
adding a census figure takes a deliberate entry. Getting *that* wrong fails
loudly, since a census figure treated as a defect blocks every push.

Checked against six synthetic cases before installing, including the one that
matters: **a counter nobody added to any list, non-zero, aborts the push.** The
hook also refuses to pass when `validate` prints no defect counters at all,
rather than reading an empty result as clean.

### Reaching the concepts showed they are too sparse to pool

The alias tables make a question reachable. `v_attribute_concept_coverage` says
whether reaching it is worth anything, and for every concept seeded the answer is
**only within a few studies**:

| concept | kinds | ml cohort | studies | studies where it varies |
|---|---|---|---|---|
| `antibiotic_exposure` | boolean | **709** (5.3%) | 5 | **5** |
| `alcohol_amount` | quantity, interval | 455 | 2 | n/a |
| `probiotic_exposure` | boolean | 344 | 1 | 1 |
| `smoking_status` | category | 259 | 3 | 2 |
| `immunosuppressant_use` | boolean | 162 | 2 | **1** |
| `antibiotic_start_date` | date | **0** | 1 | n/a |

**The best-covered concept reaches 5.3% of the modelling cohort**, and all five
of its studies are MullerMM — so **no GMrepo sample carries antibiotic exposure
at all**. A model trained on the pooled cohort cannot adjust for antibiotic use,
however cleanly the concept now decodes. `antibiotic_start_date` reaches 34
samples and none of the cohort.

`n_studies_with_both_values` is what separates sparse-but-usable from useless: a
covariate that is single-valued inside every study it appears in is collinear
with the study itself, the same trap as a disease arm drawn from one cohort. All
five antibiotic studies carry both values, so *within* those five it is a real
covariate — it is the pooling it cannot survive. `immunosuppressant_use` has two
studies and only one that varies, because `Immunosuppressive_medication` records
only "No".

**Two defects in that view were mine, and they are the pattern this sweep keeps
finding.** `MIN(value_kind)` over the ENUM reported `infant_feeding` as
`boolean` when it mixes a category key with two booleans — a column answering a
different question than its name — so it now reports every kind. And
`n_studies_with_both_values` returned **0** for a quantity, which read as "never
varies" when it only meant "no value map applies"; `smoking_exposure_index` is a
Brinkman score that plainly does vary. It is NULL there now.

### `value_numeric` held numbers that were not quantities

`sample_attributes.value_numeric` was documented as *"populated only when `value`
parses as a number ... lets a numeric trait be filtered and aggregated"*. Both
halves were false.

**A number fabricated from a string that is not one** (434 rows). `unique_id`
held `'10007539_0'` with `value_numeric = 100075390` -- the underscore dropped
and the digits run together. Not a parse, and it **collides**: `'1000753_90'`
reduces to the same number.

**A code, not a quantity** (687 rows, 2 keys). `Antibiotic_in_the_past_month._1_no_2_yes`
stores 1 and 2 as labels. Averaging labels is meaningless, and the encoding
**inverts against its neighbour** -- `1` means *no* there, while under
`Antibiotics` `1` means *yes*. Pooling the two silently flips 343 rows: the same
family of error as the LDA sign ambiguity this schema already guards against.

Both cleared. 8,048 genuine quantities kept.

**The third class is the interesting one, and it is not a deletion problem.**
13 keys mix numbers with words, and there the number is *correct* for whichever
rows have one -- what is missing is the denominator:

| key | numeric | rows | `AVG(value_numeric)` |
|---|---|---|---|
| `Stage` | 27 | 150 | **0.0** (real values run I to IVB) |
| `Antibiotics` | 74 | 518 | 0.85 |
| `Race` | 131 | 371 | 0.18 |
| `PhysicalExercise_LastWeek.h.` | 4 | 49 | **0.0** |

The activity keys are the sharpest case. Under
`PhysicalExercise_LastWeek.h.` the **only point answer is `0`**, given by 4 of 49
respondents; the other 45 answered `'1-'`, `'1-3'` or `'3+'`. So the mean is
exactly 0.0 -- not merely incomplete but **anti-correlated with the truth**,
because the rows carrying a number are precisely the people who did nothing.
`Cycling_LastWeek.h.` and `Housework_Childcare_LastWeek.h.` are identical in
shape.

Deleting those zeros would not help -- they are right. This is the same defect as
absence in the abundance matrix, where the number was never wrong and the
denominator was missing, so the fix is the same: **expose it**.
`v_sample_attribute_keys` carries `n_rows`, `n_numeric`, `numeric_fraction`,
`mixes_numbers_and_words` and `looks_like_an_interval` for all 362 keys, so a
mean can be read against what it was taken over.

**Still outstanding, and worth knowing before you use it:** antibiotic exposure
is spread across **8 keys** with three incompatible encodings --
`Antibiotics` (no/yes/0/1), `Antibiotic_in_the_past_month._1_no_2_yes` (1/2,
where 1 is *no*), `antibiotic` (No/Yes), `Antibiotics_Last3months`
(No / "Yes, 2 months ago"), `maternal_antibiotics` (0/1), plus
`antibiotic_type`, `Antibiotic_start_date` and `Antibiotic_used`. A query for
antibiotic exposure finds a fraction of the data, and pooling the coded keys
naively inverts 343 rows. Harmonising 362 source-named keys into a controlled
vocabulary is a design decision rather than a repair, so it is reported here
rather than done.

### `data_type` answered at the wrong granularity, and the stores disagreed on it

Fourth instance, on the column that governs taxonomic resolution. GMrepo records
the platform **per run** -- `AMPLICON` or `Metagenomics` in
`data/gmrepo_project_samples.csv` -- and 10 of its projects contain both, so the
loader collapsed them to `data_type = 'mixed'`:

| project | shotgun | amplicon |
|---|---|---|
| PRJNA763023 | 200 | 1,038 |
| PRJNA888860 | 230 | 264 |
| PRJNA604850 | 134 | 354 |
| PRJNA648321 | 48 | 428 |
| PRJEB42155 | 208 | 211 |
| PRJNA905444 | 88 | 236 |
| PRJNA813705 | 60 | 60 |
| PRJEB21933 | 75 | 40 |
| PRJNA730640 | 40 | 32 |
| PRJEB53891 | 20 | 20 |

**3,786 samples** -- 2,683 amplicon and 1,103 shotgun -- whose platform the
source knew for every single one.

That is not cosmetic. Platform sets the taxonomic floor: amplicon reads one
marker gene, near-identical between species of a genus, so it resolves to genus
and no further. Inside a `'mixed'` study some samples therefore carry
species-rank profiles and some cannot, with nothing recording which. Platform is
also the strongest batch covariate here, and both `loso_audit.py` and
`export_ml_matrix.py` read the **study's** `data_type` -- so for these samples
the covariate they carried was the literal word `'mixed'`.

It showed up in the committed audit: `Diabetes Mellitus, Type 1` was reported as
`platforms = 16S+mixed`, treating a summary as a third platform. It now reads
`16S+mNGS`. The aggregate platform-mixing figures are unchanged (+0.175 across
platforms against +0.060 within), because `16S+mixed` already contained a `+`
and was counted as cross-platform either way -- but the label was wrong and any
finer attribution built on it would have been too.

**And the two stores disagreed about `data_type` itself.** The mirror held four
vocabularies for two platforms:

| store | values |
|---|---|
| MySQL | `16S` (236), `mNGS` (96), `mixed` (12), `meta-analysis` (3) |
| SQLite | `AMPLICON` (191), `Metagenomics` (82), `16S rRNA` (28), `shotgun metagenomics` (9), `meta-analysis` (2) |

So `WHERE data_type = '16S'` matched **236 studies against MySQL and zero
against the mirror**. The mirror was also actively wrong rather than merely
differently spelled: it had no `'mixed'` at all, so 10 projects were labelled
`AMPLICON` when the runs the mirror actually holds from them are shotgun --
PRJNA763023 keeps exactly its 200 `Metagenomics` runs and none of its 1,038
amplicon ones.

The fix removes the class rather than the instance. `samples.sequencing_platform`
records the platform per sample, filled from the run where the source has it and
from the study only where the study names a single platform -- never inheriting
`'mixed'` or `'meta-analysis'`, which are not platforms. `studies.data_type` is
then **derived** from those samples by `derive_study_platform()` in
`gutdb/transform.py`, so `'mixed'` is now a summary of a recorded fact instead of
a value that replaced one, and neither store can drift from the other. Studies
with no samples keep their value, since an association-only study has nothing to
derive from and `meta-analysis` is true of those.

### The taxa trait vocabularies, fixed in the mirror only since they were written

Third instance of the same shape. `cleanup_gutdb.py` has normalized these
columns for as long as it has existed, but it runs **only against the SQLite
mirror**, so the two stores held different vocabularies for the same rows:

| column | MySQL distinct | mirror distinct |
|---|---|---|
| `oxygen_requirement` | 21 | 10 |
| `gram_stain` | 7 | 4 |
| `temperature_range` | 6 | 5 |
| `sporulation` | 32 | 2 |

`oxygen_requirement` is the one that costs something. `anaerobe` (576),
`Anaerobic` (50) and `Obligate anaerobe` (35) were three groups rather than one;
`aerobe` (486), `Aerobic` (161), `obligate aerobe` (8) and `Obligate aerobic`
(1) were four. Oxygen tolerance is among the few genuinely mechanistic traits
for a gut community -- the lumen is anaerobic, and a shift toward aerotolerant
taxa is a recognised dysbiosis signature -- so splitting it across case variants
quietly weakens anything that groups on it.

**`sporulation` held three different kinds of value at once:**

| what it held | taxa | what it actually is |
|---|---|---|
| `0` / `1` | 429 | a yes/no flag as integers |
| `Sporulating` / `Nonsporulating` | 183 | the same flag as words |
| `Sporulating (ascospores); reproduction also by budding` | 216 | prose about spore type |
| `Singles`, `Pairs - Singles`, `Chains - Pairs` | 128 | **cell arrangement** -- another axis |

The last group is the real defect: that text answers how cells sit next to each
other, which `cell_arrangement` already exists for. Those move rather than being
deleted, and only where `cell_arrangement` is empty -- which turns out to be **2
of the 128**. For the other 126 `cell_arrangement` already holds the same
reading, so nothing was ever stored *only* in the wrong column.

Fixing the instance without fixing the cause would invite a fourth, so the
vocabularies now live in `gutdb/transform.py` and `cleanup_gutdb.py` imports
them. `scripts/normalize_trait_vocabularies.py` applies them to MySQL. Seven of
the eight columns now hold **identical value sets** in both stores; the two
remaining differences are membership, not vocabulary (MySQL holds one virus the
mirror lacks, the mirror has `Coenocytic Hyphae` filled for three *Mucor* rows
where MySQL has it NULL).

Two things deliberately left alone. `shape` keeps `Yeast-like cell/globose to
elongate` beside `Yeast-like cell, globose to elongate`: they differ by
punctuation inside free prose, and there is no closed vocabulary for shape to
fold them onto -- inventing one would be a taxonomic judgement rather than a
repair. And any value no rule recognises keeps its cell, because the alternative
is discarding a category nobody anticipated.

### `body_site` was three vocabularies across two stores

| store | source | value | samples |
|---|---|---|---|
| MySQL | GMrepo | `NULL` | 33,630 |
| MySQL | MullerMM | `NULL` | 2,464 |
| MySQL | PubMed | `feces` | 479 |
| MySQL | MGnify | `human gastric mucosa` | 120 |
| SQLite | GMrepo | `gut` | 23,494 |

`cleanup_gutdb.py` has filled the mirror's NULLs with `'gut'` since it was
written; nothing ever did the same to MySQL. So `WHERE body_site = 'feces'`
returned 479 rows against one store and none against the other -- the same defect
shape as the `(SeqCode)` annotations, a repair applied to one store only.

`'gut'` was also the wrong word for the job: **the stomach is part of the gut**,
so it could not separate the gastric biopsies from the stool, which is the only
distinction this column is here for. The vocabulary is now the sampled material
-- `feces` (36,573) and `gastric mucosa` (120) -- in both stores.

**Why it was worth fixing.** MGYS00001188 contributes 120 gastric mucosa samples
and **no control arm of its own**, so Gastritis, Metaplasia and Stomach Neoplasms
are each modelled by comparing stomach biopsies against pooled stool controls:

| disease | CLR genus 5-fold AUC | body site alone |
|---|---|---|
| Stomach Neoplasms | **1.000** | **1.000** |
| Gastritis | **1.000** | **1.000** |
| Metaplasia | **1.000** | **1.000** |

Exactly 1.000, and a single indicator variable matches it. The two communities
share no abundant genus -- *Haematomicrobium*, *Helicobacter*, *Tsukamurella*,
*Gordonia*, *Thermus* in the mucosa against *Bacteroides*, *Phocaeicola*,
*Faecalibacterium*, *Prevotella*, *Bifidobacterium* in stool. These were never
disease models. (`Metaplasia` also mixes 30 mucosa with 10 GMrepo stool samples
inside one label, though the stool ten lack demographics and so never reach the
ML cohort.)

`body_site` is a **per-study constant everywhere** -- no source in this database
reports it per sample -- and for every source but MGnify the value is the
source's *documented scope* rather than a reading. `data/body_site_audit.csv`
records which is which in its `basis` column, and
`scripts/normalize_body_site.py` refuses to proceed if any study ever starts
holding two values. `body_site_confound()` and the exporter now warn when an arm
sits at a site with no control to compare against.

### `metabolomics_method` named a file that holds no method

All 13 metabolomics studies held the literal string `'see mtb.map.tsv'`. That
file exists -- its columns are `Compound, KEGG, HMDB,
High.Confidence.Annotation`. It is a compound-to-identifier map. **No method, no
instrument, no units.**

This column carries a standing rule: metabolite levels are the one measurement
here that cannot be pooled across studies, so check `metabolomics_method` before
combining anything. That rule was **uncheckable** -- every study gave the same
answer and the answer pointed at the wrong file.

The method is not in the Muller collection, which redistributes matrices rather
than protocols. What every cohort's `metadata.tsv` does carry is the protocol's
**DOI**, so that is what the column holds now, as `doi:<doi>`. It identifies the
protocol uniquely and makes the rule mechanical: levels are poolable only where
the value matches. **No two of the 13 cohorts share one**, which is the honest
answer to "can I combine these". `load_muller2022.py` writes the DOI going
forward instead of the placeholder.

### Two columns holding the wrong vocabulary

Both were caught by listing distinct values per column, and both are the same
mistake: a word from one axis written into a column that answers a different
question.

**`samples.qc_status`** held `'1'` (29,509), `'0'` (1,378) and `'curated'`
(2,464). The first two are GMrepo's own QC verdict, copied from its `QCStatus`
field; the third I wrote for the Muller samples, which GMrepo never assessed.
Anything filtering `qc_status = '1'` silently dropped all 2,464 of them, and
anything filtering `<> '0'` silently kept them. The provenance was never missing
either -- `studies.data_quality` already reads `curated` for all 13 Muller
studies -- so the value was redundant as well as wrong. Those samples now carry
NULL, which is what "GMrepo did not assess this" means, and the column's meaning
is written into the schema.

**`taxon_disease_associations.effect_type`** held `'LEfSe'` on 125 rows, beside
`LDA` and `log2_fold_change`. LEfSe is the method; LDA is the statistic it
produces. 124 of the 125 had no effect size at all, which is the same shape as
the 235 rows already labelled `reported`. The literature marker CSVs do declare
`effect_type=LEfSe`, and `cleanup_gutdb.py` already read that column as naming a
method -- the pipeline just copied it verbatim into a column that answers "what
is this number". `resolve_effect_type` now maps a declared method to `LDA` when a
score came with it and `reported` when none did; the method itself is not lost,
since `phenotype_comparisons.method` carries it.

The vocabulary is now five values with a clean invariant: `LDA` and
`log2_fold_change` always have an effect size, `reported`, `consensus` and
`q_value` never do.

### The 1,416 taxa with no evidence are not junk

A sweep for dead rows finds 1,416 taxa carrying no association and no abundance,
which looks like an obvious deletion. It is not one. **Every one comes from
MiMeDB, and 1,058 of them carry measured traits** -- gram stain, oxygen
requirement, sporulation, habitat, pathogenicity. They are a reference set,
loaded so those traits are already present if a study ever reports the organism.
Deleting them would remove **58% of the trait reference** (1,058 of the 1,833
taxa carrying any trait) and **35 of the 63 human_pathogen flags**.

The 358 carrying neither evidence nor traits look more deletable until you read
what they are: 62 are eukaryotes, including the plants and helminths that
`fix_curated_csv_lineage.py` documents as deliberate -- *Carica papaya*,
*Citrullus lanatus*, *Nelumbo nucifera*, and nematodes and flatworms. That script
exists because the obvious test, "this phylum looks non-microbial", is wrong
here: Streptophyta and Nematoda are *correct* for those rows. The remaining 275
are ordinary bacteria that simply have not turned up in any study yet, which is a
statement about study coverage rather than about the rows.

**The real defect was the count, not the data.** `validate` reported a single
`taxa:` figure with no indication that a quarter of it is reference material, so
anyone reading that as the number of organisms this database has *observed*
overstated it by the whole reference set. It now reports all three, so the
distinction is visible without having to know it exists:

```
taxa: 5,835              <- snapshot, 2026-10-14
taxa_observed: 4,481
taxa_reference_only: 1,354
```

Treat those as a snapshot, not a fact about the schema. The database is written
to by more than one session: the figure was 5,897 / 4,481 / 1,416 earlier the
same day, and dropped by 62 when a concurrent branch merged 62 mangled
double-genus names into the correctly named rows that already held their
evidence. Run `validate` rather than quoting a number from here.

### Backfilling the lineage that was never loaded

```bash
python scripts/backfill_lineage.py --apply
```

1,417 taxa carried no phylum and 1,563 no superkingdom, and **1,415 of those were
not dead rows** -- they had associations or abundance behind them. *Waltera* alone
accounts for 7,944 abundance rows. Most are recently described gut organisms that
arrived from GMrepo with a name and a tax ID but no lineage, so anything grouped
by phylum silently dropped them.

1,326 had a tax ID, which makes them resolvable offline from the local dump.

| | before | after |
| --- | --- | --- |
| taxa with no phylum | 1,417 | **143** |
| of those, carrying evidence | 1,415 | 141 |
| taxa with no superkingdom | 1,563 | 283 |
| genus-rank abundance mass carrying a phylum | — | **99.54%** |

**Only empty fields are filled, and that is the safety property rather than a
nicety.** Several lineages here were repaired by hand after cross-kingdom
homonyms put *Bosea* under Streptophyta and *Bogoriella* under Ascomycota; a
backfill that preferred the dump over the stored value would quietly undo that
work every time it ran. `genus` and `species` are excluded entirely -- they are
keys, so rewriting them would change identity rather than fill a gap.

The backfill surfaced one thing worth keeping: NCBI appends the nomenclatural
code to names published outside the ICNP, so the dump carries
`Methylomirabilota (SeqCode)` and nine siblings across class, order and family.
That is metadata about how a name was published, not part of the name, and it
creates a rank value of one that would split from the bare name later.
`strip_nomenclatural_code` removes it, limited to `(SeqCode)` and `(ICNP)` in
trailing position -- a parenthetical is not always an annotation, and
`Neobittarella (ex bilen et al. 2018)` is a scientific name this database
legitimately held.

Two survivors of the backfill look wrong and are not. *Bivalvulida* sits under
phylum **Cnidaria** because myxozoans are parasitic cnidarians and GMrepo
genuinely reports it in faecal samples; its `superkingdom` correctly reads
Eukaryota. 143 taxa still have no phylum because NCBI has none for them either.

### Phylum names: renames against reclassifications

```bash
python scripts/fold_phylum_synonyms.py --apply        # run this FIRST
python scripts/normalize_reclassified_phyla.py --apply
```

The two look alike in an audit and are not the same thing, and **the NCBI tax ID
tells them apart**. One node under two names is a *rename*, which folds without
anyone's approval. Two nodes is a *move between phyla*, which is a taxonomic
judgement and needs a line in `normalize_reclassified_phyla.py`'s `APPROVED` set,
because approving a pair moves every taxon that matches it.

Seven phyla were stored under both their classic and their 2021 ICNP name at
once -- 1,365 rows saying Firmicutes beside 97 saying Bacillota, 424
Bacteroidetes beside 12 Bacteroidota, 3,636 taxa in all. Nothing errored;
`GROUP BY phylum` simply reported one phylum as two, which quietly breaks the
Firmicutes/Bacteroidetes ratio. 119 of the affected rows came from
`load_muller2022.py` reading NCBI's lineage verbatim: NCBI's *scientific* name is
now the ICNP one (taxid 1239 is Bacillota, Firmicutes is the synonym), while this
schema deliberately keeps the classic name so a phylum groups as one value.

Run the fold **before** the reclassification script. That script normalizes only
NCBI's side of the comparison, so a stored post-2021 name or typo appears as a
disagreement it cannot act on; folding first leaves it looking at genuine moves,
the only thing it should judge.

**Renames** (same tax ID), folded without approval:

| | evidence | taxa |
| --- | --- | --- |
| Bacillota → Firmicutes, and 5 more of that family | one node, two names | 132 |
| Methanobacteriota → Euryarchaeota | taxid 28890 | 98 |
| Thermoproteota → Crenarchaeota | taxid 28889 | 15 |
| Pseudomonadata → Proteobacteria | a one-character slip | 1 |
| Candidatus thermoplasmatota → Thermoplasmatota | a prefix, not a phylum | 1 |

**Reclassifications** (different tax IDs), approved one at a time:

| organisms | from → to | taxa |
| --- | --- | --- |
| Methanomassiliicoccus, Thermoplasma | Euryarchaeota 28890 → Thermoplasmatota 2283796 | 7 |
| Gonapodya, Monoblepharis | Chytridiomycota 4761 → Monoblepharomycota 3705299 | 4 |
| *Rhizophagus irregularis* | Mucoromycota 1913637 → Glomeromycota 214504 | 2 |
| *Sphaerobacter thermophilus* | Chloroflexi 32061 → Thermomicrobiota 3027942 | 2 |
| *Rhodothermus marinus* | Bacteroidetes 976 → Rhodothermota 1853220 | 2 |

Three further disagreements were neither: **cross-kingdom homonyms** surviving in
the SQLite mirror, which is built from CSVs and predates the repair
`fix_curated_csv_lineage.py` made in MySQL. *Bosea* sat under Streptophyta
because a plant genus shares the name, *Bogoriella* under Ascomycota because a
lichen does, and *Bacteroides pectinophilus* under Bacteroidetes -- NCBI brackets
that name because the organism is actually a Lachnospiraceae. MySQL had all three
right, so the mirror was synced from it rather than given `APPROVED` entries:
approving "Streptophyta → Proteobacteria" as a rule would move any plant that
came after.

Distinct phylum values **86 → 76**, phyla split under two names **7 → 0**,
unapproved disagreements **26 → 1**. The one remaining is not a phylum problem: a
MiMeDB row named `Sphingomonas sphingobacterium paucimobilis`, two organisms
concatenated, whose tax ID resolves to *Sphingobacterium paucimobilis* while its
stored lineage says Proteobacteria. It carries no associations and no abundances.

### Why the species layer stops where it does

Species coverage reads as 34% of all samples, which sounds like a gap worth
closing. Measured against a denominator that can physically carry a species
call it is **90.4%** (12,566 of 13,907 profiled shotgun samples). The other
19,899 samples are 16S, where one marker gene cannot separate species in a
genus, so no amount of loading will move them. GMrepo's shotgun samples are at
**100%**.

The entire shortfall is 1,341 Muller shotgun samples, and the collection's
`species.tsv` tables cannot close it honestly. `scripts/fetch_muller2022.py
--species` will fetch them (221MB for the four cohorts that ship one;
YACHIDA_CRC_2019 is shotgun but ships none, so 294 of the 1,341 have no source
at all), and the measurements that follow are the reason they are not loaded:

| species label | distinct | abundance mass |
| --- | --- | --- |
| maps to a species already in `taxa` | 1,969 | **50.2%** |
| a real binomial this database lacks | 12,220 | 8.2% |
| GTDB accession-derived, no published name (`4484-113 sp002050035`) | 41,966 | **41.6%** |

Loading would therefore drop and renormalize away **49.8% of every sample's
mass**, against 8.2% for the genus layer. Per sample the median retained share
is 50.0% and the 95th percentile is 71.3% -- below the *mean* of the species
profiles already here, which retain about 79.7% by GMrepo's own
`unclassified_fraction`. Renormalization hides all of this, because a profile
built on half the mass still sums to 1 and looks like every other row in
`sample_taxon_abundances`.

The usual escape hatch does not hold either. A mapping loss that falls equally
on both arms cancels in a within-study ratio, but it does not fall equally:
cases retain 4.5 to 7.8 points more mass than controls in three of the four
cohorts and less in the fourth. Renormalizing then multiplies the two arms by
different factors, a systematic directional bias of roughly 15% on every ratio
in a cohort -- the same order as the marginal signals that
`scripts/replication_audit.py` shows the disagreements live in.

So the species gap is, for this source, already closed. Moving the headline from
34% to 37% would mean adding 1,047 samples whose species profiles are built on
half the evidence of the 12,566 beside them, with nothing in the schema to say
so. Closing it properly means re-profiling raw reads against an NCBI-based
classifier, which is a different project from loading a curated table.

### Sparse per-sample metadata

Cohorts report traits that matter for confound control but that almost nothing
else reports: smoking, alcohol, blood pressure, comorbidities, surgery type. In
this collection alone that is **43,406 values over 362 distinct column names**,
and ERAWIJANTARI contributes 41 clinical columns by itself. Widening `samples`
was the wrong shape -- smoking would be populated for 1,011 of 36,693 rows, under
3%, so the honest version of that plan is 362 columns at 95-99% NULL.

`sample_attributes` holds them sparsely instead, so **absence is the absence of a
row** exactly as in `sample_taxon_abundances`. There are no empty cells by
construction, and the table is ~43k rows against that one's 3.07M. `value` keeps
the source's string; `value_numeric` is populated for the 9,169 that parse as
numbers, so a numeric trait can be aggregated without casting.

**Values are not comparable across studies, and often not commensurable.**
Smoking arrives four ways:

| attribute | cohort | n | range |
| --- | --- | --- | --- |
| `Brinkman Index` | YACHIDA_CRC_2019 | 294 | `0` .. `980` (a pack-year product) |
| `Chem ID / Smoking history` | KIM_ADENOMAS_2020 | 240 | `Missing` .. `Smoke` |
| `Tobacco_Use` | POYET_BIO_ML_2019 | 163 | `1-2 per week` .. `Yes` |
| `SmokingStatus` | ERAWIJANTARI_GC_2020 | 96 | `NotSmoking` .. `Unknown` |

Nothing is harmonized on the way in, because harmonizing means discarding the
original scale, and confound control belongs inside a study anyway -- which is
where `adjudicate_conflicts.py` already does its comparisons. Filter by study
before trusting an attribute name to mean one thing. Column names are kept
verbatim including R `make.names` artifacts (`Weight..kg.`, and `Gout...22`
beside `Gout...46`, two different gout columns the source named identically), so
a value traces back to the file it came from. 21 values (0.05%) are in-band
missing markers -- `Unknown`, `?`, `Missing` -- kept because an explicit unknown
is a recorded answer rather than a question never asked, so a consumer should
filter them.

Two things this bought immediately. The paper reports Brinkman's index differing
systematically by group, and that reproduces straight out of the table -- HS
lowest, as described. More usefully, **the MeSH labels collapse YACHIDA's
`Stage_0`/`Stage_I_II`/`Stage_III_IV` into one `Colorectal Neoplasms` arm, but
`Stage` survived as an attribute**, so a stage-resolved analysis is still
available without reloading:

```sql
SELECT stg.value AS stage, COUNT(*), ROUND(AVG(brk.value_numeric), 1)
FROM sample_attributes stg
JOIN sample_attributes brk ON brk.sample_id = stg.sample_id
                          AND brk.attribute = 'Brinkman Index'
WHERE stg.attribute = 'Stage'
GROUP BY stg.value;
```

One genus is refused outright: GTDB uses **`Copromorpha`** for a bacterium, and
the only NCBI taxid of that name (1181387) is a moth, so the domain guard in
`ensure_genera` rejects it. It is dropped *before* renormalizing -- filtering
afterwards would leave its 0.05% of mass in the denominator with no row to carry
it, and every affected sample would fail the per-rank sum-to-1 invariant.

**Metabolite levels are not comparable across studies** -- instruments, targeted
versus untargeted designs and units all differ, and nothing sums to 1, so
`studies.metabolomics_method` records the platform and pooling raw levels across
cohorts is invalid in a way that pooling relative abundances is not. Across the
13 cohorts whose `mtb.tsv` is fetched (iHMP ships a 57.6MB zip), 3,048 of 17,066
features carry a valid identifier, giving 1,206 distinct HMDB compounds of which
**334 appear in 3 or more cohorts** -- which independently reproduces the paper's
own 314-metabolite meta-analysis set. The remaining features are unannotated m/z
peaks with no cross-study identity. Upstream identifier columns are lightly
contaminated (MetaboAnalyst `METPA*` ids and KEGG `C#####` values in the HMDB
column, KEGG DRUG `D#####` ids in the KEGG column), so a loader must validate
`^HMDB\d{5,7}$` and `^C\d{5}$` rather than trust them.

### MicrobiomeHD standardized re-analysis

[Duvallet et al. 2017](https://doi.org/10.1038/s41467-017-01973-8) (PMID 29209090)
re-processed 28 case-control 16S studies through one pipeline and published the
resulting genus-level q-values. Because every dataset went through identical
processing, it is an independent check on the per-project LEfSe results this
database takes from GMrepo, rather than more of the same evidence.

```bash
python populate_database.py validate            # before
python scripts/load_microbiomehd.py             # dry run; prints the curl commands if files are missing
python scripts/load_microbiomehd.py --apply
```

Two schema additions support it: `taxon_disease_associations.p_value` and
`.q_value`. These are independent of `effect_size`, not derived from it — GMrepo
reports an LDA score and no significance, MicrobiomeHD reports an FDR-corrected
q value and no effect size — so a row may carry either, both, or neither, and
NULL means "the source did not report it" rather than "not significant". Loaded
rows use `effect_type = 'q_value'` so they never pool with LDA scores in an
`ORDER BY ABS(effect_size)` query.

Effect sizes come from file-S5, which holds `log2(mean_cases / mean_controls)`.
Three of its values are placeholders rather than measurements — the table
maximum stands in for a fold-change against a zero control mean, the minimum for
a zero case mean, and `0.0` for both means being zero — so storing the maximum
would assert a roughly 1300-fold enrichment the data cannot support. Those 59
rows keep `effect_type = 'q_value'` with a NULL `effect_size`; the 337 real
measurements become `effect_type = 'log2_fold_change'`. The sentinels are derived
from the file at load time rather than hardcoded. S5's signs agree with S1's on
all 375 shared significant rows, which is what confirms the two files are
column-aligned.

Each finding keeps one row. Because `effect_type` is part of
`uq_taxon_disease_evidence`, a row whose type changes between runs would
otherwise insert a second row beside the first, so the loader clears any other
`effect_type` from this source for the same triple before inserting.

Only the 434 rows significant at |q| < 0.05 are loaded, of 2,769 present in the
matrix. The other 2,335 are real "tested, not significant" results, but
`v_taxon_specificity` counts `COUNT(DISTINCT disease_id)` over every association
row with no direction filter, so admitting them as `direction = 'no_difference'`
would push a merely-tested genus toward `broad` or `pan_disease`. They belong in
a separate table if they are ever wanted.

38 further rows are dropped because the RDP classifier labels they carry are not
organisms (`Clostridium_IV`, `Lachnospiracea_incertae_sedis` and similar).

file-S4 is the authors' hand-curated record of what each **original paper**
reported, which makes a third view available beside the re-analysis and this
database's GMrepo evidence. It loads as `effect_type = 'reported'` on its own
comparison per study, whose `method` names the publication's own test
(`as reported (wilcoxon)`, `as reported (lefse)`, …), so it never pools with the
standardized re-analysis on the same study. That is what makes the two
comparable:

```sql
SELECT t.scientific_name, d.name,
       MAX(CASE WHEN a.effect_type <> 'reported' THEN a.direction END) AS reanalysis,
       MAX(CASE WHEN a.effect_type =  'reported' THEN a.direction END) AS as_published
FROM taxon_disease_associations a
JOIN taxa t ON t.id = a.taxon_id
JOIN diseases d ON d.id = a.disease_id
JOIN phenotype_comparisons c ON c.id = a.comparison_id
JOIN studies s ON s.id = c.study_id
WHERE a.source_database = 'MicrobiomeHD'
GROUP BY t.id, d.id, s.id
HAVING reanalysis IS NOT NULL AND as_published IS NOT NULL;
```

Species-level rows load as well as genus-level ones, since `taxa` holds
(genus, species). That matters more than the count suggests: it is what brings in
*Fusobacterium nucleatum* in colorectal cancer at q = 1.3e-05, along with
*Porphyromonas asaccharolytica*, *Peptostreptococcus stomatis* and
*Leptotrichia hofstadii* — findings with no genus-level equivalent in the file.
Subspecies collapse to their species, keeping the stronger q, because
`s__nucleatum;sb__polymorphum` and `s__nucleatum;sb__nucleatum` are one organism
under a (genus, species) key.

**S4 is notes, not a matrix, and only 38 of its 1,027 rows load — under 4%.**
Treat it as a sample of the literature, not a summary of it. The loader reports
the attrition every run: 509 rows sit at a rank above genus or are OTU-level
(`taxa.genus` is NOT NULL, so a family- or phylum-level finding has no key at
all), 210 belong to studies with no counterpart in the re-analysis, 167 name a
family with an empty `g__`, 70 report no parseable q value, 22 are above
threshold, 6 are species rows with an empty `s__`, 3 are subspecies collapsed
into their species, and 1 is a disease-versus-disease comparison.

Threshold notation (`<0.05`, `<0.001`) was measured and deliberately not parsed:
every such row also fails another rule, so interpreting it recovers nothing.

file-S2 is the meta-analysis's own conclusion: a genus is listed for a disease
when it is significant at q < 0.05 in the same direction in at least two of that
disease's datasets, for the five diseases with at least three datasets. It is the
complement of S3 — consistent *within* one disease rather than *across* several —
and because it names a (genus, disease) pair it is an association rather than a
taxon annotation. It belongs to no single cohort, so it is attributed to the
paper itself, which enters as study `PMID29209090` with `data_type =
'meta-analysis'`, under `effect_type = 'consensus'` with no effect size or q
value, since a consensus has neither. 71 calls load: 58 for *C. difficile*
infection, 5 each for obesity and colorectal cancer, 3 for IBD.

### Rejecting classifier labels

RDP and SILVA emit labels that are not organisms, and every file here carries
some: `Clostridium_IV`, `Lachnospiracea_incertae_sedis`,
`Escherichia/Shigella` (two genera the classifier could not separate),
`Ruminococcus2` (a reference-database cluster), `Clostridium_sensu_stricto`.
One rule governs all of them, and it distinguishes creating from referencing:

- A genus the database already holds is always accepted. Nothing is invented,
  and the existing row keeps whatever identification it was loaded with.
  *Raoultella* is why this matters — NCBI synonymised it into *Klebsiella*, so it
  is no longer a genus-rank name in the dump, but it is a real organism already
  present here with a tax ID, and refusing it would drop a real finding.
- A genus not yet present must be a genus-rank name under Bacteria or Archaea in
  the NCBI dump, and must not match the classifier-label pattern.

Rejections are printed on every run rather than counted silently.

Study ids in S4 are abbreviated differently from the re-analysis
(`ra_littman`, `ibd_hut`, `mhe_zhang`) and the repository documents no mapping,
so only exact and unambiguous-prefix matches are used. The rest are skipped
rather than matched on a shared author name — attaching a publication's claim to
the wrong study is worse than omitting it. Genus names are checked against the
NCBI dump before becoming taxa, which is how the source's `Peptosreptococcus`
typo is caught; it is reported and skipped, never silently corrected.

## Data notes

- MiMeDB taxonomic and metabolic text is normalized before matching. Species are matched case-insensitively by `(genus, species)`.
- The derivation rules follow the supplied workflow: facultative anaerobes and conflicting pathways become `mixed`; anaerobes default to `fermenter`; aerobes/microaerophiles default to `respirator`; high-specificity food-source terms take precedence.
- The included IBS marker CSV was transcribed from the supplied PDF. Confirm spellings and values against the live GMrepo project export before treating it as publication-ready evidence.
- GMrepo currently documents downloadable project/run data and REST endpoints for phenotype and abundance access. Marker imports are intentionally file-based so the pipeline does not depend on an undocumented web endpoint.
