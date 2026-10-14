# Gut microbiome database population pipeline

A curated and normalized gut microbiome dataset integrating microbial taxonomy, metabolic traits, oxygen requirements, succinate production, and disease associations, with a focus on IBS-related microbiome analysis.

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
