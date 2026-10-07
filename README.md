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

- **Within-study, never pooled.** Pooling cases and controls across cohorts
  agrees with the curated associations only 62% of the time, which is what
  comparing across different protocols and sequencing runs produces.
- **Absence is zero, not missing.** An undetected taxon has no row in
  `sample_taxon_abundances`. Averaging only the rows that exist compares the
  samples where a taxon was abundant against the samples where it was abundant.
  The denominator is every sample in the arm profiled at that taxon's rank.

**Only case/control contrasts vote.** 2,340 of the 16,889 associations compare
one disease against another rather than against Health, and 745 of the
replicated pairs mix the two kinds. `v_taxon_disease_evidence` votes across both
as though they answered the same question, and it counts enriched and depleted
studies with two separate `COUNT(DISTINCT study_id)` expressions, so a study
holding both contrasts votes on both sides at once. *Faecalibacterium
prausnitzii* in ulcerative colitis shows what that costs: four studies report it
enriched in UC *versus Crohn disease*, two report it depleted *versus Health*,
and the pooled view calls it enriched in UC -- reversing one of the most
replicated findings in the field. The adjudicator rebuilds the association vote
from vs-Health comparisons only, one vote per study, and stores the old pooled
numbers beside it in `data/adjudication_audit.csv` for comparison.
`assoc_contrast_scope` records which contrasts a pair actually had.

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
```

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
- **`Study.Group` needs a hand-written MeSH mapping per cohort**, and several are
  undecodable from the data alone: `0`/`1` in SINHA_CRC_2016, `D`/`H`/`C` in
  MARS_IBS_2020, `MP`/`HS` in YACHIDA_CRC_2019. `Adenoma` is absent from
  `diseases` and is needed for KIM_ADENOMAS_2020.

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
