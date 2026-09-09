# Gut-microbiome–disease database: queryable SQLite build

This repo's pipeline targets MySQL 8 (`gutdb/config.py`, `gutdb/db.py`). This build
is a server-free alternative: it loads the same source CSVs into the same
normalized schema as a single-file SQLite database, so the data can be queried
with nothing running.

It reads the repo's CSVs and `gutdb/transform.py`; it writes only new files and
changes nothing in the MySQL pipeline.

## Files

| File | What it is |
| --- | --- |
| `gutdb/schema_sqlite.sql` | `gutdb/schema.sql` translated to SQLite. Table names, column names, keys and view semantics preserved 1:1; only engine syntax changes (`AUTO_INCREMENT`→`AUTOINCREMENT`, `ENUM`→`TEXT`+`CHECK`, boolean expressions→`CASE`, integer division→`*1.0`) — plus `p_value`/`q_value` added and five always-empty columns omitted, both noted in the file header and below. |
| `scripts/build_gutdb_sqlite.py` | Loader. Imports `gutdb.transform` rather than re-implementing normalization, so taxon keys, phylum folding, effect-size quantization and direction derivation match what the MySQL loads produce. Idempotent — rebuilds the file from scratch. |
| `scripts/cleanup_gutdb.py` | Back-fill and normalization pass; see the section at the end. Idempotent. |
| `scripts/query_gutdb.py` | The seven queries below; writes one CSV per query. |
| `docs/gutdb_overview.png` | Two-panel summary of what landed. |
| `docs/gutdb_cleanup_coverage.png` | What the cleanup pass filled. |
| `docs/redundant_columns.csv` | Column-redundancy audit: 22 verdicts with the evidence behind each and a recommendation. |
| `docs/column_profile.csv` | All 92 columns — percent empty, distinct values, dominant value and its share. |

Outputs (`gut_microbiome.sqlite`, `q1…q7_*.csv`, `coverage_before_after.csv`,
`cleanup_log.txt`) are written to the working directory and are gitignored —
they are regenerable from the three scripts.

Rebuild, clean and re-query, from the repo root:

```bash
python scripts/build_gutdb_sqlite.py   # ~35 s
python scripts/cleanup_gutdb.py        # ~11 s
python scripts/query_gutdb.py          # ~30 s
```

The scripts locate the repo from their own path; set `GUTDB_REPO` to override
that, `GUTDB_SQLITE` to write the database elsewhere. The cleanup pass needs the
NCBI taxonomy dump unzipped into `./taxdump` (`names.dmp`, `nodes.dmp`) — from
`https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip`, or point `GUTDB_TAXDUMP`
at an existing copy. Without it the pass still runs and simply fills less.

## What loaded

| Table | Rows |
| --- | --- |
| `taxa` | 4,374 |
| `diseases` | 161 |
| `studies` | 312 |
| `phenotype_comparisons` | 378 |
| `taxon_disease_associations` | 16,390 |
| `samples` | 33,630 |
| `sample_taxon_abundances` | 1,243,990 |

Sources: the curated taxon table (`gut_microbiome_new_full.csv`, 2,112 rows), the
full GMrepo association export (16,073 rows), the supplementary IBS-project
export, the transcribed `PRJNA705217` marker file, 14 per-disease literature
marker files (318 association rows in total), the project-samples metadata file,
and the species-level abundance matrix. Every source gets a row in
`ingestion_runs` with read/inserted/skipped counts.

## Decisions worth knowing about

**Comparison identity excludes `method`.** A `phenotype_comparisons` row is keyed
on study + the two phenotype groups + effect type. `method` is deliberately not
part of the key: `gmrepo_PRJNA705217_ibs_markers.csv` omits the column while the
full GMrepo export records `LEfSe` for the same rows, so keying on it split one
comparison in two and double-counted all 34 associations. With the key as built,
`PRJNA705217` has one comparison holding 55 associations, and 32 of the marker
file's 34 rows are correctly recognized as already present.

**Two significance columns added.** `taxon_disease_associations.p_value` and
`.q_value` are new. The literature marker CSVs carry these values and the MySQL
schema has nowhere to put them, so they were being discarded. This is the only
schema departure.

**The `Spore` column carries three different things.** In
`gut_microbiome_new_full.csv` it holds a 0/1 sporulation flag, free-text
sporulation phrases, *and* cell-arrangement values (`Singles`, `Pairs - Chains`,
…). Each is routed to the column that describes it — the arrangement values go
to `cell_arrangement`, not into `sporulation`. Unrecognized values (`Unknown`)
stay NULL.

**`energy_mode` needed a synonym map.** The CSV uses `respirer` (159 rows); the
schema's enum says `respirator`. Mapped. Values outside the enum fall back to
`N/A`.

**Genus lineage backfilled from species.** Genus-rank taxa are created by the
association exports, which ship no lineage columns, so `phylum` came out NULL for
all 789 of them — and `v_taxon_specificity` groups by phylum. Each genus now
inherits the modal lineage of the curated species indexed under it (640 fields
filled across 4 ranks; 482 of 789 genera now have a phylum). Genera with no
curated species stay NULL rather than being guessed at.

**One taxon collapse in the source.** `PRJNA993675` reports
*Lacrimispora amygdalina* under two NCBI tax IDs (253257 and 759821) with
identical LDA scores — a synonym pair. Because `taxa` is keyed on
(genus, species), both fold into one row and the duplicate evidence rows are
dropped; `ncbi_tax_id` keeps 253257.

## Integrity checks (all clean)

`PRAGMA foreign_key_check` returns nothing; no orphaned association, sample or
abundance rows; no `relative_abundance` outside [0, 1]; no `direction = 'marker'`
rows left unresolved; one signed effect lacks enrichment context; 5 samples have
no disease label. All 12,049 samples with abundance data sum to exactly 1.0
within the species rank, so `v_abundance_coverage` flags none of them.

## Queries

1. **`q1_content_summary`** — associations by source, effect type and method.
   16,072 GMrepo LEfSe rows (1,297 taxa, 142 diseases, 273 studies) plus 318
   literature rows (174 taxa, 14 diseases). GMrepo skews depleted
   (8,993 vs 7,079); mean |LDA| 3.22 against 4.31 for the literature set, which
   is what you would expect if published markers are the ones that survived a
   significance filter.
2. **`q2_disease_rollup`** — 156 diseases ranked by association count, with study
   and sample counts. Crohn disease (2,195), ulcerative colitis (1,555) and
   colorectal neoplasms (1,447) hold 32% of the database.
3. **`q3_pan_disease_taxa`** — the 587 taxa reported for ≥5 diseases. Across all
   taxa, `v_taxon_specificity` calls 377 single-disease, 418 narrow, 393 broad and
   143 pan-disease. The pan-disease set are general dysbiosis indicators, so they
   will not discriminate between diseases in a classifier.
4. **`q4_replicated_evidence`** — 2,131 (taxon, disease) pairs seen in ≥2 studies.
   59% are unanimous in direction; 875 have at least one study disagreeing, and
   414 sit at agreement ≤0.6, i.e. a near-tie. Sorted conflicts-first, so the top
   of this file is the list of claims that need adjudication before any of them is
   used as a feature.
5. **`q5_ibs_markers`** — the IBS panel: 255 associations, 218 taxa, 4 projects,
   141 depleted vs 114 enriched. Strongest effects are depletions of *Alistipes*
   (LDA −4.79), *Bacteroides* (−4.67) and *Ruminococcus* (−4.31) — but
   *Bacteroides* is +4.66 in `PRJNA637763`, an example of the conflicts in q4.
6. **`q6_cohort_readiness`** — per-study covariate and feature completeness from
   `v_ml_cohort_summary`, 273 studies.
7. **`q7_association_vs_abundance`** — the cross-check that joins the two halves of
   the database: for each species-level pair with a consensus direction and ≥20
   case and ≥20 health samples, the mean relative abundance in cases against
   health. Only 62% of 2,670 pairs move in the direction the curated association
   claims (68% of enrichments, 57% of depletions). The comparison pools samples
   across studies and ignores batch effects, so it is a smell test rather than a
   validation — but a coin-flip-ish depletion rate is a reason to prefer
   within-study effect estimates.

## Not done

- No MiMeDB enrichment. `mimedb_microbes_v2.csv` is in the repo and would fill
  gram stain, oxygen tolerance and the pathogen flag for many taxa; the flag is
  left NULL rather than defaulted to false, as `gutdb/transform.py` intends.
- Abundances are species-rank only in the source file, so `v_abundance_genus` is
  empty. Genus-level abundance would have to be aggregated from species or
  re-fetched from GMrepo.
- No statistical testing on q7 — differences are means, not tested contrasts.

## Cleanup and back-fill pass (`cleanup_gutdb.py`)

Run after the build: `python scripts/build_gutdb_sqlite.py && python scripts/cleanup_gutdb.py`.
Idempotent — a second run reports zero fills, zero normalizations, zero merges.
It only ever writes into empty cells or rewrites a value onto its canonical form;
no association, sample or abundance row is invented, and row counts are unchanged
apart from the three merged disease rows.

**70,648 cells filled across 25 columns**, plus 3 duplicate disease rows merged.
Per-column detail is in `coverage_before_after.csv`; the full log is `cleanup_log.txt`.

### 1. Duplicate disease concepts merged

Three names were the same MeSH concept spelled two ways, so each disease's evidence
was split across two rows:

| dropped | kept | MeSH | rows repointed |
|---|---|---|---|
| `Colitis Ulcerative` | `Colitis, Ulcerative` | D003093 | 147 samples |
| `Crohn's disease` | `Crohn Disease` | D003424 | 235 samples |
| `Hereditary Nonpolyposis Colorectal Cancer` | `Colorectal Neoplasms, Hereditary Nonpolyposis` | D003123 | 18 samples |

The duplicates carried no association rows but did carry 400 samples between them,
so the abundance side of these three diseases was detached from the association
side — `Colorectal Neoplasms, Hereditary Nonpolyposis` had associations and zero
samples while its duplicate held 18 samples. Merging is why q7 now reports 2,678
association/abundance pairs instead of 2,670. 158 diseases remain (was 161).

### 2. Taxon traits and lineage backfilled

Sources, in priority order — a local source is never overridden by a lower one:

1. `gutdb/final_taxa.tsv` (curated: gram stain, shape, oxygen, pH, spore, mobility,
   flagella, membranes, metabolism, energy source, habitat, temperature, lineage)
2. `mimedb_microbes_v2.csv` (MiMeDB: same traits plus `human_pathogen`, NCBI tax id)
3. `gut_microbiome_925_oxygen_pH_spore.csv` (oxygen/pH/spore for 925 species)
4. NCBI taxonomy `taxdump` (`names.dmp` + `nodes.dmp`) for `ncbi_tax_id`, `phylum`,
   `class_name`, `order_name`, `family`, `superkingdom`

Species rows match on the (genus, species) key produced by `gutdb.transform.taxon_key`.
A genus row is filled only where every species of that genus in the sources agrees on
the value — a unanimity vote, so a genus is never assigned a trait its members dispute.
`human_pathogen` is excluded from that vote: the source records pathogenicity only as
a positive assertion, so one pathogenic species would otherwise brand its whole genus.

Weighted by evidence rather than by row, the effect is largest where it matters:
phylum went from 58% to 99.6% of association rows, gram stain from 0% to 56%,
oxygen requirement from 30% to 62%.

### 3. Vocabularies normalized

Values already in the database were rewritten onto one spelling per concept:
`Anaerobe`/`anaerobic`/`Obligate anaerobe` → `anaerobe`, `Facultatively anaerobe` →
`facultative anaerobe`, and so on — oxygen requirement went from 25 distinct values
to 10, sporulation from 14 to 2. Under a controlled vocabulary an unlisted value is
kept but lowercased, so a category nobody anticipated survives while case variants
stop being separate values (`nanaerobe`, `aerotolerant` are still there).
`superkingdom` had collected three vocabularies — MiMeDB domains, the curated table's
kingdoms (`Eubacteria`, `Fungi`), and NCBI's post-2024 kingdom names (`Bacillati`,
`Pseudomonadati`) — and now folds onto `Bacteria` / `Archaea` / `Eukaryota` (3 values,
99.5% of taxa). The `Spore` column's mixed content was split at build time:
sporulation yes/no stays in `sporulation`, morphology text goes to `cell_arrangement`,
which is left as free descriptive text (41 values) since it is not a controlled field.

### 4. Sample and comparison context

- `samples.body_site` filled for all 33,630 rows — every source row is a stool
  sample, recorded in the GMrepo export but dropped by the loader.
- `phenotype_comparisons.method` completed for the last 45 comparisons, and their
  `notes` now record the source file and phenotype pair.

### What is still empty, and why

- **Sample demographics** — `bmi` (33,332 of 33,630 rows), `sex` (21,005),
  `age_years` (20,629), `country` (66), plus `qc_status` (2,743). Blank in the
  GMrepo export itself; filling them needs a GMrepo API pull, not a local source.
- **5 samples with no disease** (`PRJNA769284`) — that project contains both
  `Prostatic Neoplasms` (57) and `Health` (26) samples, so no project-level label can
  be assigned honestly. Left NULL.
- **317 taxa without an NCBI tax id** (220 genus, 97 species) and 23 without a domain
  — names the taxdump cannot resolve, mostly metagenome-derived or renamed clades
  (`Bulleidia p-1630-c5`, `Streptococcus luteciae`, `Eubacterium biforme`).
- **Traits for 1,691–3,286 taxa** — absent from all three trait sources, and their
  genus has no unanimous value either.
- **`ph_preference`** — the sources only ever record `alkaliphile` or `acidophile`;
  a missing value does not mean neutrophile.
- **`human_pathogen`** — 63 positives, the rest NULL by design
  (`gutdb.transform.normalize_pathogen_flag` never writes 0, which would fabricate
  negative evidence for taxa that simply were not annotated).
- **`p_value` / `q_value` / `effect_size`** — blank where the source is blank; the
  GMrepo export carries no significance columns.
- **`studies.title` for the 39 `PubMed` studies** — `project_title` was dropped from the
  marker CSVs; see the last section. `project_accession` (the PMID) still identifies them.

## Columns dropped from the mirror

A column-redundancy audit (`docs/redundant_columns.csv`) found five columns empty in
every row of every table this build can populate. They are gone from
`gutdb/schema_sqlite.sql` and from the SQLite file; `gutdb/schema.sql`,
`gutdb/pipeline.py` and `scripts/load_biomapai_study.py` are untouched, because four
of the five are live on the MySQL side.

| Column | Why it was always empty here |
|---|---|
| `samples.subject_id` | Written only by `scripts/load_biomapai_study.py` (`PRJNA1125469`, 2–3 timepoints per subject), a study this mirror does not load |
| `samples.timepoint` | Same loader, same reason |
| `samples.unclassified_fraction` | `pipeline.sync_gmrepo_abundances` computes it and writes it straight to MySQL, but it is not among that function's CSV export fieldnames, so no value reaches this build |
| `sample_taxon_abundances.detection_threshold` | `pipeline.py` inserts it from a source field the GMrepo export does not carry |
| `ingestion_runs.error_message` | Written only by `gutdb/pipeline.py` on the MySQL side. This mirror's `Loader.finish_run` takes no error argument and hardcodes `status = 'completed'`, so nothing local could ever fill it. All 20 run rows were `completed` with the column NULL. A load that dies part-way still shows up: its row stays at `status = 'running'`, which the `CHECK` constraint keeps alongside `failed` |

Consequences: `v_abundance_genus` and `v_abundance_species` no longer select
`subject_id`; `v_abundance_coverage` no longer reports `unclassified_fraction` (it was
NULL in all 12,049 rows). Row counts are unchanged in all 8 tables, `PRAGMA
foreign_key_check` and `integrity_check` are clean, all 7 query outputs return the same
row counts as before, and `VACUUM` took the file from 71.7 MB to 65.0 MB. Restore any
column here if its source becomes loadable — the header of `schema_sqlite.sql` records
where each one comes from.

## Columns dropped from the literature marker CSVs

The same audit checked the 15 hand-curated marker files (14 `data/literature_*_markers.csv`
plus `data/gmrepo_PRJNA705217_ibs_markers.csv`, 352 rows). Every comparison in them is
disease-versus-healthy, which made four columns pure restatement, and `project_title`
was dropped alongside them:

| Column | Evidence over the 352 rows |
|---|---|
| `phenotype_b` | equals `disease_name` in all 352 |
| `phenotype_b_mesh_id` | equals `mesh_id` in all 318 rows of the 14 files that carry both; the IBS file carries neither |
| `phenotype_a` | the constant `Health` in all 352 |
| `phenotype_a_mesh_id` | the constant `D006262` (Health) in all 352 |
| `project_title` | dropped by request; not redundant — see the consequence below |

1,730 cells removed. `positive_enriched_in` and `negative_enriched_in` were left alone:
`positive_enriched_in` matches `phenotype_b` in only 34 of 352 rows, so those two carry
the direction of each association, not a copy of the arms.

`Loader.comparison` in `scripts/build_gutdb_sqlite.py` and the comparison upsert in
`gutdb/pipeline.py` now default the control arm to `Health`/`D006262` and take the case
arm from `disease_name`/`mesh_id`. Neither default fires on the full GMrepo export, which
carries real phenotype pairs and supplies them; `comparison_key` is built from the same
two names as before, so its values are unchanged.

**Consequence of dropping `project_title`:** it was the only source of `studies.title`
for literature studies, and no other column reconstructs it (`project_id` is a bare PMID,
`citation` gives author/year/journal, `project_description` is a methods summary). Titles
are now NULL for the 39 `PubMed` studies, and `PRJNA705217` falls back to the title
`pipeline.py` synthesizes for GMrepo projects — `GMrepo curated Health vs. Irritable
Bowel Syndrome` in place of `Gut metagenomes of patients with post-infection IBS`.
`studies.title` is 273 of 312 non-null, down from 312.

Verified by building the database twice from the same loader, once with the original CSVs
and once with the trimmed ones: all 8 tables byte-identical except those 40 title cells,
same 4,374 / 161 / 312 / 378 / 16,390 / 33,630 / 1,243,990 row counts, and after the
cleanup pass the file matches the previously committed database in every table but
`studies.title`. `PRAGMA integrity_check` and `foreign_key_check` clean, all 7 queries
return their prior row counts.
