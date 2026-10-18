# BMI comparison without discarding incomplete demographics

The default model excludes BMI. A paired alternative adds BMI and a BMI-missing
indicator. Both use CLR microbiome features, age, an age-missing indicator, and
sex with an unknown category. The comparison uses exactly the same samples and
held-out studies for both versions. Source tables and CSVs are never changed.

Install the optional modelling dependencies if needed:

```bash
python3 -m pip install -r requirements-ml.txt
```

Run IBS versus Health on MySQL:

```bash
python3 scripts/compare_bmi_models.py
```

Run the identical analysis code against the SQLite mirror:

```bash
python3 scripts/compare_bmi_models.py --sqlite gut_microbiome.sqlite
```

Other disease targets and experimental settings can be supplied explicitly:

```bash
python3 scripts/compare_bmi_models.py "Crohn Disease" \
  --rank genus --prevalence 0.10 --output-dir results/crohn_bmi_run_1
```

The default output is `results/bmi_comparison/<store>/<disease>_<rank>/`:

- `report.md`: paired per-study scores, missingness retained, and limitations.
- `comparison.json`: the pinned census, selection counts, settings, software
  versions, masking audit, per-fold metrics and imputer values, and held-out
  predictions keyed by run accession.
- `models_and_inputs.joblib`: both full-cohort refits and the exact selected
  metadata and abundance inputs. These inputs pin the experiment even when
  agents later modify the live database. Binary model bundles are gitignored.

Existing output directories are never overwritten. Use a new `--output-dir`
for another run. At least three eligible studies are required. If fewer exist,
the command exits with status 2 and writes a `not_scored_insufficient_studies`
report without AUCs or a model bundle. A missing study in the mirror is not
silently synthesized from MySQL, and the two stores' model results are not
treated as interchangeable.

## Selection and preprocessing

Each store is opened in a read-only transaction, and counts are pinned before
reading the cohort. MySQL uses a repeatable-read consistent snapshot; SQLite
uses a read transaction with `query_only`. SQLite is never rewritten, so there
is no destructive write needing a backup in this workflow.

Use only fecal profiles, exclude explicit QC failures, and require at least
three detected taxa at the chosen rank (configurable with `--min-taxa`). NULL
QC is retained where sources supply usable profiles without that field; this
is not a claim that those samples passed GMrepo QC. Genus and species layers
are never mixed or changed. Absent taxon rows become zero within an observed
rank profile, and prevalence denominators include every training sample.

Select the lexicographically first eligible run for each known subject within
a study. Without a subject ID, fall back to the biological sample accession,
then the run accession. This is a deterministic representative, not necessarily
a baseline visit. Unknown subject IDs remain a limitation: distinct sample
accessions may still belong to one person. Nothing is deleted from the stores.

Require at least five cases and five controls per study after that selection.
Entire studies are held out with LeaveOneGroupOut. Only studies meeting the
both-arm requirement participate; this avoids pooling controls from unrelated
cohorts. The legacy `ml_ready` flags are not used because they require complete
sex and age and would remove the samples this experiment needs to retain.

For every fold, fit the taxon prevalence filter using training data only;
renormalize within the retained rank, replace zeros multiplicatively, and
compute CLR. Fit numeric medians, scaling, and sex encoding on the training
fold only. All-missing numeric training features use zero plus an explicit
missingness indicator; their zero is an internal computational placeholder,
never written as a source measurement. Profiles with no mass in the training
taxon vocabulary receive neutral CLR values plus a flag and are counted in the
report. Both models retain those same test rows.

When `data/gmrepo_demographic_provenance.csv` exists, flag-matched uncertain
values are masked in the model inputs only: untraced historical values,
inferred sample links, and censored ages represented as exact values. Masking
requires the run, field and current value to match the report, so an unrelated
later correction is not blanked. This report covers only the earlier cleanup's
filled cells and is not a complete provenance certification of all covariates.

## Interpretation

The estimator is fixed L2 logistic regression (`C=1`, balanced class weights).
No hyperparameter search or LDA association features are used. Primary scores
are paired AUCs averaged with equal weight per held-out study. Balanced
accuracy and Brier scores are also saved. The BMI model adds its missingness
indicator as well as its numeric value, so any improvement might reflect
reporting patterns rather than BMI biology.

Keep the no-BMI model as the default. A small difference across a handful of
studies is exploratory, not evidence of a clinically useful model or a stable
BMI benefit. Choosing features using these folds requires a separate external
test set for confirmation. Saved full-cohort refits are research models; their
training scores are not used as performance estimates.

References: [scikit-learn leakage guidance](https://scikit-learn.org/stable/common_pitfalls.html#data-leakage),
[imputation](https://scikit-learn.org/stable/modules/impute.html), and
[LeaveOneGroupOut](https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.LeaveOneGroupOut.html).
