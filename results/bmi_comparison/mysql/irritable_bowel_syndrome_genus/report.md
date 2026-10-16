# Irritable Bowel Syndrome: paired BMI comparison

Store: mysql. Snapshot: 2026-10-16T15:52:36.419056+00:00.

887 representatives: 574 cases and 313 controls, 5 held-out studies.

Retained missing demographics: sex 613, age 613, BMI 762.

| Held-out study | Test cases / controls | Observed test BMI | AUC without BMI | AUC with BMI | Change |
|---|---:|---:|---:|---:|---:|
| MullerMM_MARS_IBS_2020 | 51 / 24 | 75 / 75 | 0.654 | 0.651 | -0.003 |
| PRJEB37924 | 338 / 164 | 0 / 502 | 0.600 | 0.602 | +0.002 |
| PRJNA637763 | 85 / 26 | 0 / 111 | 0.466 | 0.460 | -0.006 |
| PRJNA682378 | 71 / 78 | 0 / 149 | 0.452 | 0.445 | -0.006 |
| PRJNA705217 | 29 / 21 | 50 / 50 | 0.494 | 0.493 | -0.002 |

Equal-study mean AUC: **0.533 without BMI**, **0.530 with BMI** (change -0.003).

Default remains the model without BMI. This is an exploratory comparison, not a clinical validation.

## Method and limits

- Paired exploratory comparison; study averages weight studies equally.
- No rows are excluded because sex, age, or BMI is missing.
- Only same-study case/control fecal cohorts meeting the minimum arm sizes are scored.
- Taxon prevalence, medians, categorical encoding and scaling are fitted on training folds only.
- The BMI version adds both BMI and its missingness indicator; any change is not necessarily a biological BMI effect.
- Entirely missing training numeric features receive zero plus an explicit missing indicator.
- Unknown subject IDs prevent verifying independence of every participant.
- Absent retained-taxon mass uses neutral CLR plus a flag and is counted per fold.
- Model selection on these folds requires a new independent test set for confirmation.

All database access was read only; no tables, source values, ranks, or QC labels were modified.
The saved bundle contains the exact selected inputs and both full-cohort refits. Reported scores come exclusively from held-out predictions in comparison.json.
