# CKD demographic cleanup

Applied 2026-10-16 to project PRJNA562327, using the pinned source checks in
`results/ckd_link_verification/verification.json`.

The original audit contained 489 clinical-to-run mappings. Only 57 HzCKD source
identifiers match public ENA aliases after case/separator normalization. The
273 HC and 159 CKD mappings relied on unverified clinical-to-library bridges,
even though their demographic values match the publication.

## Changes

- Active `data/ckd_study_demographics_audit.csv`: 489 → 57 rows. The surviving
  join basis explicitly describes the source-identifier match.
- MySQL `gut_microbiome.samples` and `data/gmrepo_project_samples.csv`: cleared
  sex, age_years and bmi for the 432 withdrawn mappings (1,296 cells per store).
- SQLite `samples`: those 432 runs already had NULL demographics and remain
  NULL. Filled 57 verified sexes, 57 ages and 41 BMIs; 16 source-missing BMIs
  remain NULL. Keys are run accessions plus checked project identity, not IDs.
- Refreshed active demographic provenance: withdrawn fills no longer appear.
  Historical logs and originals are retained in the backup, not active tables.
- No samples, abundance observations, disease associations or taxa were deleted.
  No new database tables were created. This was a CKD-specific repair, not a
  claim that all other projects' demographics have been independently verified.

Backup directory:
`data/archive/ckd_verified_cleanup_20261016T164854514525Z/`.
It contains the pre-write SQLite snapshot, original CSVs, the verification
manifest, targeted MySQL rows, and a step-by-step recovery report. Final
machine-readable changes and before/after census are in
`results/ckd_link_verification/cleanup.json`.

## Safe reruns

`python3 scripts/harvest_ckd_study.py` is now a verification-gated dry run against
both stores and the active CSVs. It uses the hash-checked original audit in the
cleanup backup, never guesses library prefixes, fills verified gaps only, and
preserves conflicting independent measurements. `--apply` makes a new backup
before writing. The original unsafe prefix-inference importer is retired.

The underlying command also accepts an explicit source snapshot:

```sh
python3 scripts/repair_ckd_demographics.py --audit-source data/archive/ckd_verified_cleanup_20261016T164854514525Z/ckd_study_demographics_audit.csv
```

Withdrawn values must not be restored from historical CSVs without an explicit
participant-to-library/BioSample/run crosswalk. A matching numerical suffix,
matching population distributions, or an old `documented` label is not proof.
Fresh independently verified demographic recovery remains allowed.
