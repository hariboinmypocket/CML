# CKD demographic-link verification

Verification snapshot: 2026-10-16T16:39:45.391272+00:00. Source audit SHA-256: `e9c1dc167449763cc6b8ab715bfce57f1a4f56cbae3cfaf1c3f76410b4a04c0a`.

All 489 audit records' demographic values match Supporting Excel 1. ENA confirms all recorded library-to-run links.

For the 273 previously inferred HC mappings, exact-run SRA metadata confirms the same libraries and BioSample accessions, but none independently reports age, sex or BMI. No complete clinical HC identifier matched a public sample alias.

The original 21 supplementary spreadsheets and two supporting PDFs were inspected. Supporting Excel 1 has a clinical sample_number field, but no crosswalk was found from it to sequencing accessions. ENA contains overlapping healthy-control numbering series Mobio_ZZ_HC_1..160 and zzHC_1..320. Matching suffixes or the larger range does not verify participant identity.

**Verdict: the 273 HC participant-to-run assignments remain unverified.** Do not promote them to confirmed demographics or copy them into SQLite on this evidence. An explicit author-provided participant-to-BioSample/library/run crosswalk is needed.

The other 216 records were previously labelled documented. Only the 57 HzCKD aliases match the source clinical IDs after case/separator normalization. The 159 CKD records also rely on an institutional-prefix assumption; the old label is not independent proof.

No database tables, sample values, audit rows, or mapping-status labels were changed.

Sources: [article](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7578882/fullTextXML), [Supporting Excel 1](https://rcastoragev2.blob.core.windows.net/5de81978ee2b164cb591ed3fede2bb56/ADVS-7-2001936-s003.xls), [ENA run metadata](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA562327&result=read_run&fields=run_accession,sample_accession,secondary_sample_accession,experiment_accession,study_accession,library_name,sample_alias,sample_title,experiment_alias,experiment_title,study_alias&format=tsv&download=false). Full per-record checks and source hashes are in verification.json.
