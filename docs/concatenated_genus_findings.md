# Concatenated two-genus taxa rows: analysis before the repair

Status: **the taxa repair is applied** (`scripts/fix_concatenated_genus_taxa.py
--apply`, 2026-10-14: 62 merged, 14 renamed, 3 survivors renamed, 167 trait
values moved, 0 skipped, audit in `data/concatenated_genus_audit.csv`). **The
loader is NOT fixed yet** -- `gutdb/transform.py` has the replacement resolver
`resolve_source_binomial`, but `_collapse_mimedb` still calls the old path, so
the next MiMeDB ingest re-creates these names. The loader section below is the
live record for that remaining step.

## Scope confirmed: 69 rows, and a defensible superset

The detection rule from the brief (rank='species', `species` contains a space,
first token is a genus scientific name in the local taxdump and differs from the
`genus` column, skipping `Candidatus ...`) flags **70** rows. One is the false
positive the brief names -- `Kocuria marina subsp. indica` (id 1529911) -- and it
is separated from the real cases by a check worth keeping: `genus + ' ' + first
token` ("Kocuria marina") resolves in the taxdump, so the first token is this
row's own epithet and the rest is infraspecific. No other row in the table
behaves that way, leaving exactly **69**.

A second, better-founded definition of scope came out of the loader work: the
rows the *fixed* loader would no longer produce. Computed from the 2,648 MiMeDB
source rows, that set contains the 69 plus classes the taxdump-genus rule cannot
see:

- 4 rows where the repeated genus arrived behind the `Candidatus` prefix, so the
  stored name duplicates one genus instead of naming two:
  `Methanomassiliicoccus candidatus methanomassiliicoccus intestinalis` (2170),
  `Stoquefichus candidatus stoquefichus massiliensis` (2252),
  `Blastococcus candidatus blastococcus massiliensis` (3426),
  `Nitrososphaera candidatus nitrososphaera gargensis` (3731).
- rows whose stale genus is a *synonym* genus name, absent from the taxdump as a
  genus-rank scientific name: `Actinoplanes ampullariella violaceochromogenes`
  (3455), `Propionibacterium pseudopropionibacterium propionicum` (3468).
- rows that are fabricated names but not two-*genus* names, e.g.
  `Unclassified burkholderiales burkholderiales bacterium` (2815) and
  `Lymphocryptovirus human gammaherpesvirus 4` (3185). These were deliberately
  excluded from the loader fix (see guards below), so they stay out of scope and
  the loader keeps producing them -- called out rather than silently dropped.

Repairing the loader without repairing a row it would re-key means the next
ingest inserts a fresh duplicate beside the repaired row, so the repair set and
the loader change have to be decided together.

## The brief's premise about dependents does not hold

> Many of these rows DO carry associations and abundances, unlike the single
> orphan already fixed.

Measured against the live database: **all 69 carry zero** `sample_taxon_abundances`,
zero `taxon_disease_associations` and zero `taxon_disease_adjudication` rows.
Verified both by subquery across the set and by direct per-row counts on
2219, 3679, 3071, 2152, 3391, 3685, 3225, 2251.

What is true is the splitting the brief describes: **52 of the 69** have a
correctly-named row that holds the evidence (`Phocaeicola bacteroides vulgatus`
beside row 1569 `Phocaeicola vulgatus` with 10,263 abundances; `Blautia
ruminococcus gnavus` beside row 3826 `Mediterraneibacter gnavus` with 8,637).
So the repair is a merge-or-rename decision on `taxa` only: no evidence gets
relabelled, no abundance summing is exercised, no rank group is renormalised.
The summing path from `fix_candidatus_taxa.py` should still be implemented --
the database changes hourly and a dependent could appear between the dry run and
`--apply` -- but on today's data it will not fire.

## Per-row classification

Resolution is by tax ID, as the precedent requires: strain tax IDs walk up to
the species node, and the node's scientific name is the target. 67 of 69 resolve
cleanly and the repaired name reproduces the organism MiMeDB's own `name` field
describes. Two guards found every case needing a hand decision:

1. **Same-species-node invariant** (for merges: the canonical row's tax ID must
   resolve to the same species node as the mangled row's). Violations: 3319 and
   3630 only, both because a *canonical* row carries a tax ID absent from the
   taxdump (1328388, 1232666 -- deleted strain nodes).
2. **Epithet drift** (stored epithet vs NCBI's epithet for the tax ID, which is
   the F. prausnitzii trap `merge_duplicate_taxa.py`'s OVERRIDES exist for).
   Three rows: 2861, 3391, 3681.

The five hand decisions:

| row | reading |
| --- | --- |
| 2861 `Propionibacterium porphyromonas gingivalis` | MiMeDB record MMDBm0002385's `name` ("Propionibacterium freudenreichii subsp. shermanii CIRM-BIA1") and `ncbi_tax_id` (754252) agree; only `species` says *Porphyromonas gingivalis*. Its traits (Gram **positive**, facultative, bacilli) are P. freudenreichii's, not P. gingivalis's, which independently shows the `species` cell is the corrupt one. Merge into row 1608 `Propionibacterium freudenreichii`. **Do not** touch P. gingivalis. |
| 3391 `Actinomyces schaalia odontolytica` | tax 411466 is strain ATCC 17982, which NCBI has moved to *Schaalia dentiphila* (species 3050224), while the species *S. odontolytica* (1660, type strain ATCC 17929) still stands and MiMeDB's `name` and `species` both say odontolytica. Exactly F. prausnitzii: retaxid to the species the row names (1660) and merge into row 1773 `Schaalia odontolytica` (1,503 abundances). Traits came from ATCC 17982, now a different species, so **suppress the trait backfill** for this row. |
| 3630 `Staphylococcus mammaliicoccus sciuri` | tax 1232666 is not in the taxdump; MiMeDB's `name` ("Staphylococcus sciuri subsp. sciuri Z8") and `species` ("Mammaliicoccus sciuri") agree on the epithet. Merge into row 1306 `Mammaliicoccus sciuri` (tax 1296). |
| 3681 `Lachnoclostridium clostridium methoxybenzovorans` | not drift: names.dmp records *Clostridium methoxybenzovorans* as a synonym of 69825 *Lacrimispora indolis*, and strain SR3 (1286820) sits under it. Merge into row 1213, as 3684 does. |
| 3319 `Hafnia klebsiella pneumoniae` | canonical row 1370 carries dead strain node 1328388; names match exactly and the mangled row's 1328337 resolves to species 573. Merge. |

Two further rows need no override but do need the plan built **target-first**:
`Bacillus alkalihalobacillus pseudofirmus` (2139, strain tax 398511) and
`Alkalihalobacillus bacillus pseudofirmus` (2630, species tax 79885) resolve to
the same target with no existing row to merge into, so they collide with each
other. 2630 already carries the species tax ID and should be the survivor.

And a rule, not an override, covers two more: where no row holds the target key
but exactly one row holds the target **tax ID** under a stale name, that row is
the same organism and the mangled row merges into it, with the survivor renamed
to NCBI's current name -- `merge_duplicate_taxa.py`'s canonicalisation. That is
3649 into 1248 (`Lactobacillus hilgardii` -> `Lentilactobacillus hilgardii`,
tax 1588) and 3717 into 157 (`Selenomonas ruminantium` ->
`Pseudoselenomonas ruminantium`, tax 971). Both survivors have no dependents.

## The SQLite mirror holds none of these rows

The mirror is built from the curated CSV and GMrepo, not MiMeDB: 0 rows with
`source_database='MiMeDB'`, and none of the mangled keys present. The repair is
MySQL-only. Keep the mirror pass anyway (keyed on `(genus_key, species_key)`,
never `taxa.id`) so it reports "absent" rather than assuming.

## Loader fix: what the data forced

The brief's suggestion -- prefer the genus the `species`/`name` field implies
when that field is a full binomial -- is the right shape, but two premises in it
do not survive contact with the source:

- **MiMeDB does not consistently keep the current name in `species`.** It is
  current in `Lactobacillus` + `Limosilactobacillus reuteri`, and *stale* in
  `Phocaeicola` + `Bacteroides vulgatus`, `Lactiplantibacillus` +
  `Lactobacillus plantarum`, `Mediterraneibacter` + `Ruminococcus torques`.
  MMDBm0000292 is the reverse of the brief's model: genus `Fannyhessea` is
  current, `species` "atopobium vaginae" is stale.
- So no offline rule picks the *current* genus. The achievable goal is weaker and
  still worth having: **never store a name for an organism that does not exist.**
  Taking the binomial as written always yields a real name -- sometimes a
  synonym, which `merge_duplicate_taxa.py` already folds by tax ID. The
  alternative (keep the `genus` column, drop the foreign token) fabricates:
  `Clavispora auris`, `Alfipia birgiae`, `Nialla hominis`, `Brochotrix
  thermosphacta`. Rejected for that reason.

Verified over all 2,648 source rows. Unguarded, the rule changes 194 keys and
has two unwanted effects: it would start ingesting ~40 rows the loader currently
drops (`if genus:` skips rows with an empty `genus` column), and it would turn
order and family names into genera (`Burkholderiales bacterium`, `Lachnospiraceae
bacterium`) and `Human gammaherpesvirus 4` into genus "Human". Guarded -- keep
the `if genus:` skip, require the first token to be alphabetic, >= 4 letters and
not a rank suffix (-aceae/-ales/-ineae/-idae/-mycetes/-bacteria...), require the
second token to be a lowercase epithet that is not `bacterium`/`archaeon`/
`symbiont` and contains neither "virus" nor "phage", and route `Candidatus`
through `parse_scientific_name` -- it changes **143 of 2,648** keys and nothing
else. Of the 115 distinct new keys, **84 land on an existing row** (an update, no
new row) and **0 are genuinely new organisms**.

The residue to be honest about: **31 new keys land on no row while their tax ID
is already held by another row**, so a post-repair ingest would add that many
stale-synonym rows (`Bacteroides vulgatus` for tax 821 beside `Phocaeicola
vulgatus`). That is name-variant drift the existing dedup script is built to
clear, not fabricated corruption -- but it is a recurrence, and the README does
not currently document `merge_duplicate_taxa.py` as a post-ingest step. Two ways
to close it, both unimplemented: have `upsert_taxon` update the single existing
row that already carries the row's `ncbi_tax_id` instead of inserting beside it
(cheap, DB-only, but does not help strain tax IDs), or let the loader resolve
names from `GUTDB_TAXDUMP` when present (matches the canon, but makes ingest
non-reproducible between machines with and without the dump).

## Still to do

1. Point `_collapse_mimedb` at `resolve_source_binomial` (both call sites -- it
   derives the key once for grouping and again per group). `upsert_taxon`
   re-derives the epithet with `extract_species_epithet`, which is idempotent on
   an already-resolved pair, so it needs no change.
2. Decide what to do about the 31 keys that would land on no row while their tax
   ID is already held: either have `upsert_taxon` update the single row already
   carrying the tax ID, or document `merge_duplicate_taxa.py` as a post-ingest
   step in the README, which it currently is not.
3. One trait moved by the merge is wrong at source and now sits on a row that
   matters: MiMeDB records `gram: Positive` for both of its *Bacteroides
   vulgatus* records (MMDBm0000056/57), so row 1569 `Phocaeicola vulgatus`
   (10,263 abundances) now reads Gram-positive in a phylum where 107 rows read
   Negative and MiMeDB's own sibling records (massiliensis, plebeius, dorei,
   coprophilus) all say Negative. The value is MiMeDB's claim, faithfully moved
   by the backfill rule, but relocating it from an orphan to an
   evidence-carrying row is what makes it visible. Row 216 *Parabacteroides
   distasonis* carries the same error independently of this repair. Worth a
   gram-stain-against-phylum plausibility check rather than a one-cell edit.

Baseline before the repair: taxa 5,897; abundance_rows 3,073,876; associations
16,889; `abundance_rank_groups_not_summing_to_1` 0; `disagrees_unapproved` 0.
After: taxa **5,835** (62 merged away), abundances and associations unchanged at
3,073,876 and 16,889 (no row in scope carried any), rank groups not summing to 1
still **0**, `disagrees_unapproved` still **0** in both stores, and re-running
the detection finds nothing -- the repair is idempotent.

Row 3230 was already `Sphingobacterium paucimobilis` in the live database before
this work: the precedent fix was applied from another worktree, whose branch is
not an ancestor of this one.
