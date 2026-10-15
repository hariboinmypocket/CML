#!/usr/bin/env python3
"""Give every sample the platform that produced it, not its study's summary.

studies.data_type answers at the wrong granularity. GMrepo records the platform
per RUN -- 'AMPLICON' or 'Metagenomics' in data/gmrepo_project_samples.csv -- and
10 of its projects contain both, so the loader collapsed them to 'mixed':

    PRJNA763023   200 shotgun  1038 amplicon
    PRJNA888860   230           264
    PRJNA604850   134           354
    PRJNA648321    48           428
    PRJEB42155    208           211
    PRJNA905444    88           236
    PRJNA813705    60            60
    PRJEB21933     75            40
    PRJNA730640    40            32
    PRJEB53891     20            20

3,786 samples -- 1,103 shotgun and 2,683 amplicon -- whose platform the source
knew for every single one.

WHY THE STUDY-LEVEL VALUE IS NOT GOOD ENOUGH. Platform sets the taxonomic floor:
amplicon sequencing reads one marker gene, which is near-identical between
species of a genus, so it resolves to genus and no further. Inside a 'mixed'
study some samples therefore carry species-rank profiles and some cannot, with
nothing recording which. Platform is also the strongest batch covariate in this
database -- scripts/loso_audit.py measures platform mixing at +0.175 optimism
against +0.060 for the rest -- and both that audit and export_ml_matrix.py read
the STUDY's data_type, so for these 3,786 samples the batch covariate they
carried was the literal word 'mixed', which is not a batch.

Filled from two sources, finest first, and never guessed:

  1. The run's own data_type in data/gmrepo_project_samples.csv, matched on
     run_accession and folded through normalize_platform, which already maps
     GMrepo's vocabulary ('AMPLICON' -> '16S', 'Metagenomics' -> 'mNGS').
  2. The study's data_type, but ONLY where it names a single platform. A study
     reading 'mixed' or 'meta-analysis' tells a sample nothing, so those stay
     NULL rather than inheriting a word that is not a platform.

Both stores. Idempotent. Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import csv
import os
import sqlite3
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.transform import (NOT_A_PLATFORM, derive_study_platform,
                             normalize_platform)

DRY_RUN = "--apply" not in sys.argv
SQLITE = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")
SOURCE = Path("data/gmrepo_project_samples.csv")
# NOT_A_PLATFORM and derive_study_platform come from gutdb.transform so the two
# stores and the loaders share one definition.


def platform_by_run() -> dict[str, str]:
    """run_accession -> normalized platform, from GMrepo's per-run export."""
    out: dict[str, str] = {}
    if not SOURCE.exists():
        print(f"   {SOURCE} not present; per-run platforms unavailable")
        return out
    with SOURCE.open(newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            run = (row.get("run_id") or "").strip()
            platform = normalize_platform(row.get("data_type"))
            if run and platform not in NOT_A_PLATFORM:
                out[run] = platform
    return out


def fill(label, cur, placeholder, by_run) -> Counter:
    counts: Counter = Counter()
    cur.execute("""SELECT s.id, s.run_accession, st.data_type
                   FROM samples s JOIN studies st ON st.id = s.study_id""")
    rows = cur.fetchall()
    updates: list[tuple[str, int]] = []
    for sample_id, run, study_type in rows:
        platform = by_run.get((run or "").strip())
        if platform:
            counts["from_run"] += 1
        else:
            study_platform = normalize_platform(study_type)
            if study_platform in NOT_A_PLATFORM:
                counts["left_null"] += 1
                continue
            platform = study_platform
            counts["from_study"] += 1
        updates.append((platform, sample_id))

    if not DRY_RUN and updates:
        cur.executemany(f"UPDATE samples SET sequencing_platform = {placeholder} "
                        f"WHERE id = {placeholder}", updates)
    print(f"   {label}: {counts['from_run']} from the run, "
          f"{counts['from_study']} from a single-platform study, "
          f"{counts['left_null']} left NULL")
    return counts


def reconcile_study_data_type(label, cur, placeholder) -> Counter:
    """Derive studies.data_type from the samples, so the two cannot disagree.

    The mirror's studies.data_type held FOUR vocabularies for two platforms --
    'AMPLICON' (191) beside '16S rRNA' (28), 'Metagenomics' (82) beside 'shotgun
    metagenomics' (9) -- while MySQL held the normalized '16S'/'mNGS'. So
    `WHERE data_type = '16S'` returned 236 studies against one store and ZERO
    against the other, on the column that governs taxonomic resolution.

    The mirror was also actively wrong, not merely differently spelled: it had no
    'mixed' at all, so each of the 10 projects holding both platforms was labelled
    as one of them.

    Deriving it removes the class of problem rather than this instance. A study
    whose samples all share a platform gets that platform; one whose samples
    disagree gets 'mixed', which is now a summary of a recorded fact rather than
    a value that replaced one. Studies with no samples are left alone -- an
    association-only study has nothing to derive from, and 'meta-analysis' is a
    true statement about those.
    """
    counts: Counter = Counter()
    cur.execute("""SELECT st.id, st.data_type, s.sequencing_platform
                   FROM studies st
                   JOIN samples s ON s.study_id = st.id
                   WHERE s.sequencing_platform IS NOT NULL
                   GROUP BY st.id, st.data_type, s.sequencing_platform""")
    seen: dict[int, tuple[str, set]] = {}
    for study_id, stored, platform in cur.fetchall():
        seen.setdefault(study_id, (stored, set()))[1].add(platform)
    for study_id, (stored, platforms) in seen.items():
        derived = derive_study_platform(platforms)
        if not derived:
            continue
        if derived == stored:
            counts["already_right"] += 1
            continue
        counts[f"{stored} -> {derived}"] += 1
        if not DRY_RUN:
            cur.execute(f"UPDATE studies SET data_type = {placeholder} "
                        f"WHERE id = {placeholder}", (derived, study_id))
    # A study with no samples has nothing to DERIVE from, but its stored value
    # can still be the raw source vocabulary, and normalizing that is not a
    # derivation. The mirror held AMPLICON (59 studies), '16S rRNA' (28) and
    # 'shotgun metagenomics' (9) this way, all of them sample-less, so
    # `WHERE data_type = '16S'` matched 122 studies there against 236 in MySQL.
    # normalize_platform leaves 'mixed' and 'meta-analysis' alone, so nothing
    # gains a platform it has no evidence for.
    cur.execute("""SELECT st.id, st.data_type FROM studies st
                   WHERE NOT EXISTS (SELECT 1 FROM samples s
                                     WHERE s.study_id = st.id)""")
    for study_id, stored in cur.fetchall():
        folded = normalize_platform(stored)
        if not folded or folded == stored:
            continue
        counts[f"{stored} -> {folded} (no samples)"] += 1
        if not DRY_RUN:
            cur.execute(f"UPDATE studies SET data_type = {placeholder} "
                        f"WHERE id = {placeholder}", (folded, study_id))

    print(f"   {label}: studies.data_type derived from the samples")
    for key, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"      {key:34} {n:>5}")
    return counts


def report(label, cur) -> None:
    cur.execute("""SELECT st.data_type, s.sequencing_platform, COUNT(*)
                   FROM samples s JOIN studies st ON st.id = s.study_id
                   GROUP BY st.data_type, s.sequencing_platform
                   ORDER BY st.data_type, COUNT(*) DESC""")
    print(f"   {label}: study data_type -> per-sample platform")
    for study_type, platform, n in cur.fetchall():
        flag = "   <- recovered" if study_type == "mixed" and platform else ""
        print(f"      {str(study_type):16} {str(platform):10} {n:>7}{flag}")


def main() -> int:
    by_run = platform_by_run()
    print(f"per-run platforms available from the source: {len(by_run)}")

    print("\nmysql:")
    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    fill("mysql", cur, "%s", by_run)
    if not DRY_RUN:
        conn.commit()
    reconcile_study_data_type("mysql", cur, "%s")
    if not DRY_RUN:
        conn.commit()
    report("mysql", cur)
    conn.close()

    if Path(SQLITE).exists():
        print("\nsqlite mirror:")
        sq = sqlite3.connect(SQLITE)
        sc = sq.cursor()
        fill("sqlite", sc, "?", by_run)
        if not DRY_RUN:
            sq.commit()
        reconcile_study_data_type("sqlite", sc, "?")
        if not DRY_RUN:
            sq.commit()
        report("sqlite", sc)
        sq.close()
    else:
        print(f"\nsqlite mirror {SQLITE} not present; skipped")

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
