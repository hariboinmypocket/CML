#!/usr/bin/env python3
"""Record which scale each cohort's metabolite levels are on.

sample_metabolite_levels mixes scales that differ in KIND, not magnitude, and
nothing recorded which:

    cohort                             levels      neg     zero          max
    YACHIDA_CRC_2019                  130,242        0   78,400     665,180
    KIM_ADENOMAS_2020                  89,760   46,514      513         3.3
    ERAWIJANTARI_GASTRIC_CANCER_2020   48,960        0   33,677     634,438
    SINHA_CRC_2016                     46,374   23,553      281         7.6
    POYET_BIO_ML_2019                  35,594        0        0   3.3e9
    KOSTIC_INFANTS_DIABETES_2015       13,471        0       13   3.07e10
    KANG_AUTISM_2017                    2,596        0      407       461.4

Two cohorts arrive already log- or z-transformed -- SINHA holds 23,553 negative
levels against 22,540 positive ones with a mean of -0.04 -- while the other
eleven are non-negative with maxima spanning nine orders of magnitude.

THE ZERO IS THE ACTUAL HAZARD. On a non-negative scale a zero means the compound
was not detected: 78,400 of YACHIDA's 130,242 levels, 60% of them. On a centred
scale a zero sits mid-distribution and means "at this cohort's average", which is
the OPPOSITE of absence. The same stored number meant opposite things across
studies with nothing to tell them apart, and the standing rule to normalize
within study does not cover it -- a uniform log1p over this table returns NaN on
70,067 rows and takes the log of a log on 136,134 more.

Derived from each study's own level range rather than asserted, because a
negative level cannot be a concentration and is therefore positive evidence of
an upstream transform. The absence of negatives is NOT evidence of a unit, only
of non-negativity, which is exactly what the other label says. The source records
units nowhere; studies.metabolomics_method carries the protocol DOI.

MySQL only; the SQLite mirror has no metabolite tables.
Idempotent. Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.transform import derive_metabolite_scale, zero_means_absence

DRY_RUN = "--apply" not in sys.argv


def main() -> int:
    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()

    cur.execute("""SELECT st.id, st.project_accession, st.metabolite_scale,
                          MIN(l.level), MAX(l.level), COUNT(*),
                          SUM(l.level = 0), SUM(l.level < 0)
                   FROM studies st
                   JOIN samples s ON s.study_id = st.id
                   JOIN sample_metabolite_levels l ON l.sample_id = s.id
                   GROUP BY st.id, st.project_accession, st.metabolite_scale""")
    rows = cur.fetchall()
    print(f"studies carrying metabolite levels: {len(rows)}\n")
    print(f"{'cohort':42} {'scale':14} {'levels':>8} {'zeros':>7} "
          f"{'zero means':>12}")

    changed = 0
    for study_id, accession, stored, lo, hi, n, zeros, negatives in sorted(
            rows, key=lambda r: -r[5]):
        scale = derive_metabolite_scale(lo, hi)
        meaning = "not detected" if zero_means_absence(scale) else "mid-scale"
        print(f"{accession[:42]:42} {scale:14} {n:>8} {zeros:>7} {meaning:>12}")
        if scale and scale != stored:
            changed += 1
            if not DRY_RUN:
                cur.execute("UPDATE studies SET metabolite_scale = %s WHERE id = %s",
                            (scale, study_id))

    if not DRY_RUN:
        conn.commit()
    print(f"\n{changed} studies updated")

    # The count that says why this matters: zeros whose meaning was ambiguous.
    cur.execute("""SELECT st.metabolite_scale, COUNT(*) n, SUM(l.level = 0) zeros
                   FROM sample_metabolite_levels l
                   JOIN samples s ON s.id = l.sample_id
                   JOIN studies st ON st.id = s.study_id
                   GROUP BY st.metabolite_scale""")
    print("\nlevels by scale, and how many are zero:")
    for scale, n, zeros in cur.fetchall():
        meaning = ("not detected" if zero_means_absence(str(scale))
                   else "at the cohort mean")
        print(f"   {str(scale):14} {n:>7} levels, {zeros:>7} zeros "
              f"-> read as {meaning}")
    conn.close()

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
