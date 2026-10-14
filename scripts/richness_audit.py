#!/usr/bin/env python3
"""How many parts does each profile actually have, and does that track the label?

613 (sample, rank) profiles in this database hold two taxa or fewer: 269 with a
single genus, 316 with two genera, 28 with two species. All 612 of them that
appear in the GMrepo CSVs under data/ match the source row for row, and GMrepo's
own numbers for them sum to exactly 1.0, so GMrepo reports these as COMPLETE
profiles rather than truncated ones. Nothing was lost at load and there is
nothing here to delete.

What is wrong is that nothing stops them reaching a model, and the shape of the
problem is not the one the single-taxon sweep found. Those 269 were failed
classifications -- mean unclassified_fraction 0.64 against 0.20 database-wide.
The two-taxon cohort sits at 0.27, which is close to normal, and the reason is
that thinness here is mostly a property of the STUDY, not of the sample:

  PRJNA1031953   101 samples   median  1 genus    max  4
  PRJNA832282     82 samples   median  2 genera   max  6
  PRJDB6133       58 samples   median  3 genera   max  7
  PRJNA422961    118 samples   median  3 genera   max  9
  PRJNA923180     27 samples   median  3 genera   max  7
  PRJNA889687    130 samples   median  3 genera   max 12

Six studies, 516 samples, in which NO sample clears 15 genera. A per-sample
filter cannot see this, because inside these studies a thin profile is not an
outlier -- it is the whole cohort.

THE HARM IS A CONFOUND THAT LOOKS LIKE A RESULT. Profile depth differs between
the arms of a cohort, and when a disease is carried by one thin study while its
controls are pooled from a hundred rich ones, depth predicts the label directly.
Neuromyelitis Optica is the clean case: its entire case arm is PRJNA422961 at a
median of 3 genera, its controls sit at 45, and a classifier on CLR genus
features scores AUC 0.995 -- while profile richness ALONE, a single integer that
is not a biological measurement, scores 0.982. The model is a richness detector.

This is NOT what scripts/loso_audit.py measures. Leave-one-study-out needs a
held-out study carrying both arms, and a single-study disease has none, so these
cohorts are exactly the ones LOSO silently declines to score and whose
random-split numbers therefore go unchallenged.

So this script reports the diagnostic rather than deleting anything:
richness_auc, the AUC of profile richness used by itself as the only predictor.
Read it as the share of a cohort's separability that needs no biology. Where a
study contributes both arms the same number is computed within study, which
removes the pooling and shows whether the depth difference is really about the
disease.

Writes data/richness_audit.csv.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.diagnostics import USABLE_GENERA, rank_auc, richness_confound

OUT = Path("data/richness_audit.csv")
ARMS_OUT = Path("data/richness_arm_imbalance.csv")
MIN_ARM = 20            # smaller than this and the AUC is noise
RANK = "genus"


def load(conn) -> pd.DataFrame:
    return pd.read_sql_query(
        f"""SELECT v.sample_id, v.phenotype, v.project_accession, v.data_type,
                   g.n_genera
            FROM v_ml_ready_samples v
            JOIN (SELECT a.sample_id, COUNT(*) AS n_genera
                  FROM sample_taxon_abundances a
                  JOIN taxa t ON t.id = a.taxon_id
                  WHERE t.taxonomic_rank = '{RANK}'
                  GROUP BY a.sample_id) g ON g.sample_id = v.sample_id
            WHERE v.ml_ready_deduplicated = 1""", conn)


def thin_profiles(conn) -> pd.DataFrame:
    """One row per (rank, richness); the unit counted is the profile, not the row.

    Grouping first and counting after is the whole point: a two-taxon profile is
    two rows in sample_taxon_abundances, so counting rows reports every such
    profile twice and turns 316 into 632.
    """
    return pd.read_sql_query(
        """SELECT taxonomic_rank, n_taxa, COUNT(*) AS n_profiles
           FROM (SELECT a.sample_id, t.taxonomic_rank, COUNT(*) AS n_taxa
                 FROM sample_taxon_abundances a
                 JOIN taxa t ON t.id = a.taxon_id
                 GROUP BY a.sample_id, t.taxonomic_rank
                 HAVING COUNT(*) <= 2) p
           GROUP BY taxonomic_rank, n_taxa
           ORDER BY taxonomic_rank, n_taxa""", conn)


def per_study(frame: pd.DataFrame) -> pd.DataFrame:
    g = frame.groupby("project_accession")["n_genera"]
    out = pd.DataFrame({"n_samples": g.size(), "median_genera": g.median(),
                        "max_genera": g.max()})
    return out[out["max_genera"] < USABLE_GENERA].sort_values("n_samples",
                                                              ascending=False)


def arm_imbalance(frame: pd.DataFrame) -> list[dict]:
    """Per-study depth ratios: the confound that no split can remove.

    See gutdb/diagnostics.py for why PRJNA769052 (2 genera against its own 50
    controls at 43) is a different kind of finding from PRJNA788674 (15 against
    30, in acute norovirus).
    """
    rows = []
    for (study, data_type), grp in frame.groupby(["project_accession", "data_type"],
                                                 dropna=False):
        for phenotype in sorted(set(grp["phenotype"]) - {"Health"}):
            arm = grp[grp["phenotype"].isin([phenotype, "Health"])]
            report = richness_confound(arm, richness="n_genera")
            for entry in report["per_study_depth"]:
                rows.append({"phenotype": phenotype, "data_type": data_type, **entry})
    return rows


def per_disease(frame: pd.DataFrame) -> list[dict]:
    """One row per disease: pooled and within-study richness AUC.

    Controls are the whole pooled Health arm, because that is what a model built
    off v_ml_ready_samples is given. The within-study figure then says whether
    the depth difference survives dropping the pooling.
    """
    health = frame[frame["phenotype"] == "Health"]
    rows = []
    for phenotype, arm in frame[frame["phenotype"] != "Health"].groupby("phenotype"):
        if len(arm) < MIN_ARM or len(health) < MIN_ARM:
            continue
        cohort = pd.concat([arm, health])
        report = richness_confound(cohort, richness="n_genera")
        rows.append({
            "phenotype": phenotype,
            "cases": len(arm),
            "controls": len(health),
            "case_studies": arm["project_accession"].nunique(),
            "median_genera_case": report["median_richness_case"],
            "median_genera_control": report["median_richness_control"],
            "richness_gap": (report["median_richness_case"]
                             - report["median_richness_control"]),
            "richness_auc_pooled": report["richness_auc_pooled"],
            "richness_auc_within_study": report["richness_auc_within_study"],
            "studies_with_both_arms": report["studies_with_both_arms"],
            "single_study_case_arm": int(arm["project_accession"].nunique() == 1),
        })
    return rows


def main() -> int:
    conn = connect(Settings.from_env(".env"))

    print("thin profiles, by rank (source-faithful: see the module docstring)")
    thin = thin_profiles(conn)
    for _, r in thin.iterrows():
        print(f"   {r['taxonomic_rank']:8} {int(r['n_taxa'])} taxa  "
              f"{int(r['n_profiles']):>5} profiles")
    print(f"   total {int(thin['n_profiles'].sum())} profiles")

    frame = load(conn)
    print(f"\nml_ready_deduplicated samples with a {RANK} profile: {len(frame)}")

    degenerate = per_study(frame)
    print(f"\nstudies where NO sample reaches {USABLE_GENERA} genera "
          f"({len(degenerate)} studies, {int(degenerate['n_samples'].sum())} "
          f"modelling-ready samples)")
    for acc, r in degenerate.iterrows():
        print(f"   {acc:16} {int(r['n_samples']):>5} samples  "
              f"median {r['median_genera']:>4.0f}  max {int(r['max_genera']):>3}")

    rows = per_disease(frame)
    rows.sort(key=lambda r: -(r["richness_auc_pooled"] or 0))
    print(f"\nrichness as the only predictor, per disease ({len(rows)} cohorts)")
    print(f"   {'phenotype':38} {'pooled':>7} {'within':>7} {'gap':>6} {'studies':>8}")
    for r in rows[:15]:
        within = r["richness_auc_within_study"]
        print(f"   {r['phenotype'][:38]:38} {r['richness_auc_pooled']:>7.3f} "
              f"{(f'{within:.3f}' if within is not None else 'n/a'):>7} "
              f"{r['richness_gap']:>+6.0f} {r['case_studies']:>8}")

    arms = [r for r in arm_imbalance(frame)
            if r["depth_ratio"] is not None and r["depth_ratio"] <= 0.5]
    arms.sort(key=lambda r: r["depth_ratio"])
    print(f"\nstudies whose own arms differ in depth by more than half "
          f"({len(arms)} of them) -- the confound no split can remove")
    print(f"   {'study':16} {'phenotype':30} {'case':>10} {'control':>8} {'ratio':>6}")
    for r in arms:
        print(f"   {r['project_accession']:16} {r['phenotype'][:30]:30} "
              f"{r['median_richness_case']:>6.0f} (n={r['n_case']:>3}) "
              f"{r['median_richness_control']:>8.0f} {r['depth_ratio']:>6.2f}")

    flagged = [r for r in rows if (r["richness_auc_pooled"] or 0) >= 0.75]
    single = [r for r in flagged if r["single_study_case_arm"]]
    print(f"\n{len(flagged)} cohorts where richness alone scores 0.75 or better; "
          f"{len(single)} of those have a single-study case arm, which is the "
          f"combination LOSO cannot score")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).sort_values("richness_auc_pooled", ascending=False
                                   ).to_csv(OUT, index=False)
    pd.DataFrame(arm_imbalance(frame)).sort_values("depth_ratio"
                                                   ).to_csv(ARMS_OUT, index=False)
    print(f"\nwrote {OUT} ({len(rows)} rows) and {ARMS_OUT}")
    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
