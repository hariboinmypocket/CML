#!/usr/bin/env python3
"""Compare a no-BMI baseline with adding BMI, holding out entire studies.

No database writes. Default: IBS, genus, MySQL. --sqlite reads the mirror using
the same cohort and model code. Missing age/sex/BMI are retained and imputed
only inside training folds. Existing complete-demographics ML views are unused.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import re
import sqlite3
import sys

import joblib
import numpy as np
import pandas as pd
import sklearn

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.modeling import (InsufficientStudiesError, evaluate_bmi, load_model_cohort,
                            mask_uncertain_demographics, query_frame)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("disease", nargs="?", default="Irritable Bowel Syndrome")
    parser.add_argument("--env-file", default=".env")
    parser.add_argument("--sqlite", type=Path, metavar="PATH")
    parser.add_argument("--rank", choices=["genus", "species"], default="genus")
    parser.add_argument("--prevalence", type=float, default=0.10)
    parser.add_argument("--min-taxa", type=int, default=3)
    parser.add_argument("--min-arm", type=int, default=5)
    parser.add_argument("--min-studies", type=int, default=3)
    parser.add_argument("--provenance", type=Path,
                        default=ROOT / "data/gmrepo_demographic_provenance.csv")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    is_sqlite = args.sqlite is not None
    store = "sqlite" if is_sqlite else "mysql"
    slug = re.sub(r"[^a-z0-9]+", "_", args.disease.lower()).strip("_")
    out = args.output_dir or ROOT / "results/bmi_comparison" / store / f"{slug}_{args.rank}"
    # Require a fresh output directory; do not overwrite a previous experiment.
    if out.exists():
        parser.error(f"Output directory already exists: {out}; choose --output-dir")
    if is_sqlite:
        connection = sqlite3.connect(args.sqlite.resolve().as_uri() + "?mode=ro", uri=True)
        connection.execute("PRAGMA query_only=ON")
        connection.execute("BEGIN")
    else:
        connection = connect(Settings.from_env(args.env_file))
        connection.start_transaction(readonly=True, consistent_snapshot=True,
                                     isolation_level="REPEATABLE READ")
    snapshot = {"store": store, "read_started_utc": datetime.now(timezone.utc).isoformat()}
    try:
        if not is_sqlite:
            snapshot["database_snapshot_utc"] = str(query_frame(connection, "SELECT UTC_TIMESTAMP(6)").iloc[0, 0])
        snapshot["census"] = {table: int(query_frame(connection, f"SELECT COUNT(*) n FROM {table}",
                                                    sqlite=is_sqlite).iloc[0, 0])
                              for table in ["samples", "studies", "taxa", "sample_taxon_abundances"]}
        print("Pinned read-only snapshot:", json.dumps(snapshot), flush=True)
        meta, raw, selection = load_model_cohort(connection, args.disease, sqlite=is_sqlite,
                                                rank=args.rank, min_taxa=args.min_taxa)
    finally:
        connection.rollback()
        connection.close()
    meta, uncertainty = mask_uncertain_demographics(meta, args.provenance)
    print(f"Eligible profiles before study-arm selection: {len(meta):,}; "
          f"studies: {meta.project_accession.nunique()}", flush=True)

    def progress(row):
        print(f"Held out {row['held_out_study']}: "
              f"AUC no BMI {row['auc_without_bmi']:.3f}; with BMI {row['auc_with_bmi']:.3f}; "
              f"observed test BMI {row['test_bmi_observed']}/{row['n_test']}", flush=True)

    try:
        result, bundle = evaluate_bmi(meta, raw, disease=args.disease, min_arm=args.min_arm,
                                       min_studies=args.min_studies, min_prevalence=args.prevalence,
                                       progress=progress)
    except InsufficientStudiesError as error:
        out.mkdir(parents=True, exist_ok=False)
        failure = {"status": "not_scored_insufficient_studies", "disease": args.disease,
                   "snapshot": snapshot, "selection": selection,
                   "eligible_studies": error.eligible_studies,
                   "required_studies": error.required, "reason": str(error),
                   "uncertainty_masking": uncertainty}
        (out / "comparison.json").write_text(json.dumps(failure, indent=2) + "\n")
        (out / "report.md").write_text(
            f"# {args.disease}: BMI comparison not scored\n\n"
            f"Store: {store}. Read snapshot: {snapshot['read_started_utc']}.\n\n"
            f"{error}. The minimum was not relaxed.\n\n"
            "Eligible studies: " + ", ".join(error.eligible_studies) + ".\n\n"
            "No model was fitted and no AUC is reported. No database records were modified.\n")
        print(f"NOT SCORED: {error}. Snapshot and reason saved to {out}", flush=True)
        return 2
    digest = hashlib.sha256()
    for frame in [bundle["metadata"].reset_index(drop=True), bundle["abundances"].reset_index(drop=True)]:
        digest.update(pd.util.hash_pandas_object(frame, index=False).to_numpy().tobytes())
        digest.update(json.dumps(list(frame.columns)).encode())
    result.update(status="completed", snapshot=snapshot, selection=selection, uncertainty_masking=uncertainty,
                  input_sha256=digest.hexdigest(),
                  settings={"rank": args.rank, "prevalence": args.prevalence,
                            "min_taxa": args.min_taxa, "min_arm": args.min_arm,
                            "min_studies": args.min_studies},
                  software={"python": platform.python_version(), "sklearn": sklearn.__version__,
                            "numpy": np.__version__, "pandas": pd.__version__},
                  interpretation=[
                      "Paired exploratory comparison; study averages weight studies equally.",
                      "No rows are excluded because sex, age, or BMI is missing.",
                      "Only same-study case/control fecal cohorts meeting the minimum arm sizes are scored.",
                      "Taxon prevalence, medians, categorical encoding and scaling are fitted on training folds only.",
                      "The BMI version adds both BMI and its missingness indicator; any change is not necessarily a biological BMI effect.",
                      "Entirely missing training numeric features receive zero plus an explicit missing indicator.",
                      "Unknown subject IDs prevent verifying independence of every participant.",
                      "Absent retained-taxon mass uses neutral CLR plus a flag and is counted per fold.",
                      "Model selection on these folds requires a new independent test set for confirmation.",
                  ])
    out.mkdir(parents=True, exist_ok=False)
    (out / "comparison.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    bundle.update(snapshot=snapshot, settings=result["settings"], input_sha256=result["input_sha256"])
    joblib.dump(bundle, out / "models_and_inputs.joblib", compress=3)
    s = result["summary"]
    report = [f"# {args.disease}: paired BMI comparison", "",
              f"Store: {store}. Snapshot: {snapshot['read_started_utc']}.", "",
              f"{s['samples']:,} representatives: {s['cases']:,} cases and {s['controls']:,} controls, "
              f"{s['held_out_studies']} held-out studies.", "",
              f"Retained missing demographics: sex {s['missing_sex_retained']:,}, "
              f"age {s['missing_age_retained']:,}, BMI {s['missing_bmi_retained']:,}.", "",
              "| Held-out study | Test cases / controls | Observed test BMI | AUC without BMI | AUC with BMI | Change |",
              "|---|---:|---:|---:|---:|---:|"]
    for r in result["folds"]:
        report.append(f"| {r['held_out_study']} | {r['test_cases']} / {r['test_controls']} | "
                      f"{r['test_bmi_observed']} / {r['n_test']} | {r['auc_without_bmi']:.3f} | "
                      f"{r['auc_with_bmi']:.3f} | {r['auc_change_with_bmi']:+.3f} |")
    report += ["", f"Equal-study mean AUC: **{s['macro_study_auc_without_bmi']:.3f} without BMI**, "
               f"**{s['macro_study_auc_with_bmi']:.3f} with BMI** "
               f"(change {s['macro_auc_change_with_bmi']:+.3f}).", "",
               "Default remains the model without BMI. This is an exploratory comparison, not a clinical validation.", "",
               "## Method and limits", ""]
    report.extend("- " + text for text in result["interpretation"])
    report += ["", "All database access was read only; no tables, source values, ranks, or QC labels were modified.",
               "The saved bundle contains the exact selected inputs and both full-cohort refits. "
               "Reported scores come exclusively from held-out predictions in comparison.json.", ""]
    (out / "report.md").write_text("\n".join(report))
    print(json.dumps(s, indent=2), flush=True)
    print(f"Saved experiment: {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
