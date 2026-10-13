#!/usr/bin/env python3
"""Export a model-ready CLR feature matrix for one disease against Health.

    python scripts/export_ml_matrix.py "Crohn Disease"
    python scripts/export_ml_matrix.py "Crohn Disease" --rank=species --prevalence=0.2

Writes data/ml_matrix_<disease>_<rank>.csv with one row per sample, the label and
covariates first, then one CLR column per taxon, plus a .json of how it was made.

Sample selection uses v_ml_ready_samples' ml_ready_deduplicated flag, so every
row has sex, age and a feature vector, and repeated runs of one subject are
collapsed to the first. Training on a sample that appears twice leaks an
individual across the train/test split, which flatters a model without improving
it.

The CLR columns are produced by gutdb.compositional, in the order that keeps them
meaningful: filter taxa by prevalence, renormalize, replace zeros
multiplicatively, then transform. The matrix is written out rather than stored
because a CLR value only means something relative to the feature set it was
computed over -- see that module's docstring.

THE ROWS COME FROM MANY STUDIES AND THAT IS A CONFOUNDER, not a detail. A Crohn
matrix draws on 105 of them, each with its own protocol, sequencing run and
recruitment, and cohort membership is correlated with the label because a study
contributes mostly cases or mostly controls. A model given these features and a
random train/test split will partly learn which study a sample came from and
report a flattering score for it. project_accession is written into the matrix
for that reason: group the cross-validation on it, or carry it as a batch
covariate. Everything else in this database that compares arms does so within a
study for the same reason.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.compositional import clr_matrix
from gutdb.config import Settings
from gutdb.db import connect

OUT_DIR = Path("data")


def arg(flag: str, default):
    for item in sys.argv[1:]:
        if item.startswith(f"--{flag}="):
            return item.split("=", 1)[1]
    return default


def main() -> int:
    positional = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not positional:
        print(__doc__)
        return 2
    disease = positional[0]
    rank = str(arg("rank", "genus"))
    prevalence = float(arg("prevalence", 0.10))

    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    cur.execute(
        """SELECT sample_id, phenotype, sex, age_years, bmi, country,
                  project_accession, data_type
           FROM v_ml_ready_samples
           WHERE ml_ready_deduplicated = 1 AND phenotype IN (%s, 'Health')""",
        (disease,))
    rows = cur.fetchall()
    cur.close()
    if not rows:
        print(f"no deduplicated ML-ready samples for {disease!r}")
        conn.close()
        return 1

    import pandas as pd
    meta = pd.DataFrame(rows, columns=["sample_id", "phenotype", "sex", "age_years",
                                       "bmi", "country", "project_accession",
                                       "data_type"]).set_index("sample_id")
    n_case = int((meta["phenotype"] != "Health").sum())
    n_ctrl = int((meta["phenotype"] == "Health").sum())
    print(f"{disease}: {n_case:,} cases, {n_ctrl:,} controls, "
          f"{meta['project_accession'].nunique()} studies")
    if n_case == 0 or n_ctrl == 0:
        print("need both arms; nothing to export")
        conn.close()
        return 1

    matrix, provenance = clr_matrix(conn, rank, list(meta.index),
                                    min_prevalence=prevalence)
    conn.close()
    if matrix.empty:
        print("no taxa survived the prevalence filter")
        return 1

    cur2 = None
    names = {}
    conn2 = connect(Settings.from_env(".env"))
    c = conn2.cursor()
    ph = ",".join(["%s"] * len(matrix.columns))
    c.execute(f"SELECT id, scientific_name FROM taxa WHERE id IN ({ph})",
              tuple(matrix.columns))
    names = dict(c.fetchall())
    c.close(); conn2.close()
    matrix.columns = [names.get(t, f"taxon_{t}") for t in matrix.columns]

    out = meta.join(matrix, how="inner")
    slug = re.sub(r"[^a-z0-9]+", "_", disease.lower()).strip("_")
    path = OUT_DIR / f"ml_matrix_{slug}_{rank}.csv"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(path)

    per_study = meta.groupby("project_accession")["phenotype"].apply(
        lambda v: float((v != "Health").mean()))
    single_arm = int(((per_study == 0) | (per_study == 1)).sum())
    provenance.update({
        "disease": disease, "n_cases": n_case, "n_controls": n_ctrl,
        "studies": int(meta["project_accession"].nunique()),
        "studies_contributing_one_arm_only": single_arm,
        "selection": "v_ml_ready_samples.ml_ready_deduplicated = 1",
        "column_order": "metadata, then one CLR column per taxon",
        "warning": ("rows pool multiple studies and study is confounded with the "
                    "label; group cross-validation on project_accession or carry "
                    "it as a batch covariate"),
    })
    path.with_suffix(".json").write_text(json.dumps(provenance, indent=2) + "\n")

    print(f"   taxa kept at prevalence >= {prevalence:.0%}: {provenance['taxa_kept']}")
    print(f"   median mass retained by that filter: "
          f"{provenance['median_mass_retained_by_filter']:.1%}")
    print(f"   matrix: {out.shape[0]} samples x {matrix.shape[1]} CLR features")
    print(f"   !! {single_arm} of {meta['project_accession'].nunique()} studies "
          f"contribute only one arm; study is confounded with the label.")
    print("      Group cross-validation on project_accession.")
    print(f"   wrote {path} and {path.with_suffix('.json').name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
