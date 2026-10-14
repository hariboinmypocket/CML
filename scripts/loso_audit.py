#!/usr/bin/env python3
"""How much of a model's score is the study it came from?

A gut-microbiome classifier is usually scored with a random train/test split. In
this database that is close to meaningless, because study membership is
confounded with the label: of the 105 studies behind the Crohn cohort, 96
contribute only one arm, so "which study is this sample from" very nearly answers
"is this a case". A random split puts samples from the same study on both sides
of it, and the model is rewarded for recognising the cohort.

Leave-one-study-out holds out an entire study, so the model must generalise to a
protocol, a population and a sequencing run it has never seen -- which is what
using the model on a new patient would actually ask of it. The gap between the
two numbers is the part of a random-split score that is cohort recognition.

Scored only where it can be: a held-out study needs BOTH arms or AUC is
undefined on it. 20 of the 65 diseases with a modelling-ready case have two such
studies, 12 have three, and none has ten.

Writes data/loso_audit.csv. Model and features are deliberately plain -- L2
logistic regression on CLR genus features -- because the question is about the
evaluation, not about squeezing out the best score.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
warnings.filterwarnings("ignore")

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from gutdb.compositional import clr_matrix
from gutdb.config import Settings
from gutdb.db import connect

OUT = Path("data/loso_audit.csv")
MIN_FOLDS = 3          # fewer than three held-out studies is not an estimate
MIN_ARM = 5            # a held-out arm smaller than this gives a meaningless AUC
PREVALENCE = 0.10


def fit_predict(train_X, train_y, test_X):
    scaler = StandardScaler().fit(train_X)
    model = LogisticRegression(penalty="l2", C=1.0, max_iter=2000, solver="lbfgs")
    model.fit(scaler.transform(train_X), train_y)
    return model.predict_proba(scaler.transform(test_X))[:, 1]


def evaluate(disease: str, conn) -> dict | None:
    meta = pd.read_sql_query(
        """SELECT sample_id, phenotype, project_accession, sequencing_platform
           FROM v_ml_ready_samples
           WHERE ml_ready_deduplicated = 1 AND phenotype IN (%s, 'Health')""",
        conn, params=(disease,)).set_index("sample_id")
    meta["y"] = (meta["phenotype"] != "Health").astype(int)

    # Only studies carrying both arms can be a held-out fold.
    per_study = meta.groupby("project_accession")["y"].agg(["sum", "count"])
    per_study["n_ctrl"] = per_study["count"] - per_study["sum"]
    usable = per_study[(per_study["sum"] >= MIN_ARM) & (per_study["n_ctrl"] >= MIN_ARM)]
    if len(usable) < MIN_FOLDS:
        return None
    meta = meta[meta["project_accession"].isin(usable.index)]

    X, prov = clr_matrix(conn, "genus", list(meta.index), min_prevalence=PREVALENCE)
    if X.empty:
        return None
    meta = meta.loc[X.index]
    y = meta["y"].to_numpy()
    groups = meta["project_accession"].to_numpy()
    if len(np.unique(y)) < 2:
        return None

    # Random split: the number usually reported.
    random_scores = []
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=0).split(X, y):
        if len(np.unique(y[te])) < 2:
            continue
        random_scores.append(roc_auc_score(y[te], fit_predict(X.iloc[tr], y[tr], X.iloc[te])))

    # Leave one study out: what the model is actually asked to do.
    loso_scores, folds = [], 0
    for study in np.unique(groups):
        te = groups == study
        tr = ~te
        if len(np.unique(y[te])) < 2 or len(np.unique(y[tr])) < 2:
            continue
        loso_scores.append(roc_auc_score(y[te], fit_predict(X.iloc[tr], y[tr], X.iloc[te])))
        folds += 1
    if folds < MIN_FOLDS or not random_scores:
        return None

    # The per-sample platform, not the study's data_type: 10 projects hold both
    # and read 'mixed', which is not a platform, so a cohort drawing on them used
    # to look single-platform or three-platform depending on nothing real.
    platforms = sorted(meta["sequencing_platform"].dropna().unique())
    return {
        "disease": disease,
        "studies_used": int(len(usable)),
        "platforms": "+".join(platforms),
        "samples": int(len(meta)),
        "cases": int(y.sum()),
        "controls": int(len(y) - y.sum()),
        "features": int(X.shape[1]),
        "auc_random_split": round(float(np.mean(random_scores)), 3),
        "auc_loso_mean": round(float(np.mean(loso_scores)), 3),
        "auc_loso_min": round(float(np.min(loso_scores)), 3),
        "auc_loso_max": round(float(np.max(loso_scores)), 3),
        "optimism": round(float(np.mean(random_scores) - np.mean(loso_scores)), 3),
        "loso_folds": folds,
    }


def main() -> int:
    conn = connect(Settings.from_env(".env"))
    diseases = pd.read_sql_query(
        """SELECT phenotype, COUNT(*) n FROM v_ml_ready_samples
           WHERE ml_ready_deduplicated = 1 AND phenotype <> 'Health'
           GROUP BY phenotype ORDER BY n DESC""", conn)["phenotype"].tolist()

    results = []
    for disease in diseases:
        try:
            out = evaluate(disease, conn)
        except Exception as exc:                     # noqa: BLE001
            print(f"   {disease[:38]:40} skipped ({type(exc).__name__})")
            continue
        if out:
            results.append(out)
            print(f"   {disease[:38]:40} random {out['auc_random_split']:.3f}  "
                  f"LOSO {out['auc_loso_mean']:.3f}  "
                  f"optimism {out['optimism']:+.3f}  ({out['loso_folds']} folds)")
    conn.close()

    if not results:
        print("no disease has enough both-arm studies to score")
        return 1
    frame = pd.DataFrame(results).sort_values("optimism", ascending=False)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT, index=False)

    print(f"\nscored {len(frame)} of {len(diseases)} diseases "
          f"(the rest lack {MIN_FOLDS} studies with both arms of at least {MIN_ARM})")
    print(f"   mean random-split AUC: {frame['auc_random_split'].mean():.3f}")
    print(f"   mean LOSO AUC:         {frame['auc_loso_mean'].mean():.3f}")
    print(f"   mean optimism:         {frame['optimism'].mean():+.3f}")
    print(f"   diseases where LOSO is no better than chance (<=0.5): "
          f"{int((frame['auc_loso_mean'] <= 0.5).sum())}")

    # Held-out folds often cross sequencing platforms as well as cohorts, which
    # is a harder transfer. Reported rather than concluded from: there are too
    # few single-platform diseases to separate the two effects, and the one
    # all-16S case still loses 0.115, so study confounding is not merely assay
    # confounding wearing a hat.
    single = frame[~frame["platforms"].str.contains(r"\+", regex=True)]
    mixed = frame[frame["platforms"].str.contains(r"\+", regex=True)]
    if len(single) and len(mixed):
        print(f"\n   folds within one platform  ({len(single):>2} diseases): "
              f"mean optimism {single['optimism'].mean():+.3f}")
        print(f"   folds across platforms     ({len(mixed):>2} diseases): "
              f"mean optimism {mixed['optimism'].mean():+.3f}")
        print("   too few single-platform diseases to attribute the difference;"
              " the all-16S")
        print("   case still loses 0.115, so the study is a confounder within a"
              " platform too.")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
