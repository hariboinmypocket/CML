"""Paired BMI experiments with study-held-out, training-only preprocessing.

Read either store without changing it. Database IDs are local join keys only;
model features use stable genus/species keys. Missing demographics never select
the cohort. Statistical imputations exist only inside fitted model pipelines.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.compose import ColumnTransformer
from sklearn.exceptions import ConvergenceWarning
from sklearn.impute import MissingIndicator, SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, OneHotEncoder, StandardScaler
from sklearn.utils.validation import check_is_fitted

from .compositional import clr, multiplicative_replacement


class InsufficientStudiesError(ValueError):
    def __init__(self, eligible_studies, required):
        self.eligible_studies = list(eligible_studies)
        self.required = required
        super().__init__(f"Only {len(self.eligible_studies)} eligible both-arm studies; need {required}")


class TrainingCLR(TransformerMixin, BaseEstimator):
    """Select taxa on training prevalence, then transform each sample separately.

    An abundance absent from a profiled sample is zero. A profile with no mass
    in the retained vocabulary gets a neutral CLR vector plus an explicit flag;
    it is counted in the report rather than removed from one model's test set.
    """
    def __init__(self, min_prevalence=0.10):
        self.min_prevalence = min_prevalence

    def fit(self, X, y=None):
        if not 0 < self.min_prevalence <= 1:
            raise ValueError("min_prevalence must be in (0, 1]")
        values = np.asarray(X, dtype=float)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError("Abundances must be finite and non-negative")
        self.n_features_in_ = values.shape[1]
        self.feature_names_in_ = np.asarray(X.columns, dtype=object)
        self.keep_ = (values > 0).mean(axis=0) >= self.min_prevalence
        if self.keep_.sum() < 2:
            raise ValueError("Fewer than two training-prevalent taxa; CLR comparison unsupported")
        self.selected_taxa_ = self.feature_names_in_[self.keep_]
        return self

    def transform(self, X):
        check_is_fitted(self, "keep_")
        # ColumnTransformer passes the training input vocabulary in fixed order.
        frame = X.loc[:, self.selected_taxa_].astype(float)
        if not np.isfinite(frame.to_numpy()).all() or (frame < 0).any().any():
            raise ValueError("Abundances must be finite and non-negative")
        mass = frame.sum(axis=1)
        empty = mass <= 0
        composition = frame.div(mass.where(~empty, 1), axis=0)
        composition.loc[empty, :] = 1 / frame.shape[1]
        transformed = clr(multiplicative_replacement(composition))
        return np.column_stack([transformed.to_numpy(), empty.astype(float)])

    def get_feature_names_out(self, input_features=None):
        check_is_fitted(self, "keep_")
        return np.append(self.selected_taxa_, "no_retained_taxon_mass")


def _categorical_object(X):
    return np.asarray(X, dtype=object)


def bmi_pipeline(taxon_columns, *, include_bmi=False, min_prevalence=0.10):
    numeric = ["age_years"] + (["bmi"] if include_bmi else [])
    preprocessing = ColumnTransformer([
        ("microbiome", Pipeline([
            ("clr", TrainingCLR(min_prevalence)),
            ("scale", StandardScaler()),
        ]), list(taxon_columns)),
        ("numeric", Pipeline([
            ("impute", SimpleImputer(strategy="median", keep_empty_features=True)),
            ("scale", StandardScaler()),
        ]), numeric),
        # Always emit indicators, even if a value first becomes missing in test.
        ("missingness", MissingIndicator(features="all"), numeric),
        ("sex", Pipeline([
            ("object", FunctionTransformer(_categorical_object, feature_names_out="one-to-one")),
            ("impute", SimpleImputer(strategy="constant", fill_value="unknown",
                                      keep_empty_features=True)),
            ("encode", OneHotEncoder(handle_unknown="ignore", sparse_output=False)),
        ]), ["sex"]),
    ], remainder="drop")
    return Pipeline([
        ("preprocess", preprocessing),
        ("classifier", LogisticRegression(C=1.0, class_weight="balanced",
                                           solver="lbfgs", max_iter=3000)),
    ])


def query_frame(connection, sql, params=(), *, sqlite=False):
    cursor = connection.cursor()
    try:
        cursor.execute(sql.replace("%s", "?") if sqlite else sql, tuple(params))
        return pd.DataFrame(cursor.fetchall(), columns=[d[0] for d in cursor.description])
    finally:
        cursor.close()


def load_model_cohort(connection, disease, *, sqlite=False, rank="genus", min_taxa=3):
    """Read one-rank profiles and metadata from the caller's pinned snapshot."""
    if rank not in {"genus", "species"} or min_taxa < 2:
        raise ValueError("Choose genus/species and min_taxa >= 2")
    if sqlite:
        columns = {r[1] for r in connection.execute("PRAGMA table_info(samples)")}
    else:
        columns = set(query_frame(connection, "SHOW COLUMNS FROM samples")["Field"])
    subject = "s.subject_id" if "subject_id" in columns else "NULL"
    meta = query_frame(connection, f"""SELECT s.id AS sample_id, s.run_accession,
        st.project_accession, s.gmrepo_sample_id AS biological_sample_id,
        {subject} AS subject_id, d.name AS phenotype,
        s.sex, s.age_years, s.bmi, s.body_site, s.qc_status
        FROM samples s JOIN studies st ON st.id=s.study_id
        JOIN diseases d ON d.id=s.disease_id
        WHERE d.name IN (%s, 'Health')
          AND st.id IN (SELECT sx.study_id FROM samples sx
                         JOIN diseases dx ON dx.id=sx.disease_id WHERE dx.name=%s)
        ORDER BY st.project_accession, s.run_accession""", (disease, disease), sqlite=sqlite)
    selection = {"candidate_rows": len(meta), "rank": rank, "min_taxa": min_taxa,
                 "subject_ids_available_in_store": "subject_id" in columns}
    meta = meta[meta.body_site.eq("feces") & ~meta.qc_status.astype(str).eq("0")].copy()
    selection["fecal_nonfailed_rows"] = len(meta)
    sample_ids = meta.sample_id.tolist()
    pieces = []
    for start in range(0, len(sample_ids), 1000):
        batch = sample_ids[start:start + 1000]
        placeholders = ",".join(["%s"] * len(batch))
        pieces.append(query_frame(connection, f"""SELECT a.sample_id,
            t.genus_key, t.species_key, a.relative_abundance
            FROM sample_taxon_abundances a JOIN taxa t ON t.id=a.taxon_id
            WHERE t.taxonomic_rank=%s AND a.sample_id IN ({placeholders})
              AND a.relative_abundance > 0""", (rank, *batch), sqlite=sqlite))
    if not pieces:
        raise ValueError("No case/control fecal samples for this disease")
    long = pd.concat(pieces, ignore_index=True)
    if long.empty:
        raise ValueError(f"No {rank} profiles for this disease")
    long["taxon_key"] = ["taxon:" + json.dumps([g, s], ensure_ascii=False, separators=(",", ":"))
                         for g, s in zip(long.genus_key, long.species_key)]
    richness = long.groupby("sample_id").taxon_key.nunique()
    meta = meta[meta.sample_id.isin(richness[richness >= min_taxa].index)].copy()
    selection["profile_eligible_rows_before_dedup"] = len(meta)
    # Lexicographically first eligible run, independent of either store's IDs.
    sample_key = meta.biological_sample_id.fillna("").astype(str)
    subject_key = meta.subject_id.fillna("").astype(str)
    meta["person_key"] = np.where(subject_key.ne(""), "subject:" + subject_key,
                                  np.where(sample_key.ne(""), "sample:" + sample_key,
                                           "run:" + meta.run_accession))
    meta = meta.drop_duplicates(["project_accession", "person_key"]).set_index("sample_id")
    selection["representative_rows"] = len(meta)
    selection["representatives_without_known_subject_id"] = int(meta.subject_id.isna().sum())
    selected = long[long.sample_id.isin(meta.index)]
    if selected.duplicated(["sample_id", "taxon_key"]).any():
        raise ValueError("Duplicate stable taxon keys within a sample; repair identity before modelling")
    raw = selected.pivot(index="sample_id", columns="taxon_key", values="relative_abundance")
    raw = raw.reindex(index=meta.index, columns=sorted(raw.columns)).fillna(0).astype(float)
    for field in ("age_years", "bmi"):
        meta[field] = pd.to_numeric(meta[field], errors="coerce").astype(float)
    meta["sex"] = meta.sex.map(lambda x: str(x).strip().lower() if pd.notna(x) else np.nan)
    meta.loc[~meta.sex.isin(["male", "female"]), "sex"] = np.nan
    meta["sex"] = meta.sex.astype(object)
    meta.loc[~meta.age_years.between(0, 110), "age_years"] = np.nan
    meta.loc[~meta.bmi.between(0, 200, inclusive="right"), "bmi"] = np.nan
    return meta, raw, selection


def mask_uncertain_demographics(meta, provenance_path: str | Path):
    """Mask flagged existing imports for this experiment, not in the database."""
    meta = meta.copy()
    path = Path(provenance_path)
    counts = Counter()
    if not path.exists():
        return meta, {"report_available": False, "masked_cells": {}}
    report = pd.read_csv(path, dtype=str).fillna("")
    for row in report.itertuples(index=False):
        if not row.quality_flags or row.field not in {"sex", "age_years", "bmi"}:
            continue
        mask = meta.run_accession.eq(row.run_id)
        if row.field == "sex":
            mask &= meta.sex.eq(row.value)
        else:
            mask &= np.isclose(meta[row.field], float(row.value), rtol=0, atol=1e-6,
                               equal_nan=False)
        counts[row.field] += int(mask.sum())
        meta.loc[mask, row.field] = np.nan
    return meta, {"report_available": True, "report_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                  "masked_cells": dict(counts)}


def evaluate_bmi(meta, raw, *, disease, min_arm=5, min_studies=3, min_prevalence=0.10,
                 progress=None):
    """Evaluate both versions on identical leave-one-study-out folds."""
    if min_arm < 1 or min_studies < 3:
        raise ValueError("Require min_arm >= 1 and at least three studies")
    y_all = meta.phenotype.eq(disease).astype(int)
    counts = y_all.groupby(meta.project_accession).agg(["sum", "count"])
    eligible = counts[(counts["sum"] >= min_arm) & ((counts["count"] - counts["sum"]) >= min_arm)]
    rejected = counts.loc[~counts.index.isin(eligible.index)]
    if len(eligible) < min_studies:
        raise InsufficientStudiesError(eligible.index, min_studies)
    meta = meta[meta.project_accession.isin(eligible.index)].copy()
    raw = raw.loc[meta.index]
    X = pd.concat([raw, meta[["sex", "age_years", "bmi"]]], axis=1)
    y = meta.phenotype.eq(disease).astype(int).to_numpy()
    groups = meta.project_accession.to_numpy()
    rows, predictions = [], []
    for train, test in LeaveOneGroupOut().split(X, y, groups):
        study = str(groups[test][0])
        if set(groups[train]) & set(groups[test]):
            raise AssertionError("Study leakage")
        row = {"held_out_study": study, "n_train": len(train), "n_test": len(test),
               "test_cases": int(y[test].sum()), "test_controls": int((1-y[test]).sum()),
               "train_bmi_observed": int(meta.iloc[train].bmi.notna().sum()),
               "test_bmi_observed": int(meta.iloc[test].bmi.notna().sum()),
               "test_sex_missing": int(meta.iloc[test].sex.isna().sum()),
               "test_age_missing": int(meta.iloc[test].age_years.isna().sum())}
        probabilities = {}
        for name, include_bmi in (("without_bmi", False), ("with_bmi", True)):
            model = bmi_pipeline(raw.columns, include_bmi=include_bmi, min_prevalence=min_prevalence)
            with warnings.catch_warnings():
                # An unconverged fold must not be silently reported as a score.
                warnings.simplefilter("error", ConvergenceWarning)
                model.fit(X.iloc[train], y[train])
            p = model.predict_proba(X.iloc[test])[:, 1]
            probabilities[name] = p
            row[f"auc_{name}"] = float(roc_auc_score(y[test], p))
            row[f"balanced_accuracy_{name}"] = float(balanced_accuracy_score(y[test], p >= 0.5))
            row[f"brier_{name}"] = float(brier_score_loss(y[test], p))
            pre = model.named_steps["preprocess"]
            transform = pre.named_transformers_["microbiome"].named_steps["clr"]
            row["taxa_selected_on_training"] = len(transform.selected_taxa_)
            row["test_profiles_without_retained_mass"] = int(
                (raw.iloc[test].loc[:, transform.selected_taxa_].sum(axis=1) <= 0).sum())
            row[f"training_numeric_medians_{name}"] = (
                pre.named_transformers_["numeric"].named_steps["impute"].statistics_.tolist())
        row["auc_change_with_bmi"] = row["auc_with_bmi"] - row["auc_without_bmi"]
        rows.append(row)
        for i, index in enumerate(test):
            predictions.append({"run_accession": str(meta.iloc[index].run_accession),
                                "held_out_study": study, "label": int(y[index]),
                                "without_bmi": float(probabilities["without_bmi"][i]),
                                "with_bmi": float(probabilities["with_bmi"][i])})
        if progress:
            progress(row)
    summary = {"disease": disease, "samples": len(meta), "cases": int(y.sum()),
               "controls": int((1-y).sum()), "held_out_studies": len(rows),
               "missing_sex_retained": int(meta.sex.isna().sum()),
               "missing_age_retained": int(meta.age_years.isna().sum()),
               "missing_bmi_retained": int(meta.bmi.isna().sum()),
               "excluded_studies_without_minimum_both_arms": rejected.index.tolist()}
    for name in ("without_bmi", "with_bmi"):
        summary[f"macro_study_auc_{name}"] = float(np.mean([r[f"auc_{name}"] for r in rows]))
    summary["macro_auc_change_with_bmi"] = (
        summary["macro_study_auc_with_bmi"] - summary["macro_study_auc_without_bmi"])
    summary["default_model"] = "without_bmi"
    # Refit on the complete eligible training cohort for future research use.
    models = {}
    for name, include_bmi in (("without_bmi", False), ("with_bmi", True)):
        models[name] = bmi_pipeline(raw.columns, include_bmi=include_bmi,
                                    min_prevalence=min_prevalence).fit(X, y)
    return {"summary": summary, "folds": rows, "held_out_predictions": predictions}, {
        "models": models, "metadata": meta, "abundances": raw, "disease": disease,
        "feature_keys": raw.columns.tolist(), "label_definition": {0: "Health", 1: disease},
    }
