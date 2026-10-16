"""Prevent leakage, complete-case selection, and mixing taxonomy ID spaces."""
import sqlite3
import unittest

import numpy as np
import pandas as pd

from gutdb.modeling import (InsufficientStudiesError, TrainingCLR, bmi_pipeline,
                            evaluate_bmi, load_model_cohort)


class PreprocessingTests(unittest.TestCase):
    def test_taxon_selection_does_not_see_test_only_taxa(self):
        train = pd.DataFrame({"a": [0.8, 0.2], "b": [0.2, 0.8], "test_only": [0, 0]})
        transformer = TrainingCLR(0.5).fit(train)
        test = pd.DataFrame({"a": [0], "b": [0], "test_only": [1]})
        self.assertEqual(transformer.selected_taxa_.tolist(), ["a", "b"])
        transformed = transformer.transform(test)
        np.testing.assert_allclose(transformed, [[0, 0, 1]], atol=1e-12)
        self.assertEqual(transformer.selected_taxa_.tolist(), ["a", "b"])

    def test_baseline_is_independent_of_bmi_and_keeps_missing_covariates(self):
        X = pd.DataFrame({"a": [.9, .8, .2, .1], "b": [.1, .2, .8, .9],
                          "age_years": [20., 40., np.nan, np.nan],
                          "sex": ["female", np.nan, "male", np.nan],
                          "bmi": [np.nan] * 4})
        model = bmi_pipeline(["a", "b"]).fit(X, [0, 0, 1, 1])
        expected = model.predict_proba(X)
        changed = X.copy()
        changed["bmi"] = [10., 200., 80., 40.]
        np.testing.assert_allclose(model.predict_proba(changed), expected)
        test = X.copy()
        test["age_years"] = 1000.
        model.predict_proba(test)
        imputer = model.named_steps["preprocess"].named_transformers_["numeric"].named_steps["impute"]
        np.testing.assert_allclose(imputer.statistics_, [30.])
        alternative = bmi_pipeline(["a", "b"], include_bmi=True).fit(X, [0, 0, 1, 1])
        self.assertTrue(np.isfinite(alternative.predict_proba(X)).all())
        self.assertEqual(alternative.predict_proba(X).shape, (4, 2))

    def test_whole_missing_demographic_studies_are_retained_in_paired_folds(self):
        records, abundances = [], []
        for study in ["A", "B", "C"]:
            for i in range(8):
                label = i >= 4
                records.append(dict(run_accession=f"{study}_{i}", project_accession=study,
                                    phenotype="Disease" if label else "Health",
                                    sex=np.nan, age_years=np.nan, bmi=np.nan))
                abundances.append([0.8 if label else .2, .2 if label else .8])
        meta = pd.DataFrame(records)
        raw = pd.DataFrame(abundances, columns=["a", "b"])
        result, bundle = evaluate_bmi(meta, raw, disease="Disease", min_arm=2)
        self.assertEqual(result["summary"]["samples"], 24)
        self.assertEqual(result["summary"]["missing_age_retained"], 24)
        self.assertEqual(len(result["held_out_predictions"]), 24)
        self.assertEqual(len({r["run_accession"] for r in result["held_out_predictions"]}), 24)
        self.assertEqual({r["held_out_study"] for r in result["folds"]}, {"A", "B", "C"})
        for fold in result["folds"]:
            self.assertEqual((fold["n_train"], fold["n_test"]), (16, 8))
            self.assertAlmostEqual(fold["auc_without_bmi"], fold["auc_with_bmi"])
        self.assertTrue(bundle["metadata"].age_years.isna().all())

    def test_insufficient_studies_do_not_produce_a_score(self):
        meta = pd.DataFrame({"phenotype": ["Disease", "Health"] * 2,
                             "project_accession": ["A", "A", "B", "B"]})
        with self.assertRaises(InsufficientStudiesError) as caught:
            evaluate_bmi(meta, pd.DataFrame(index=meta.index), disease="Disease", min_arm=1)
        self.assertEqual(caught.exception.eligible_studies, ["A", "B"])


def fixture(offset):
    c = sqlite3.connect(":memory:")
    c.executescript("""
        CREATE TABLE samples(id INTEGER, run_accession TEXT, study_id INTEGER,
            disease_id INTEGER, gmrepo_sample_id TEXT, subject_id TEXT,
            sex TEXT, age_years REAL, bmi REAL, body_site TEXT, qc_status TEXT);
        CREATE TABLE studies(id INTEGER, project_accession TEXT);
        CREATE TABLE diseases(id INTEGER, name TEXT);
        CREATE TABLE taxa(id INTEGER, genus_key TEXT, species_key TEXT, taxonomic_rank TEXT);
        CREATE TABLE sample_taxon_abundances(sample_id INTEGER, taxon_id INTEGER,
            relative_abundance REAL);
        INSERT INTO studies VALUES (1,'PROJECT');
        INSERT INTO diseases VALUES (1,'Disease'),(2,'Health');
    """)
    for i, name in enumerate(["a", "b", "c", "species_a"], 1):
        c.execute("INSERT INTO taxa VALUES (?,?,?,?)", (offset+i, name, "",
                  "species" if i == 4 else "genus"))
    # Two runs from the same physical sample; choose RUN_A, never by local ID.
    for i, run in enumerate(["RUN_B", "RUN_A"], 1):
        c.execute("INSERT INTO samples VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                  (offset+i, run, 1, 1, "SAMPLE", None, None, None, None, "feces", "1"))
        for j in range(1, 5):
            c.execute("INSERT INTO sample_taxon_abundances VALUES (?,?,?)",
                      (offset+i, offset+j, 1/3 if j < 4 else 1.0))
    return c


class StoreIdentityTests(unittest.TestCase):
    def test_stable_taxon_keys_and_rank_selection_ignore_local_ids(self):
        frames = []
        for offset in [0, 1000]:
            c = fixture(offset)
            try:
                meta, raw, selection = load_model_cohort(c, "Disease", sqlite=True)
            finally:
                c.close()
            self.assertEqual(meta.run_accession.tolist(), ["RUN_A"])
            self.assertTrue(meta.age_years.isna().all())
            self.assertEqual(raw.shape, (1, 3))
            self.assertEqual(selection["profile_eligible_rows_before_dedup"], 2)
            frames.append(raw.reset_index(drop=True))
        pd.testing.assert_frame_equal(*frames)


if __name__ == "__main__":
    unittest.main()
