"""Tests for the richness confound diagnostic.

The property under test is that pooled and within-study AUC disagree exactly
when the depth difference comes from pooling. That disagreement is the whole
output of gutdb/diagnostics.py, so it is what has to be pinned down: a version
that computed both the same way would have reported Neuromyelitis Optica's 0.982
as a within-study effect and called it biology.
"""
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.diagnostics import (body_site_confound, rank_auc,
                               richness_confound)


def cohort(rows):
    return pd.DataFrame(rows, columns=["phenotype", "project_accession",
                                       "n_genera_profiled"])


class RankAucTests(unittest.TestCase):
    def test_perfect_separation_is_one(self):
        self.assertEqual(rank_auc(pd.Series([1, 1, 0, 0]), pd.Series([4, 3, 2, 1])), 1.0)

    def test_reversed_separation_is_zero(self):
        self.assertEqual(rank_auc(pd.Series([1, 1, 0, 0]), pd.Series([1, 2, 3, 4])), 0.0)

    def test_all_ties_are_chance(self):
        """An integer predictor is full of ties and must not be scored as signal."""
        self.assertEqual(rank_auc(pd.Series([1, 1, 0, 0]), pd.Series([5, 5, 5, 5])), 0.5)

    def test_one_empty_arm_scores_nothing_rather_than_guessing(self):
        self.assertIsNone(rank_auc(pd.Series([1, 1, 1]), pd.Series([3, 2, 1])))

    def test_index_is_not_assumed_to_be_a_range(self):
        labels = pd.Series([1, 1, 0, 0], index=[11, 22, 33, 44])
        self.assertEqual(rank_auc(labels, pd.Series([4, 3, 2, 1], index=labels.index)),
                         1.0)


class RichnessConfoundTests(unittest.TestCase):
    def test_pooling_artefact_is_visible_as_a_split_between_the_two_numbers(self):
        """The Neuromyelitis Optica shape: a thin case study, rich pooled controls.

        Every sample in study THIN is thin, cases and controls alike, so inside
        that study richness says nothing. Pooled against rich controls from
        elsewhere it separates the arms completely.
        """
        rows = [("D", "THIN", r) for r in (2, 3, 2, 3, 4)]
        rows += [("Health", "THIN", r) for r in (2, 3, 4, 2, 3)]
        rows += [("Health", f"RICH{i}", 40 + i) for i in range(10)]
        out = richness_confound(cohort(rows))
        self.assertGreater(out["richness_auc_pooled"], 0.8)
        self.assertAlmostEqual(out["richness_auc_within_study"], 0.5, delta=0.1)
        self.assertEqual(out["studies_with_both_arms"], 1)

    def test_a_real_within_study_effect_survives(self):
        """The Caliciviridae shape: the gap holds up inside the study itself."""
        rows = [("D", "ONE", r) for r in (10, 12, 14, 11, 13)]
        rows += [("Health", "ONE", r) for r in (28, 30, 32, 29, 31)]
        out = richness_confound(cohort(rows))
        self.assertGreater(out["richness_auc_pooled"], 0.9)
        self.assertGreater(out["richness_auc_within_study"], 0.9)

    def test_depth_ratio_flags_a_study_that_split_its_own_arms(self):
        """The PRJNA769052 shape, which holding the study out does not remove."""
        rows = [("D", "ONE", r) for r in (1, 2, 2, 3, 2)]
        rows += [("Health", "ONE", r) for r in (40, 43, 45, 42, 44)]
        out = richness_confound(cohort(rows))
        ratio = out["per_study_depth"][0]["depth_ratio"]
        self.assertLess(ratio, 0.1)

    def test_balanced_depth_reports_chance_and_no_flag(self):
        rows = [("D", "ONE", r) for r in (40, 42, 44, 41, 43)]
        rows += [("Health", "ONE", r) for r in (41, 43, 40, 44, 42)]
        out = richness_confound(cohort(rows))
        self.assertAlmostEqual(out["richness_auc_pooled"], 0.5, delta=0.15)
        self.assertAlmostEqual(out["per_study_depth"][0]["depth_ratio"], 1.0,
                               delta=0.1)

    def test_an_arm_too_small_to_score_is_left_out_of_the_within_figure(self):
        """A two-sample arm cannot support an AUC, and must not contribute one."""
        rows = [("D", "SMALL", 2), ("Health", "SMALL", 40)]
        rows += [("D", "BIG", r) for r in (5, 6, 7, 8, 9)]
        rows += [("Health", "BIG", r) for r in (41, 42, 43, 44, 45)]
        out = richness_confound(cohort(rows))
        self.assertEqual(out["studies_with_both_arms"], 1)
        self.assertEqual([d["project_accession"] for d in out["per_study_depth"]],
                         ["BIG"])

    def test_direction_is_fixed_so_the_two_numbers_are_comparable(self):
        """Fewer taxa must mean 'case' in both figures, or they cannot be compared."""
        rows = [("D", "ONE", r) for r in (2, 3, 2, 3, 2)]
        rows += [("Health", "ONE", r) for r in (40, 41, 42, 43, 44)]
        out = richness_confound(cohort(rows))
        self.assertEqual(out["richness_auc_pooled"], 1.0)
        self.assertEqual(out["richness_auc_within_study"], 1.0)
        self.assertLess(out["median_richness_case"], out["median_richness_control"])


def sites(rows):
    return pd.DataFrame(rows, columns=["phenotype", "body_site"])


class BodySiteConfoundTests(unittest.TestCase):
    def test_unmatched_site_is_flagged_and_scores_one(self):
        """The MGYS00001188 shape: biopsies as cases, stool as controls."""
        out = body_site_confound(sites(
            [("D", "gastric mucosa")] * 3 + [("Health", "feces")] * 4))
        self.assertEqual(out["sites_in_one_arm_only"], ["feces", "gastric mucosa"])
        self.assertEqual(out["shared_sites"], [])
        self.assertEqual(out["site_auc"], 1.0)
        self.assertEqual(out["n_case_at_unmatched_site"], 3)

    def test_a_matched_cohort_raises_nothing(self):
        out = body_site_confound(sites(
            [("D", "feces")] * 3 + [("Health", "feces")] * 3))
        self.assertEqual(out["sites_in_one_arm_only"], [])
        self.assertEqual(out["shared_sites"], ["feces"])
        self.assertIsNone(out["site_auc"])

    def test_a_site_present_in_both_arms_is_not_counted_as_exclusive(self):
        """Sites in both arms can be adjusted for; only one-sided ones cannot."""
        out = body_site_confound(sites(
            [("D", "feces")] * 3 + [("D", "gastric mucosa")] * 2
            + [("Health", "feces")] * 3 + [("Health", "gastric mucosa")] * 2))
        self.assertEqual(out["sites_in_one_arm_only"], [])
        self.assertEqual(out["shared_sites"], ["feces", "gastric mucosa"])
        self.assertEqual(out["n_case_at_unmatched_site"], 0)

    def test_a_partially_unmatched_cohort_scores_below_one(self):
        """Half the cases at a matched site: body site cannot separate alone."""
        out = body_site_confound(sites(
            [("D", "gastric mucosa")] * 3 + [("D", "feces")] * 3
            + [("Health", "feces")] * 6))
        self.assertEqual(out["sites_in_one_arm_only"], ["gastric mucosa"])
        self.assertLess(out["site_auc"], 1.0)
        self.assertGreater(out["site_auc"], 0.5)
        self.assertEqual(out["n_case_at_unmatched_site"], 3)


if __name__ == "__main__":
    unittest.main()
