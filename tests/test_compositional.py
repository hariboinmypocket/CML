"""Tests for the CLR pipeline.

These assert mathematical properties rather than remembered outputs: a
compositional transform is only trustworthy if it has the invariances it claims,
and those are cheap to state exactly.
"""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.compositional import clr, multiplicative_replacement


def frame(rows):
    return pd.DataFrame(rows, index=range(len(rows)), columns=range(len(rows[0])))


class MultiplicativeReplacementTests(unittest.TestCase):
    def test_rows_still_sum_to_one(self):
        out = multiplicative_replacement(frame([[0.5, 0.5, 0.0, 0.0]]))
        self.assertAlmostEqual(float(out.iloc[0].sum()), 1.0, places=12)

    def test_ratios_between_observed_parts_are_preserved(self):
        """The reason not to just add a pseudocount to everything.

        A composition carries only ratios, so a zero-replacement that changes the
        ratios of the parts that WERE measured has destroyed the data it was
        meant to protect.
        """
        out = multiplicative_replacement(frame([[0.6, 0.3, 0.1, 0.0, 0.0]]))
        self.assertAlmostEqual(out.iloc[0, 0] / out.iloc[0, 1], 0.6 / 0.3, places=12)
        self.assertAlmostEqual(out.iloc[0, 1] / out.iloc[0, 2], 0.3 / 0.1, places=12)

    def test_naive_pseudocount_would_not_preserve_them(self):
        naive = np.array([0.6, 0.3, 0.1, 0.0, 0.0]) + 1e-3
        naive = naive / naive.sum()
        self.assertNotAlmostEqual(naive[0] / naive[1], 0.6 / 0.3, places=4)

    def test_zeros_become_positive(self):
        out = multiplicative_replacement(frame([[0.9, 0.1, 0.0]]))
        self.assertGreater(out.iloc[0, 2], 0.0)

    def test_delta_adapts_to_the_sample(self):
        """A shallow sample's detection limit is not a deep sample's."""
        out = multiplicative_replacement(frame([[0.5, 0.5, 0.0], [0.999, 0.001, 0.0]]))
        self.assertGreater(out.iloc[0, 2], out.iloc[1, 2])

    def test_row_of_all_zeros_is_nan_not_invented(self):
        out = multiplicative_replacement(frame([[0.0, 0.0, 0.0]]))
        self.assertTrue(out.iloc[0].isna().all())


class ClrTests(unittest.TestCase):
    def test_each_row_sums_to_zero(self):
        """The defining property: CLR removes the constant-sum constraint."""
        out = clr(frame([[0.5, 0.3, 0.2], [0.1, 0.1, 0.8]]))
        for i in range(len(out)):
            self.assertAlmostEqual(float(out.iloc[i].sum()), 0.0, places=12)

    def test_invariant_to_rescaling_a_row(self):
        """Counts or proportions give the same answer -- that is the point."""
        a = clr(frame([[0.5, 0.3, 0.2]]))
        b = clr(frame([[500.0, 300.0, 200.0]]))
        np.testing.assert_allclose(a.to_numpy(), b.to_numpy(), atol=1e-12)

    def test_rejects_zeros_rather_than_silently_producing_inf(self):
        with self.assertRaises(ValueError):
            clr(frame([[0.5, 0.5, 0.0]]))

    def test_not_invariant_to_the_part_set(self):
        """Why CLR is not stored as a column.

        Adding one part changes every other value, so any stored CLR is wrong for
        a modeller who filters the feature set -- which is the first thing a
        modeller does.
        """
        three = clr(frame([[0.5, 0.3, 0.2]])).to_numpy()[0]
        four = clr(frame([[0.5, 0.3, 0.2, 0.001]])).to_numpy()[0][:3]
        self.assertFalse(np.allclose(three, four, atol=0.01))

    def test_pipeline_end_to_end_on_a_sparse_row(self):
        raw = frame([[0.7, 0.2, 0.1, 0.0, 0.0, 0.0]])
        out = clr(multiplicative_replacement(raw))
        self.assertAlmostEqual(float(out.iloc[0].sum()), 0.0, places=12)
        # The observed parts keep their order; the replaced zeros sit below them.
        self.assertGreater(out.iloc[0, 0], out.iloc[0, 1])
        self.assertGreater(out.iloc[0, 2], out.iloc[0, 3])


class DegenerateRowTests(unittest.TestCase):
    def test_clr_of_a_single_part_is_zero(self):
        """Why a one-taxon profile has to be excluded rather than transformed.

        With one part, that part IS the geometric mean, so its CLR is zero no
        matter what the abundance was. 269 samples in this database hold exactly
        one genus.
        """
        out = clr(frame([[1.0]]))
        self.assertAlmostEqual(float(out.iloc[0, 0]), 0.0, places=12)

    def test_two_parts_carry_information(self):
        out = clr(frame([[0.9, 0.1]]))
        self.assertGreater(out.iloc[0, 0], out.iloc[0, 1])
        self.assertAlmostEqual(float(out.iloc[0].sum()), 0.0, places=12)


if __name__ == "__main__":
    unittest.main()
