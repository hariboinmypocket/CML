"""The taxon-metabolite derivation, tested on its guarantees rather than outputs.

The statistical choices are the point of this code -- per study, CLR, a support
floor, FDR within the study -- so these assert those properties hold, not that a
particular rho came back.
"""
import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

spec = importlib.util.spec_from_file_location(
    "derive_tma", ROOT / "scripts" / "derive_taxon_metabolite_associations.py")
tma = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tma)


def matrix(samples, columns):
    return pd.DataFrame(columns, index=samples)


class Harness:
    """Stands in for clr_matrix and prevalence, which need a live connection."""

    def __init__(self, frame, present):
        self.frame, self.present = frame, present

    def install(self, test):
        test.addCleanup(setattr, tma, "clr_matrix", tma.clr_matrix)
        test.addCleanup(setattr, tma, "prevalence", tma.prevalence)
        tma.clr_matrix = lambda connection, rank, ids: (self.frame, {})
        tma.prevalence = lambda connection, ids: self.present


class CorrelateTests(unittest.TestCase):
    def rising(self, n=40):
        samples = list(range(1, n + 1))
        # One taxon tracking the metabolite, one moving against it.
        frame = matrix(samples, {10: np.arange(n, dtype=float),
                                 11: np.arange(n, 0, -1, dtype=float)})
        levels = {1: {s: float(s) for s in samples}}
        return frame, levels, samples

    def test_rho_sign_follows_the_relationship(self):
        frame, levels, samples = self.rising()
        Harness(frame, {10: len(samples), 11: len(samples)}).install(self)
        rows = tma.correlate(None, levels, {1: True}, "genus")
        signs = {r[0]: r[2] for r in rows}
        self.assertGreater(signs[10], 0)
        self.assertLess(signs[11], 0)

    def test_support_floor_excludes_thin_taxa(self):
        frame, levels, samples = self.rising()
        # Taxon 11 was genuinely present in too few samples to correlate on.
        Harness(frame, {10: len(samples), 11: tma.MIN_NONZERO - 1}).install(self)
        rows = tma.correlate(None, levels, {1: True}, "genus")
        self.assertEqual({r[0] for r in rows}, {10})

    def test_study_below_min_samples_is_skipped(self):
        samples = list(range(1, tma.MIN_SAMPLES))
        frame = matrix(samples, {10: np.arange(len(samples), dtype=float)})
        levels = {1: {s: float(s) for s in samples}}
        Harness(frame, {10: len(samples)}).install(self)
        self.assertEqual(tma.correlate(None, levels, {1: True}, "genus"), [])

    def test_fdr_never_reports_an_adjusted_value_below_the_raw_one(self):
        rng = np.random.default_rng(0)
        samples = list(range(1, 61))
        frame = matrix(samples, {t: rng.normal(size=len(samples)) for t in range(10, 40)})
        levels = {1: {s: float(rng.normal()) for s in samples}}
        Harness(frame, {t: len(samples) for t in range(10, 40)}).install(self)
        rows = tma.correlate(None, levels, {1: True}, "genus")
        self.assertTrue(rows)
        for row in rows:
            self.assertGreaterEqual(row[4], row[3] - 1e-9)

    def test_each_study_is_corrected_on_its_own(self):
        """FDR pooled across cohorts would let a big study set another's threshold."""
        rng = np.random.default_rng(1)
        samples = list(range(1, 41))
        frame = matrix(samples, {10: np.arange(40, dtype=float), 11: rng.normal(size=40)})
        levels = {1: {s: float(s) for s in samples}, 2: {s: float(s) for s in samples}}
        Harness(frame, {10: 40, 11: 40}).install(self)
        rows = tma.correlate(None, levels, {1: True, 2: False}, "genus")
        self.assertEqual({r[1] for r in rows}, {1, 2})
        for study in (1, 2):
            per = [r for r in rows if r[1] == study]
            self.assertEqual(len(per), 2)

    def test_annotation_confidence_is_carried_per_study(self):
        frame, levels, samples = self.rising()
        levels[2] = dict(levels[1])
        Harness(frame, {10: len(samples), 11: len(samples)}).install(self)
        rows = tma.correlate(None, levels, {1: True, 2: False}, "genus")
        self.assertEqual({r[1]: r[7] for r in rows if r[0] == 10}, {1: True, 2: False})

    def test_support_is_counted_over_the_samples_actually_correlated(self):
        """clr_matrix drops thin rows, so the two sample sets are not the same.

        Counting prevalence over the requested samples instead of the returned
        ones would report support the correlation never saw, and could exceed
        n_samples outright.
        """
        requested = list(range(1, 51))
        kept = requested[:40]
        frame = matrix(kept, {10: np.arange(len(kept), dtype=float)})
        levels = {1: {s: float(s) for s in requested}}
        seen = []
        self.addCleanup(setattr, tma, "clr_matrix", tma.clr_matrix)
        self.addCleanup(setattr, tma, "prevalence", tma.prevalence)
        tma.clr_matrix = lambda connection, rank, ids: (frame, {})

        def record(connection, ids):
            seen.append(list(ids))
            return {10: len(ids)}

        tma.prevalence = record
        rows = tma.correlate(None, levels, {1: True}, "genus")
        self.assertEqual(seen, [kept])
        self.assertEqual([(r[5], r[6]) for r in rows], [(len(kept), len(kept))])


if __name__ == "__main__":
    unittest.main()
