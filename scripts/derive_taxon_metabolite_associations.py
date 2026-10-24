#!/usr/bin/env python3
"""Associate taxon abundance with a measured metabolite level, within each study.

Succinate by default. The computation is deliberately per-study and rank-based,
because the levels are not comparable across studies: for succinate the study
means span ten orders of magnitude (0.06 to 3.0e+08) and two cohorts are stored
centred, so they carry negative values. Pooling them would give a confident
number about nothing. Spearman within a study is invariant to that -- any
monotone rescaling leaves the ranks alone -- so the scale difference stops
mattering once nothing crosses a study boundary.

One row per (taxon, metabolite, study). The combination across studies is a
view, not a stored number, so the per-study evidence stays inspectable and a
single cohort driving a result is visible rather than averaged away.

The abundances are CLR transformed first, which is not a refinement but the
difference between a result and an artifact. Relative abundances in one sample
sum to 1, so when a few dominant genera rise with succinate every other part
must fall, and correlating the raw proportions manufactures negative
associations. Measured on this data: raw proportions give 463 negative against
126 positive below FDR 0.05, 79% negative. The same test on CLR values gives
268 negative against 319 positive, 46%. The skew was the simplex, not biology.

CLR also fixes what absence means here. Zero-filling a raw vector treats "not
detected" as a number on the same scale as a measurement; the CLR path routes
absences through multiplicative_replacement instead, which keeps the ratios
between the parts that were observed. Support is therefore reported separately
as n_nonzero -- how many samples the taxon was genuinely present in -- because
after replacement nothing is zero any more and n_samples alone would hide a
correlation resting on a handful of real observations.

Preview by default; --apply writes.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.compositional import clr_matrix
from gutdb.config import Settings
from gutdb.db import connect

# A correlation over a handful of non-zero points is noise with a p-value.
MIN_SAMPLES = 20
MIN_NONZERO = 10
FDR = 0.05


def metabolite_by_study(connection, metabolite_id: int):
    """Per study: the metabolite level by sample, and the annotation confidence."""
    cursor = connection.cursor()
    try:
        cursor.execute(
            """SELECT s.study_id, m.sample_id, m.level, m.high_confidence_annotation
               FROM sample_metabolite_levels m
               JOIN samples s ON s.id = m.sample_id
               WHERE m.metabolite_id = %s
                 AND EXISTS (SELECT 1 FROM sample_taxon_abundances a
                             WHERE a.sample_id = m.sample_id)""",
            (metabolite_id,))
        levels: dict[int, dict[int, float]] = {}
        confidence: dict[int, bool] = {}
        for study_id, sample_id, level, high in cursor.fetchall():
            levels.setdefault(study_id, {})[sample_id] = float(level)
            confidence[study_id] = confidence.get(study_id, True) and bool(high)
    finally:
        cursor.close()
    return levels, confidence


def prevalence(connection, sample_ids: list[int]) -> dict[int, int]:
    """Samples each taxon was actually present in, before any replacement."""
    cursor = connection.cursor()
    try:
        cursor.execute(
            "SELECT taxon_id, COUNT(*) FROM sample_taxon_abundances "
            "WHERE relative_abundance > 0 AND sample_id IN (%s) GROUP BY taxon_id"
            % ",".join(["%s"] * len(sample_ids)), sample_ids)
        return dict(cursor.fetchall())
    finally:
        cursor.close()


def correlate(connection, levels, confidence, rank: str) -> list[tuple]:
    rows = []
    for study_id, by_sample in levels.items():
        if len(by_sample) < MIN_SAMPLES:
            continue
        # clr_matrix drops its own thin rows, so the level vector is taken in
        # the order of the matrix it returns rather than of the request.
        matrix, _ = clr_matrix(connection, rank, sorted(by_sample))
        if matrix.empty or len(matrix) < MIN_SAMPLES:
            continue
        metabolite = np.array([by_sample[s] for s in matrix.index])
        present = prevalence(connection, list(matrix.index))
        found = []
        for taxon_id in matrix.columns:
            nonzero = int(present.get(taxon_id, 0))
            if nonzero < MIN_NONZERO:
                continue
            rho, p = stats.spearmanr(matrix[taxon_id].to_numpy(), metabolite)
            if not np.isfinite(rho) or not np.isfinite(p):
                continue
            found.append((taxon_id, float(rho), float(p), len(matrix), nonzero))
        if not found:
            continue
        # Benjamini-Hochberg within the study, over the taxa actually tested
        # there. Thousands of raw p-values across 11 cohorts would otherwise
        # hand back a page of discoveries that are the multiplicity.
        q = stats.false_discovery_control([f[2] for f in found], method="bh")
        for (taxon_id, rho, p, n, nonzero), q_value in zip(found, q):
            rows.append((taxon_id, study_id, round(rho, 4), p, float(q_value),
                         n, nonzero, confidence[study_id]))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env-file", default=".env")
    parser.add_argument("--kegg-id", default="C00042", help="default C00042, succinate")
    parser.add_argument("--rank", default="genus",
                        help="the abundance layer; the metabolomics cohorts carry only genus")
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    connection = connect(Settings.from_env(args.env_file))
    try:
        cursor = connection.cursor()
        cursor.execute("SELECT id, name FROM metabolites WHERE kegg_id = %s", (args.kegg_id,))
        found = cursor.fetchone()
        cursor.close()
        if not found:
            raise SystemExit(f"No metabolite with kegg_id {args.kegg_id}")
        metabolite_id, name = found
        print(f"{name} (metabolite_id {metabolite_id})")

        levels, confidence = metabolite_by_study(connection, metabolite_id)
        print(f"studies with both the metabolite and profiles: {len(levels)}")
        rows = correlate(connection, levels, confidence, args.rank)
        studies = len({r[1] for r in rows})
        taxa = len({r[0] for r in rows})
        passing = [r for r in rows if r[4] < FDR]
        print(f"tested {len(rows):,} (taxon, study) pairs: {taxa} taxa across {studies} studies")
        print(f"below FDR {FDR}: {len(passing):,} "
              f"({sum(1 for r in passing if r[2] > 0):,} positive, "
              f"{sum(1 for r in passing if r[2] < 0):,} negative)")
        low = {r[1] for r in rows if not r[7]}
        if low:
            print(f"{len(low)} study/studies carry a low-confidence annotation upstream; "
                  "kept and flagged, not dropped")
        if not args.apply:
            print("Preview only. Use --apply to store.")
            return
        cursor = connection.cursor()
        try:
            cursor.executemany(
                """INSERT INTO taxon_metabolite_associations
                   (taxon_id, metabolite_id, study_id, rho, p_value, q_value,
                    n_samples, n_nonzero, high_confidence_annotation)
                   VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
                   ON DUPLICATE KEY UPDATE rho=VALUES(rho), p_value=VALUES(p_value),
                     q_value=VALUES(q_value), n_samples=VALUES(n_samples),
                     n_nonzero=VALUES(n_nonzero),
                     high_confidence_annotation=VALUES(high_confidence_annotation)""",
                [(r[0], metabolite_id, r[1], r[2], r[3], r[4], r[5], r[6], r[7])
                 for r in rows])
            connection.commit()
            cursor.execute("SELECT COUNT(*) FROM taxon_metabolite_associations "
                           "WHERE metabolite_id = %s", (metabolite_id,))
            print(f"stored rows for this metabolite: {cursor.fetchone()[0]:,}")
        finally:
            cursor.close()
    finally:
        connection.close()


if __name__ == "__main__":
    main()
