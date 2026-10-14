"""Compositional transforms for the abundance matrix.

Relative abundances are compositional: they carry only ratios, they sum to 1
within a rank, and the usual statistics do not apply to them. Two taxa cannot
both rise, so an increase in one forces an apparent decrease in the others, and a
model trained on raw proportions learns that artefact. The centred log-ratio
transform removes it by expressing each part against the sample's own geometric
mean.

NOTHING HERE IS STORED IN THE DATABASE, and that is a deliberate decision rather
than work left undone.

  CLR is defined relative to the set of parts included. Adding or removing one
  taxon changes the geometric mean and therefore every other value in the sample.
  For a toy composition of three parts the CLR is [0.476, -0.035, -0.441]; add a
  fourth part at 0.1% and the same three become [1.910, 1.400, 0.994]. A stored
  column would be silently wrong the moment anyone filtered the feature set,
  which is the first thing any modeller does.

  The matrix is sparse, and a transform over all of it is mostly pseudocount. A
  genus profile holds about 50 of 2,076 taxa, so 99% of the log-sum would come
  from replaced zeros: the geometric mean becomes a function of the replacement
  constant rather than of the sample. Filtering to prevalent taxa first is not
  tidying, it is what makes the number mean anything.

So the order below is load-bearing: select samples, filter taxa by prevalence,
RENORMALIZE, replace zeros, then transform. Renormalizing after filtering is what
keeps the result a composition; skipping it leaves rows summing to less than 1
and the geometric mean quietly wrong.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

RANKS = ("genus", "species")


def prevalent_taxa(connection, rank: str, sample_ids: list[int],
                   min_prevalence: float = 0.10) -> list[int]:
    """Taxon ids detected in at least `min_prevalence` of the given samples.

    Done in SQL because the point is to avoid materialising the full matrix:
    32,327 samples by 2,076 genera is 67M cells, and almost all of it is zero.
    """
    if not sample_ids:
        return []
    placeholders = ",".join(["%s"] * len(sample_ids))
    cur = connection.cursor()
    cur.execute(
        f"""SELECT a.taxon_id
            FROM sample_taxon_abundances a
            JOIN taxa t ON t.id = a.taxon_id
            WHERE t.taxonomic_rank = %s AND a.sample_id IN ({placeholders})
            GROUP BY a.taxon_id
            HAVING COUNT(DISTINCT a.sample_id) >= %s""",
        (rank, *sample_ids, max(1, int(round(min_prevalence * len(sample_ids))))))
    taxa = [r[0] for r in cur.fetchall()]
    cur.close()
    return taxa


def abundance_frame(connection, rank: str, sample_ids: list[int],
                    taxon_ids: list[int]) -> pd.DataFrame:
    """Dense samples-by-taxa frame of relative abundances, absent cells zero.

    Absence is the absence of a row in sample_taxon_abundances, so the zeros here
    are produced by the reindex rather than read from the database.
    """
    if not sample_ids or not taxon_ids:
        return pd.DataFrame(index=pd.Index(sample_ids, name="sample_id"),
                            columns=pd.Index(taxon_ids, name="taxon_id"), dtype=float)
    sp = ",".join(["%s"] * len(sample_ids))
    tp = ",".join(["%s"] * len(taxon_ids))
    frame = pd.read_sql_query(
        f"""SELECT sample_id, taxon_id, relative_abundance
            FROM sample_taxon_abundances
            WHERE sample_id IN ({sp}) AND taxon_id IN ({tp})""",
        connection, params=(*sample_ids, *taxon_ids))
    wide = frame.pivot(index="sample_id", columns="taxon_id", values="relative_abundance")
    return wide.reindex(index=sample_ids, columns=taxon_ids).fillna(0.0)


def multiplicative_replacement(frame: pd.DataFrame, delta: float | None = None
                               ) -> pd.DataFrame:
    """Replace zeros without distorting the ratios between observed parts.

    Adding a constant to every cell is the common shortcut and it is wrong: it
    changes the ratios of the parts that were actually measured, which are the
    only information a composition carries. The multiplicative replacement of
    Martin-Fernandez et al. (2003) puts `delta` in the zeros and scales the
    observed parts down by exactly the mass it introduced, leaving their ratios
    untouched and the row summing to 1.

    `delta` defaults per row to half that sample's smallest observed value, so a
    deeply sequenced sample gets a smaller replacement than a shallow one. That
    is the right behaviour: the zero means "below this sample's detection limit",
    and the limit is not the same for every sample.
    """
    values = frame.to_numpy(dtype=float, copy=True)
    out = np.empty_like(values)
    for i, row in enumerate(values):
        observed = row > 0
        if not observed.any():
            out[i] = np.nan
            continue
        d = delta if delta is not None else row[observed].min() / 2.0
        n_zero = int((~observed).sum())
        total = row.sum()
        if total <= 0:
            out[i] = np.nan
            continue
        shrink = 1.0 - (n_zero * d) / total
        if shrink <= 0:
            # More zeros than the composition can absorb at this delta; fall back
            # to a delta that leaves the observed parts half the mass.
            d = total / (2.0 * n_zero)
            shrink = 0.5
        out[i] = np.where(observed, row * shrink, d)
    return pd.DataFrame(out, index=frame.index, columns=frame.columns)


def clr(frame: pd.DataFrame) -> pd.DataFrame:
    """Centred log-ratio: ln(x_i) minus the row's mean log.

    Rows must be strictly positive, so run multiplicative_replacement first.
    Each transformed row sums to zero by construction, which is the property to
    assert in a test: it is what makes the result free of the constant-sum
    constraint the raw proportions carry.
    """
    values = frame.to_numpy(dtype=float)
    if np.nanmin(values) <= 0:
        raise ValueError(
            "CLR needs strictly positive values; call multiplicative_replacement first")
    logs = np.log(values)
    return pd.DataFrame(logs - logs.mean(axis=1, keepdims=True),
                        index=frame.index, columns=frame.columns)


def clr_matrix(connection, rank: str, sample_ids: list[int],
               min_prevalence: float = 0.10, delta: float | None = None,
               min_taxa_present: int = 2) -> tuple[pd.DataFrame, dict]:
    """Samples by taxa, CLR transformed, in the order that keeps it meaningful.

    Returns the matrix and a dict of what was done to produce it, because a CLR
    value is uninterpretable without knowing which parts it was computed over.

    `min_taxa_present` is a mathematical floor, not a quality judgement. CLR of a
    single part is identically zero -- the part IS the geometric mean -- so a row
    with one surviving taxon carries no information and its pseudocounts would
    sit at whatever delta the replacement chose. 269 samples in this database
    hold exactly one genus, faithfully: GMrepo reports one row for them and their
    mean unclassified_fraction is 0.64 against 0.20 database-wide, so
    classification largely failed rather than the loader dropping anything. 57 of
    them reach the deduplicated ML cohort, where ml_ready asks only for
    n_taxa_profiled > 0.
    Raising this above 2 is a modelling choice and is left to the caller; for
    scale, the median sample holds 62 genera and the 5th percentile holds 9.
    """
    if rank not in RANKS:
        raise ValueError(f"rank must be one of {RANKS}, got {rank!r}")
    taxa = prevalent_taxa(connection, rank, sample_ids, min_prevalence)
    raw = abundance_frame(connection, rank, sample_ids, taxa)

    kept_mass = raw.sum(axis=1)
    # Renormalize: filtering removed parts, so the rows no longer sum to 1 and
    # the geometric mean would be computed over a sub-composition pretending to
    # be a whole one.
    renormalized = raw.div(kept_mass.replace(0, np.nan), axis=0)
    usable = renormalized.dropna(how="all")

    present = (usable > 0).sum(axis=1)
    too_thin = present[present < min_taxa_present]
    usable = usable.drop(index=too_thin.index)

    replaced = multiplicative_replacement(usable, delta=delta)
    transformed = clr(replaced.dropna(how="any"))

    provenance = {
        "rank": rank,
        "samples_requested": len(sample_ids),
        "samples_returned": len(transformed),
        "taxa_kept": len(taxa),
        "min_prevalence": min_prevalence,
        "delta": "per-sample: half the smallest observed value" if delta is None else delta,
        "median_mass_retained_by_filter": round(float(kept_mass.median()), 4),
        "min_taxa_present": min_taxa_present,
        "samples_dropped_too_few_taxa": int(len(too_thin)),
    }
    return transformed, provenance
