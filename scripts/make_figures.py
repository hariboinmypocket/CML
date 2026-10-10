#!/usr/bin/env python3
"""Regenerate docs/gutdb_overview.png from the live database.

The previous copy of this figure was produced outside the repository and could
not be reproduced from it, which is how it came to show numbers from a stale
mirror (16,390 associations, 2,131 replicated pairs) alongside directions that
the 2026-10-10 fix to v_taxon_disease_evidence has since reversed. This script
exists so the figure is a function of the database rather than an artefact.

Both panels depend on direction, so both restrict to CASE/CONTROL contrasts:

  Panel a splits each disease's associations into enriched and depleted. 2,340
  of the 16,889 associations compare one disease against another, where
  "enriched in disease" means enriched relative to a different illness rather
  than to health. Counting those in a bar labelled "enriched in disease" is a
  category error, so they are excluded and the excluded count is stated.

  Panel b plots directional agreement against replication depth, reading
  agreement_ratio from the corrected view, which takes one vote per study from
  case/control contrasts only.

Writes docs/gutdb_overview.png. Pass --show-counts to print the numbers that
appear in the titles, for pasting into a caption.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.ticker  # noqa: F401
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect

OUT = Path("docs/gutdb_overview.png")
DEPLETED = "#2b6cb0"
ENRICHED = "#c0532a"
GREY = "#8a8a8a"

PANEL_A_SQL = """
SELECT d.name AS disease_name,
       SUM(a.direction = 'enriched') AS n_enriched,
       SUM(a.direction = 'depleted') AS n_depleted,
       COUNT(DISTINCT c.study_id) AS n_studies
FROM taxon_disease_associations a
JOIN diseases d ON d.id = a.disease_id
JOIN phenotype_comparisons c ON c.id = a.comparison_id
LEFT JOIN diseases pa ON pa.id = c.phenotype_a_id
LEFT JOIN diseases pb ON pb.id = c.phenotype_b_id
WHERE COALESCE(pa.name = 'Health' OR pb.name = 'Health', 0) = 1
GROUP BY d.id, d.name
ORDER BY (SUM(a.direction = 'enriched') + SUM(a.direction = 'depleted')) DESC
LIMIT 15
"""

PANEL_B_SQL = """
SELECT n_studies_case_control AS n_studies,
       agreement_ratio,
       n_studies_enriched,
       n_studies_depleted
FROM v_taxon_disease_evidence
WHERE agreement_ratio IS NOT NULL
  AND (n_studies_enriched + n_studies_depleted) > 1
"""

TOTALS_SQL = """
SELECT
  (SELECT COUNT(*) FROM taxon_disease_associations) AS n_assoc,
  (SELECT COUNT(*) FROM taxon_disease_associations a
     JOIN phenotype_comparisons c ON c.id = a.comparison_id
     LEFT JOIN diseases pa ON pa.id = c.phenotype_a_id
     LEFT JOIN diseases pb ON pb.id = c.phenotype_b_id
   WHERE COALESCE(pa.name = 'Health' OR pb.name = 'Health', 0) = 1) AS n_case_control
"""


def main() -> int:
    conn = connect(Settings.from_env(".env"))
    a = pd.read_sql_query(PANEL_A_SQL, conn)
    b = pd.read_sql_query(PANEL_B_SQL, conn)
    totals = pd.read_sql_query(TOTALS_SQL, conn).iloc[0]
    conn.close()

    a["total"] = a["n_depleted"] + a["n_enriched"]
    top_share = a["total"].head(3).sum() / totals["n_case_control"]

    n_pairs = len(b)
    unanimous = int((b["agreement_ratio"] >= 0.999).sum())
    conflicted = int(((b["n_studies_enriched"] > 0) & (b["n_studies_depleted"] > 0)).sum())
    near_tie = int((b["agreement_ratio"] <= 0.6).sum())

    plt.rcParams.update({
        "font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
        "axes.titlesize": 12, "figure.dpi": 150,
    })
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16.4, 6.5),
                                   gridspec_kw={"width_ratios": [1.05, 1]})

    order = a.iloc[::-1]
    y = np.arange(len(order))
    ax1.barh(y, order["n_depleted"], color=DEPLETED, label="depleted in disease")
    ax1.barh(y, order["n_enriched"], left=order["n_depleted"], color=ENRICHED,
             label="enriched in disease")
    ax1.set_yticks(y)
    ax1.set_yticklabels(order["disease_name"])
    for yi, (tot, ns) in enumerate(zip(order["total"], order["n_studies"])):
        ax1.text(tot + totals["n_case_control"] * 0.004, yi,
                 f"{ns} stud{'y' if ns == 1 else 'ies'}",
                 va="center", fontsize=9, color=GREY)
    ax1.set_xlabel("curated case/control marker associations")
    top3 = list(a["disease_name"].head(3))
    ax1.set_title(
        f"{top3[0]}, {top3[1]} and {top3[2]}\n"
        f"carry {top_share:.0%} of the {int(totals['n_case_control']):,} "
        f"case/control associations",
        loc="left", pad=14)
    ax1.legend(loc="lower right", frameon=False)
    ax1.margins(x=0.12)

    rng = np.random.default_rng(0)
    jitter = rng.normal(0, 0.055, len(b))
    conflict_mask = (b["n_studies_enriched"] > 0) & (b["n_studies_depleted"] > 0)
    ax2.axhline(1.0, ls="--", lw=1, color=GREY, zorder=1)
    ax2.scatter(b.loc[~conflict_mask, "n_studies"] + jitter[~conflict_mask.to_numpy()],
                b.loc[~conflict_mask, "agreement_ratio"],
                s=22, alpha=0.45, color=DEPLETED, linewidths=0, zorder=2)
    ax2.scatter(b.loc[conflict_mask, "n_studies"] + jitter[conflict_mask.to_numpy()],
                b.loc[conflict_mask, "agreement_ratio"],
                s=22, alpha=0.45, color=ENRICHED, linewidths=0, zorder=3)
    ax2.set_xscale("log")
    ax2.set_xticks([2, 3, 4, 5, 7, 10, 15])
    ax2.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax2.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax2.set_xlim(1.75, max(18, b["n_studies"].max() * 1.15))
    ax2.set_xlabel("independent case/control studies reporting the pair")
    ax2.set_ylabel("directional agreement among studies")
    ax2.text(0.015, 1.0, "  unanimous", transform=ax2.get_yaxis_transform(),
             color=GREY, fontsize=9, va="bottom", ha="left")
    ax2.text(0.98, 0.10, f"studies disagree\nin direction ({conflicted:,} pairs)",
             transform=ax2.transAxes, color=ENRICHED, fontsize=10,
             ha="right", va="bottom")
    ax2.set_title(
        f"Of {n_pairs:,} replicated pairs, {unanimous/n_pairs:.0%} are unanimous;\n"
        f"{near_tie:,} sit at or below 60% agreement",
        loc="left", pad=14)

    for ax, letter in ((ax1, "a"), (ax2, "b")):
        ax.text(-0.09 if ax is ax1 else -0.07, 1.06, letter, transform=ax.transAxes,
                fontsize=15, fontweight="bold", va="top")

    excluded = int(totals["n_assoc"]) - int(totals["n_case_control"])
    fig.text(0.5, -0.01,
             f"Case/control contrasts only: {excluded:,} of {int(totals['n_assoc']):,} "
             f"associations compare one disease against another and are excluded, "
             f"since “enriched in disease” there means enriched relative to a "
             f"different illness.",
             ha="center", fontsize=9, color=GREY)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUT, bbox_inches="tight", facecolor="white")
    print(f"wrote {OUT}")

    if "--show-counts" in sys.argv:
        print(f"\n  associations total        {int(totals['n_assoc']):,}")
        print(f"  case/control only         {int(totals['n_case_control']):,}")
        print(f"  excluded (disease v dis)  {excluded:,}")
        print(f"  replicated pairs plotted  {n_pairs:,}")
        print(f"  unanimous                 {unanimous:,} ({unanimous/n_pairs:.1%})")
        print(f"  with a directional clash  {conflicted:,}")
        print(f"  agreement <= 0.60         {near_tie:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
