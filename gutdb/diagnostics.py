"""Checks that ask whether a cohort's separability needs any biology.

Everything here answers one question: could a number that is not a biological
measurement have produced the result? The number in this module is profile
richness -- how many taxa a sample has a row for -- and it earns its own module
because it predicts the label in this database better than most taxa do.

Neuromyelitis Optica is the case that motivated it. Its entire case arm is
PRJNA422961, a study in which no sample reaches 10 genera, while its controls are
pooled from 104 studies sitting at a median of 45. A logistic regression on CLR
genus features scores AUC 0.995 on that cohort. Richness used by itself, a single
integer, scores 0.982. Almost nothing was learned about the microbiome.

Measured within study the same number falls to 0.497, which is chance: PRJNA422961's
own controls are just as thin as its cases, so the depth difference is an artefact
of pooling and not an effect of the disease. That contrast is the whole point of
reporting both numbers, and it is the same distinction that governs directions in
v_taxon_disease_evidence and concordance in q7.

Two cases do survive within study and must be read differently:

  PRJNA788674 (Caliciviridae) holds 31 cases at a median of 15 genera against its
  own 25 controls at 30. A collapse of gut diversity in acute norovirus is a
  documented effect, so a within-study gap is evidence, not a defect.

  PRJNA769052 (COVID-19) holds 37 cases at a median of 2 genera against its own
  50 controls at 43, on one platform in one project. No virus does that. This is
  the only confound shape here that NOTHING removes: a per-sample filter cannot
  see it, and leave-one-study-out holds out the study whole, so both arms leave
  together and the artefact never shows up as a fold that fails.
"""
from __future__ import annotations

import pandas as pd

# Below this, no sample in a study can support a genus-level model. It is the
# observed ceiling of the six degenerate GMrepo studies, not a derived threshold.
USABLE_GENERA = 15


def rank_auc(labels: pd.Series, score: pd.Series) -> float | None:
    """AUC from ranks, so it is usable on an integer predictor full of ties.

    Equal to the Mann-Whitney U statistic normalised by the arm sizes. Written out
    rather than imported so that reading a diagnostic never requires sklearn.
    """
    labels = pd.Series(labels).astype(int)
    score = pd.Series(score, index=labels.index)
    n_pos = int(labels.sum())
    n_neg = int(len(labels) - n_pos)
    if not n_pos or not n_neg:
        return None
    ranks = score.rank()
    return float((ranks[labels == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def richness_confound(meta: pd.DataFrame, richness: str = "n_genera_profiled",
                      label: str = "phenotype", control: str = "Health",
                      study: str = "project_accession",
                      min_arm: int = 5) -> dict:
    """How much of a two-arm cohort's separability is profile depth alone.

    `meta` needs one row per sample with a label column, a study column and a
    richness column. Returns the pooled AUC of richness as the sole predictor,
    the same figure computed inside each study that carries both arms, and the
    per-study depth ratio that catches an imbalance no split can remove.

    Richness is negated before scoring because the thin arm is the case arm: a
    figure above 0.5 must mean "fewer taxa implies case" in both numbers, or the
    pooled and within-study values would not be comparable.
    """
    y = (meta[label] != control).astype(int)
    pooled = rank_auc(y, -meta[richness])

    within, ratios = [], []
    for accession, grp in meta.groupby(study, dropna=False):
        cases = grp[grp[label] != control]
        controls = grp[grp[label] == control]
        if len(cases) < min_arm or len(controls) < min_arm:
            continue
        a = rank_auc((grp[label] != control).astype(int), -grp[richness])
        if a is not None:
            within.append(a)
        ctrl_median = float(controls[richness].median())
        ratios.append({
            "project_accession": accession,
            "n_case": len(cases),
            "n_control": len(controls),
            "median_richness_case": float(cases[richness].median()),
            "median_richness_control": ctrl_median,
            "depth_ratio": (round(float(cases[richness].median()) / ctrl_median, 3)
                            if ctrl_median else None),
        })

    return {
        "richness_auc_pooled": None if pooled is None else round(pooled, 3),
        "richness_auc_within_study": (round(sum(within) / len(within), 3)
                                      if within else None),
        "studies_with_both_arms": len(within),
        "median_richness_case": float(meta.loc[y == 1, richness].median()),
        "median_richness_control": float(meta.loc[y == 0, richness].median()),
        "per_study_depth": ratios,
    }


def body_site_confound(meta: pd.DataFrame, site: str = "body_site",
                       label: str = "phenotype", control: str = "Health") -> dict:
    """Whether a cohort's two arms were taken from different material.

    Richness at least measures something about the sample. Body site does not
    even do that: if the cases are stomach biopsies and the controls are stool,
    the contrast answers "which organ is this" and the answer happens to
    correlate perfectly with the label.

    MGYS00001188 is the case in this database. It contributes 120 gastric mucosa
    samples and NO control arm of its own, so Gastritis, Metaplasia and Stomach
    Neoplasms are each scored against pooled stool controls. All three reach an
    AUC of exactly 1.000 on CLR genus features, and body site alone -- one
    indicator variable -- also reaches exactly 1.000. The communities share no
    abundant genus: Haematomicrobium, Helicobacter, Tsukamurella, Gordonia and
    Thermus against Bacteroides, Phocaeicola, Faecalibacterium, Prevotella and
    Bifidobacterium.

    `sites_in_one_arm_only` is the number that decides it. A site present in both
    arms can be adjusted for; a site that appears only among the cases cannot,
    because there is no control at that site to compare against.
    """
    y = (meta[label] != control).astype(int)
    case_sites = set(meta.loc[y == 1, site].dropna())
    control_sites = set(meta.loc[y == 0, site].dropna())
    exclusive = (case_sites ^ control_sites)

    # Scored as "does this sample sit at a site the cases occupy", which is the
    # strongest single indicator a model could build out of body site alone.
    indicator = meta[site].isin(case_sites - control_sites).astype(int)
    return {
        "case_sites": sorted(case_sites),
        "control_sites": sorted(control_sites),
        "shared_sites": sorted(case_sites & control_sites),
        "sites_in_one_arm_only": sorted(exclusive),
        "site_auc": (None if not (case_sites - control_sites)
                     else round(rank_auc(y, indicator) or 0.5, 3)),
        "n_case_at_unmatched_site": int(
            ((y == 1) & meta[site].isin(case_sites - control_sites)).sum()),
    }
