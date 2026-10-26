"""Which traits a genus constrains, and what it means to dissent from one.

`gram_stain_against_phylum` showed that a trait contradicting its neighbours is
almost always a bad source value rather than a discovery, and that the cost of
not noticing is paid later: MiMeDB called Bacteroides vulgatus gram-positive,
the merge backfill moved that onto Phocaeicola vulgatus, and it sat on 10,263
abundance rows as the only positive among 109 Bacteroidetes.

That check covered one trait at one rank. This module generalizes it, and the
rank came out of measuring rather than out of the analogy. Share of groups with
at least 5 typed rows that reach 90% agreement, measured 2026-10-25 over 5,838
taxa:

    trait              genus  family  order  class  phylum
    ph_preference       100%     89%    85%    75%     60%
    sporulation          94%     85%    74%    72%     80%
    gram_stain           93%     87%    90%    80%     69%
    flagella_presence    92%     77%    75%    62%     64%
    mobility             87%     68%    58%    54%     50%
    temperature_range    71%     74%    74%    62%     50%
    oxygen_requirement   49%     52%    59%    57%     38%

Genus is the level that works, and two traits do not work at any level:

  * `oxygen_requirement` peaks at 59%. Oxygen tolerance genuinely varies inside
    a genus -- that variation is the aerotolerance signal this database exists
    to study -- so a check would flag biology, not defects.
  * `temperature_range` peaks at 74%, which is not enough to call a dissenter
    suspicious.

Leaving those two unchecked is the finding, not a gap. A check nobody can trust
gets ignored, and an ignored check is worse than none because it still has to be
gated on.

Two uses, one definition, so they cannot drift: validate counts dissenters, and
the merge backfills refuse to *create* one by moving a value onto a row whose
genus contradicts it.
"""
from __future__ import annotations

from typing import Any, Iterable

# Traits whose value is near-constant within a genus. See the table above.
GENUS_CONSTRAINED_TRAITS = (
    "gram_stain",
    "sporulation",
    "ph_preference",
    "mobility",
    "flagella_presence",
)

# A genus needs this many typed rows before it may outvote one of its members,
# and this share of them must agree. Deliberately the same numbers as
# gram_stain_against_phylum, so the two read the same way. The threshold also
# makes a tie unreachable: two values tied at the top means each holds at most
# half, which cannot reach 90%.
MIN_GENUS_VOTES = 5
CONSENSUS_THRESHOLD = 0.90

# Values that answer "we do not know" or "it depends". They neither vote nor get
# judged -- counting them as dissent would flag honesty.
NON_ANSWERS = ("", "variable", "uncharacterized", "unknown", "n/a")

# Rows whose dissent has been looked at and is not a defect, keyed on
# (trait, genus_key, species_key) because `taxa.id` is not portable between the
# two stores -- Aspergillus nidulans is 725 in MySQL and 724 in the mirror.
#
# This list is the point of the whole design. The honest answer to "a value that
# disagrees with its genus" is usually "fix it", but not always, and a check
# that cannot say "looked at, still correct" forces the wrong repair onto the
# interesting rows.
APPROVED_DISSENTERS: dict[tuple[str, str, str], str] = {
    # Clostridium is a gram-positive genus and these two are gram-positive
    # organisms, but both are documented as staining gram-negative, which is
    # what a stain result records. MiMeDB's own vocabulary has the category --
    # it maps "structurally positive but stains negative" onto 'variable' -- so
    # the distinction is one this schema already recognises.
    #
    # Not rewritten to 'variable', though that is arguably the better value:
    # these rows carry 4,143 and 4,336 abundance rows and 43 and 31
    # associations, the stored 'negative' came from GMrepo rather than from the
    # inconsistent MiMeDB import, and re-curating a trait on that much evidence
    # is a decision for a curator, not a side effect of adding a check. Flagged
    # here so the question is recorded rather than silently settled either way.
    ("gram_stain", "clostridium", "symbiosum"):
        "gram-positive organism that stains gram-negative; 'variable' may be "
        "the better value -- curation decision, 4,143 abundances",
    ("gram_stain", "clostridium", "innocuum"):
        "gram-positive organism that stains gram-negative; 'variable' may be "
        "the better value -- curation decision, 4,336 abundances",
}


def approved_for(trait: str) -> tuple[tuple[str, str], ...]:
    """The (genus_key, species_key) pairs excused from this trait's check."""
    return tuple(
        (genus_key, species_key)
        for (approved_trait, genus_key, species_key) in APPROVED_DISSENTERS
        if approved_trait == trait
    )


def _quoted(values: Iterable[str]) -> str:
    return ", ".join("'" + v.replace("'", "''") + "'" for v in values)


def genus_dissenter_sql(trait: str) -> str:
    """SQL counting rows whose `trait` contradicts their genus's consensus.

    Written in CTEs with no vendor functions so the same text runs against
    MySQL and the SQLite mirror; `gram_stain_against_phylum` had to be rewritten
    once for exactly that reason.
    """
    if trait not in GENUS_CONSTRAINED_TRAITS:
        raise ValueError(f"{trait} is not genus-constrained; see this module's docstring")
    non_answers = _quoted(NON_ANSWERS)
    excused = "".join(
        f"\n              AND NOT (x.genus_key = '{genus_key}' "
        f"AND x.species_key = '{species_key}')"
        for genus_key, species_key in approved_for(trait)
    )
    return f"""
            WITH votes AS (
                SELECT genus_key, LOWER(TRIM({trait})) AS v, COUNT(*) AS n
                FROM taxa
                WHERE genus_key IS NOT NULL AND TRIM(genus_key) <> ''
                  AND {trait} IS NOT NULL
                  AND LOWER(TRIM({trait})) NOT IN ({non_answers})
                GROUP BY genus_key, LOWER(TRIM({trait}))
            ),
            totals AS (
                SELECT genus_key, SUM(n) AS total, MAX(n) AS top
                FROM votes GROUP BY genus_key
            ),
            majority AS (
                SELECT v.genus_key, v.v
                FROM votes v
                JOIN totals t ON t.genus_key = v.genus_key AND v.n = t.top
                WHERE t.total >= {MIN_GENUS_VOTES}
                  AND CAST(t.top AS REAL) / t.total >= {CONSENSUS_THRESHOLD}
            )
            SELECT COUNT(*) FROM taxa x
            JOIN majority m ON m.genus_key = x.genus_key
            WHERE x.{trait} IS NOT NULL
              AND LOWER(TRIM(x.{trait})) NOT IN ({non_answers})
              AND LOWER(TRIM(x.{trait})) <> m.v{excused}
        """


def genus_consensus_checks() -> dict[str, str]:
    """One validate check per genus-constrained trait."""
    return {
        f"{trait}_contradicts_its_genus": genus_dissenter_sql(trait)
        for trait in GENUS_CONSTRAINED_TRAITS
    }


def load_genus_consensus(cursor, placeholder: str = "%s") -> dict[tuple[str, str], str]:
    """Return {(genus_key, trait): agreed value} for every genus that agrees.

    The Python half of the same rule, for callers that must decide before
    writing rather than count afterwards.
    """
    consensus: dict[tuple[str, str], str] = {}
    for trait in GENUS_CONSTRAINED_TRAITS:
        cursor.execute(
            f"SELECT genus_key, LOWER(TRIM({trait})), COUNT(*) FROM taxa "
            f"WHERE genus_key IS NOT NULL AND TRIM(genus_key) <> '' "
            f"  AND {trait} IS NOT NULL "
            f"  AND LOWER(TRIM({trait})) NOT IN ({_quoted(NON_ANSWERS)}) "
            f"GROUP BY genus_key, LOWER(TRIM({trait}))"
        )
        tallies: dict[str, dict[str, int]] = {}
        for genus_key, value, n in cursor.fetchall():
            tallies.setdefault(genus_key, {})[value] = int(n)
        for genus_key, counts in tallies.items():
            total = sum(counts.values())
            value, top = max(counts.items(), key=lambda kv: kv[1])
            if total >= MIN_GENUS_VOTES and top / total >= CONSENSUS_THRESHOLD:
                consensus[(genus_key, trait)] = value
    return consensus


def dissents_from_genus(
    consensus: dict[tuple[str, str], str],
    genus_key: Any,
    trait: str,
    value: Any,
) -> str | None:
    """The genus's agreed value when `value` contradicts it, else None.

    A trait the genus does not constrain, an unknown-ish value, or a genus with
    no consensus all return None: silence means "no objection", never "checked
    and fine".
    """
    if trait not in GENUS_CONSTRAINED_TRAITS:
        return None
    text = ("" if value is None else str(value)).strip().lower()
    if text in NON_ANSWERS:
        return None
    key = (("" if genus_key is None else str(genus_key)).strip().lower(), trait)
    agreed = consensus.get(key)
    if agreed is None or agreed == text:
        return None
    return agreed
