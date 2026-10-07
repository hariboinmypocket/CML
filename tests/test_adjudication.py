"""Tests for the vote arithmetic in scripts/adjudicate_conflicts.py.

The two functions under test are pure, and both encode a decision that a
previous version of the evidence view got wrong: a study may cast at most one
vote, and only case/control contrasts may vote at all.
"""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from adjudicate_conflicts import adjudicate, association_vote


def ab(enriched=0, depleted=0, ties=0):
    total = enriched + depleted
    ratio = (max(enriched, depleted) / total) if total else None
    return {
        "ab_enriched": enriched, "ab_depleted": depleted, "ab_ties": ties,
        "ab_studies": total, "ab_ratio": round(ratio, 3) if ratio else None,
        "ab_direction": (None if not total or enriched == depleted
                         else "enriched" if enriched > depleted else "depleted"),
    }


class AssociationVote(unittest.TestCase):
    def test_one_vote_per_study(self):
        vote = association_vote({1: {"enriched"}, 2: {"enriched"}, 3: {"depleted"}})
        self.assertEqual((vote["assoc_enriched"], vote["assoc_depleted"]), (2, 1))
        self.assertEqual(vote["assoc_consensus"], "enriched")
        self.assertEqual(vote["assoc_studies"], 3)

    def test_study_disagreeing_with_itself_votes_for_neither(self):
        """The old view counted such a study on both sides at once."""
        vote = association_vote({1: {"enriched", "depleted"}, 2: {"depleted"}})
        self.assertEqual((vote["assoc_enriched"], vote["assoc_depleted"]), (0, 1))
        self.assertEqual(vote["assoc_split_studies"], 1)
        self.assertEqual(vote["assoc_consensus"], "depleted")
        self.assertEqual(vote["assoc_agreement"], 1.0)

    def test_no_case_control_studies_is_not_a_tie(self):
        vote = association_vote({})
        self.assertEqual(vote["assoc_studies"], 0)
        self.assertIsNone(vote["assoc_consensus"])
        self.assertIsNone(vote["assoc_agreement"])

    def test_exact_tie_has_no_consensus(self):
        vote = association_vote({1: {"enriched"}, 2: {"depleted"}})
        self.assertIsNone(vote["assoc_consensus"])
        self.assertEqual(vote["assoc_agreement"], 0.5)


class Adjudicate(unittest.TestCase):
    def test_agreement_is_recorded_as_agreement(self):
        assoc = association_vote({1: {"enriched"}, 2: {"enriched"}, 3: {"enriched"}})
        self.assertEqual(adjudicate(assoc, ab(enriched=5, depleted=1)),
                         ("enriched", "association_and_abundance_agree"))

    def test_majority_association_is_not_overruled_but_is_flagged(self):
        assoc = association_vote({1: {"enriched"}, 2: {"enriched"}, 3: {"enriched"}})
        verdict, basis = adjudicate(assoc, ab(enriched=0, depleted=6))
        self.assertEqual(verdict, "enriched")
        self.assertEqual(basis, "association_consensus_abundance_disagrees")

    def test_abundance_breaks_an_association_tie(self):
        assoc = association_vote({1: {"enriched"}, 2: {"depleted"}})
        self.assertEqual(adjudicate(assoc, ab(enriched=15, depleted=5)),
                         ("enriched", "abundance_majority"))

    def test_tie_plus_tie_stays_unresolved(self):
        assoc = association_vote({1: {"enriched"}, 2: {"depleted"}})
        self.assertEqual(adjudicate(assoc, ab(enriched=3, depleted=3)),
                         (None, "unresolved_abundance_also_tied"))

    def test_tie_plus_thin_abundance_stays_unresolved(self):
        assoc = association_vote({1: {"enriched"}, 2: {"depleted"}})
        self.assertEqual(adjudicate(assoc, ab(enriched=1)),
                         (None, "unresolved_no_abundance_evidence"))

    def test_tie_plus_flat_abundance_majority_stays_unresolved(self):
        """A 6-4 abundance split is not a finding."""
        assoc = association_vote({1: {"enriched"}, 2: {"depleted"}})
        self.assertEqual(adjudicate(assoc, ab(enriched=6, depleted=4)),
                         (None, "unresolved_abundance_inconclusive"))

    def test_disease_vs_disease_only_falls_back_to_abundance(self):
        assoc = association_vote({})
        self.assertEqual(adjudicate(assoc, ab(enriched=9, depleted=1)),
                         ("enriched", "abundance_only_no_case_control_association"))

    def test_no_evidence_at_all_decides_nothing(self):
        self.assertEqual(adjudicate(association_vote({}), ab()),
                         (None, "unresolved_no_case_control_evidence"))

    def test_settled_association_without_abundance_stands_alone(self):
        assoc = association_vote({1: {"depleted"}, 2: {"depleted"}})
        self.assertEqual(adjudicate(assoc, ab(enriched=1)),
                         ("depleted", "association_consensus"))


if __name__ == "__main__":
    unittest.main()
