"""Tests for gutdb.traits: the genus-consensus rule and the SQL built from it.

The SQL is pinned against a planted dissenter rather than only asserted to
return 0 against the live data. A check that reads 0 because its query is
broken looks exactly like a check that reads 0 because the data is clean, and
this repository has already been bitten by a counter that gated nothing.
"""
from __future__ import annotations

import sqlite3
import unittest

from gutdb.traits import (APPROVED_DISSENTERS, CONSENSUS_THRESHOLD,
                          GENUS_CONSTRAINED_TRAITS, MIN_GENUS_VOTES,
                          dissents_from_genus, genus_consensus_checks,
                          genus_dissenter_sql, load_genus_consensus)


def _db(rows, trait="gram_stain"):
    """An in-memory taxa table holding (genus_key, species_key, trait) rows.

    Every constrained trait gets a column, because load_genus_consensus reads
    all of them in one pass; only `trait` is populated.
    """
    con = sqlite3.connect(":memory:")
    columns = ", ".join(f"{t} TEXT" for t in GENUS_CONSTRAINED_TRAITS)
    con.execute(
        f"CREATE TABLE taxa (id INTEGER PRIMARY KEY, scientific_name TEXT, "
        f"genus_key TEXT, species_key TEXT, {columns})")
    con.executemany(
        f"INSERT INTO taxa (scientific_name, genus_key, species_key, {trait}) "
        f"VALUES (?, ?, ?, ?)",
        [(f"{g} {s}", g, s, v) for g, s, v in rows])
    con.commit()
    return con


class PredicateTests(unittest.TestCase):
    def test_a_value_against_its_genus_reports_the_agreed_value(self):
        consensus = {("clostridium", "gram_stain"): "positive"}
        self.assertEqual(
            dissents_from_genus(consensus, "clostridium", "gram_stain", "negative"),
            "positive")

    def test_agreement_is_silent(self):
        consensus = {("clostridium", "gram_stain"): "positive"}
        self.assertIsNone(
            dissents_from_genus(consensus, "clostridium", "gram_stain", "positive"))

    def test_an_unknown_answer_is_not_dissent(self):
        """'variable' and friends answer the question; they do not contradict it."""
        consensus = {("clostridium", "gram_stain"): "positive"}
        for value in ("variable", "uncharacterized", "", None, "  "):
            self.assertIsNone(
                dissents_from_genus(consensus, "clostridium", "gram_stain", value),
                f"{value!r} should not count as dissent")

    def test_a_trait_the_genus_does_not_constrain_is_never_flagged(self):
        """oxygen_requirement peaks at 59% agreement, so it is out by design."""
        self.assertNotIn("oxygen_requirement", GENUS_CONSTRAINED_TRAITS)
        self.assertNotIn("temperature_range", GENUS_CONSTRAINED_TRAITS)
        consensus = {("clostridium", "oxygen_requirement"): "anaerobe"}
        self.assertIsNone(
            dissents_from_genus(consensus, "clostridium", "oxygen_requirement",
                                "aerobe"))

    def test_a_genus_with_no_consensus_is_silent(self):
        self.assertIsNone(
            dissents_from_genus({}, "nosuchgenus", "gram_stain", "negative"))

    def test_case_and_padding_do_not_create_dissent(self):
        consensus = {("clostridium", "gram_stain"): "positive"}
        self.assertIsNone(
            dissents_from_genus(consensus, " Clostridium ", "gram_stain", " POSITIVE "))

    def test_an_unchecked_trait_cannot_have_sql_built_for_it(self):
        with self.assertRaises(ValueError):
            genus_dissenter_sql("oxygen_requirement")

    def test_a_check_exists_for_every_constrained_trait(self):
        checks = genus_consensus_checks()
        self.assertEqual(len(checks), len(GENUS_CONSTRAINED_TRAITS))
        for trait in GENUS_CONSTRAINED_TRAITS:
            self.assertIn(f"{trait}_contradicts_its_genus", checks)


class ConsensusLoadingTests(unittest.TestCase):
    def test_a_genus_below_the_vote_floor_has_no_consensus(self):
        rows = [("smallgenus", f"sp{i}", "positive") for i in range(MIN_GENUS_VOTES - 1)]
        con = _db(rows)
        self.assertEqual(load_genus_consensus(con.cursor(), "?"), {})
        con.close()

    def test_a_genus_at_the_floor_has_one(self):
        rows = [("biggenus", f"sp{i}", "positive") for i in range(MIN_GENUS_VOTES)]
        con = _db(rows)
        self.assertEqual(
            load_genus_consensus(con.cursor(), "?").get(("biggenus", "gram_stain")),
            "positive")
        con.close()

    def test_a_split_genus_has_none(self):
        """Below the threshold nobody outvotes anybody."""
        rows = ([("split", f"a{i}", "positive") for i in range(6)]
                + [("split", f"b{i}", "negative") for i in range(4)])
        con = _db(rows)
        self.assertLess(6 / 10, CONSENSUS_THRESHOLD)
        self.assertNotIn(("split", "gram_stain"),
                         load_genus_consensus(con.cursor(), "?"))
        con.close()


class DissenterSqlTests(unittest.TestCase):
    def _count(self, con, trait="gram_stain"):
        return con.execute(genus_dissenter_sql(trait)).fetchone()[0]

    def test_a_planted_dissenter_is_found(self):
        rows = [("plantedgenus", f"sp{i}", "positive") for i in range(19)]
        rows.append(("plantedgenus", "oddone", "negative"))
        con = _db(rows)
        self.assertEqual(self._count(con), 1)
        con.close()

    def test_a_clean_genus_counts_nothing(self):
        rows = [("cleangenus", f"sp{i}", "positive") for i in range(20)]
        con = _db(rows)
        self.assertEqual(self._count(con), 0)
        con.close()

    def test_an_approved_dissenter_is_excused(self):
        """The exception list is what lets a real dissenter stay without noise."""
        trait, genus_key, species_key = next(
            k for k in APPROVED_DISSENTERS if k[0] == "gram_stain")
        rows = [(genus_key, f"sp{i}", "positive") for i in range(19)]
        rows.append((genus_key, species_key, "negative"))
        con = _db(rows, trait)
        self.assertEqual(self._count(con, trait), 0)
        # The same row without its approval would be counted, which is what
        # makes the zero above mean "excused" rather than "query broken".
        con.execute("UPDATE taxa SET species_key = 'notapproved' "
                    "WHERE species_key = ?", (species_key,))
        self.assertEqual(self._count(con, trait), 1)
        con.close()

    def test_unknown_answers_neither_vote_nor_count(self):
        rows = [("mixedgenus", f"sp{i}", "positive") for i in range(19)]
        rows += [("mixedgenus", "unknown1", "variable"),
                 ("mixedgenus", "unknown2", "uncharacterized")]
        con = _db(rows)
        self.assertEqual(self._count(con), 0)
        con.close()

    def test_every_trait_s_sql_runs(self):
        """Each generated query must parse and execute, not just gram_stain's."""
        for trait in GENUS_CONSTRAINED_TRAITS:
            con = _db([("g", f"sp{i}", "yes") for i in range(5)], trait)
            self.assertEqual(self._count(con, trait), 0)
            con.close()


if __name__ == "__main__":
    unittest.main()
