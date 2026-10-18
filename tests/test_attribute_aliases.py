"""Tests for the attribute alias seed.

These check the seed's internal consistency rather than the database, which is
what can go wrong silently here: a boolean concept whose values are not all
decoded would let a modeller pool a key with an inverted encoding, and that is
the error the tables exist to prevent. The seed lists are imported directly, so
these run without a database.
"""
import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

_spec = importlib.util.spec_from_file_location(
    "seed_attribute_aliases", ROOT / "scripts" / "seed_attribute_aliases.py")
_seed = importlib.util.module_from_spec(_spec)
sys.argv = [sys.argv[0]]          # keep the module's DRY_RUN detection quiet
_spec.loader.exec_module(_seed)

ALIASES = _seed.ALIASES
VALUE_ALIASES = _seed.VALUE_ALIASES
BY_ATTR = {a: (c, k, n) for a, c, k, n in ALIASES}
VALUES_BY_ATTR: dict[str, dict[str, str]] = {}
for attr, raw, canon in VALUE_ALIASES:
    VALUES_BY_ATTR.setdefault(attr, {})[raw] = canon


class SeedShapeTests(unittest.TestCase):
    def test_no_attribute_is_listed_twice(self):
        attrs = [a for a, *_ in ALIASES]
        self.assertEqual(len(attrs), len(set(attrs)))

    def test_every_value_alias_has_a_key_alias(self):
        """An orphan value alias would violate the foreign key at load time."""
        for attr, raw, _ in VALUE_ALIASES:
            self.assertIn(attr, BY_ATTR, f"{attr} / {raw}")

    def test_no_value_alias_is_listed_twice(self):
        pairs = [(a, r) for a, r, _ in VALUE_ALIASES]
        self.assertEqual(len(pairs), len(set(pairs)))

    def test_value_kinds_are_from_the_enum(self):
        allowed = {"boolean", "category", "quantity", "interval", "agent",
                   "date", "free_text"}
        for attr, _, kind, _ in ALIASES:
            self.assertIn(kind, allowed, attr)


class BooleanDecodingTests(unittest.TestCase):
    """The part that matters: a boolean must never be read from its raw value."""

    def test_every_boolean_key_carries_a_decoding(self):
        for attr, concept, kind, _ in ALIASES:
            if kind == "boolean":
                self.assertIn(attr, VALUES_BY_ATTR,
                              f"{attr} is boolean with no value map")

    def test_booleans_decode_only_to_yes_or_no(self):
        for attr, values in VALUES_BY_ATTR.items():
            if BY_ATTR[attr][1] != "boolean":
                continue
            for raw, canon in values.items():
                self.assertIn(canon, {"yes", "no"}, f"{attr} {raw!r} -> {canon}")

    def test_the_inverted_key_is_decoded_against_its_name(self):
        """'1' means NO under this key and YES under Antibiotics.

        Pooling them on name alone inverts 293 rows, measured. This is the single
        assertion the whole value-alias table exists to make true.
        """
        inverted = VALUES_BY_ATTR["Antibiotic_in_the_past_month._1_no_2_yes"]
        self.assertEqual(inverted["1"], "no")
        self.assertEqual(inverted["2"], "yes")
        self.assertEqual(VALUES_BY_ATTR["Antibiotics"]["1"], "yes")
        self.assertEqual(VALUES_BY_ATTR["Antibiotics"]["0"], "no")

    def test_the_other_1_no_2_yes_key_is_inverted_too(self):
        probiotic = VALUES_BY_ATTR["X_probiotic_in_the_past_2_weeks_1_no_2_yes"]
        self.assertEqual(probiotic["1"], "no")
        self.assertEqual(probiotic["2"], "yes")

    def test_a_key_whose_name_declares_a_code_is_never_read_as_a_number(self):
        """These are the keys fix_attribute_numerics.py clears value_numeric for."""
        from gutdb.transform import numeric_value_is_trustworthy
        for attr in ("Antibiotic_in_the_past_month._1_no_2_yes",
                     "X_probiotic_in_the_past_2_weeks_1_no_2_yes"):
            self.assertFalse(numeric_value_is_trustworthy(attr, "1"), attr)
            self.assertEqual(BY_ATTR[attr][1], "boolean")


class ConceptSeparationTests(unittest.TestCase):
    """A concept is the question, not the topic."""

    def test_antibiotics_keep_three_separate_concepts(self):
        concepts = {c for a, c, k, n in ALIASES if "ntibio" in a or "bx" in a}
        self.assertEqual(concepts, {"antibiotic_exposure", "antibiotic_agent",
                                    "antibiotic_start_date"})

    def test_an_agent_or_date_concept_holds_no_boolean(self):
        for attr, concept, kind, _ in ALIASES:
            if concept in ("antibiotic_agent", "antibiotic_start_date",
                           "probiotic_product"):
                self.assertNotEqual(kind, "boolean", attr)

    def test_the_false_clusters_are_absent(self):
        """Each was produced by a name heuristic and rejected by the values.

        steroids/mesalamine are distinct drugs, Breast cancer is a diagnosis,
        and the *_Level keys are five different autoantibodies.
        """
        for attr in ("steroids", "mesalamine", "Breast cancer", "GADA_Level",
                     "IA2A_Level", "IAA_Level", "AAB_Level", "ICA_Level"):
            self.assertNotIn(attr, BY_ATTR, attr)

    def test_smoking_status_and_exposure_are_different_concepts(self):
        self.assertEqual(BY_ATTR["Brinkman Index"][0], "smoking_exposure_index")
        self.assertEqual(BY_ATTR["SmokingStatus"][0], "smoking_status")
        self.assertEqual(BY_ATTR["Brinkman Index"][1], "quantity")

    def test_a_concept_may_legitimately_mix_value_kinds(self):
        """Which is why the coverage view reports every kind, not one.

        `infant_feeding` has a category key beside two booleans and
        `alcohol_amount` has a quantity beside an interval. A view taking MIN()
        over the kinds reported just 'boolean' for the first -- a column
        answering a different question than its name.
        """
        kinds = {}
        for attr, concept, kind, _ in ALIASES:
            kinds.setdefault(concept, set()).add(kind)
        self.assertEqual(kinds["infant_feeding"], {"boolean", "category"})
        self.assertEqual(kinds["alcohol_amount"], {"quantity", "interval"})
        self.assertEqual(kinds["antibiotic_exposure"], {"boolean"})

    def test_an_interval_key_is_not_called_a_quantity(self):
        """'Alcohol.drink.week.' holds ranges; its numeric subset is biased."""
        self.assertEqual(BY_ATTR["Alcohol.drink.week."][1], "interval")
        self.assertEqual(BY_ATTR["Alcohol"][1], "quantity")


if __name__ == "__main__":
    unittest.main()
