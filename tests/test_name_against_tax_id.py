"""A row's name and its tax ID must describe the same organism.

The case this exists for is already fixed in the database, so it cannot fire
there any more: `Vibrio fischeri` carried 29494, which is *Vibrio furnissii*, in
both stores. These keep the check honest about that shape, and about the much
larger set of rows it must NOT flag -- a stale-but-real name is this project's
normal state, not a defect.
"""
import importlib.util
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "check_name_against_tax_id",
    Path(__file__).resolve().parents[1] / "scripts" / "check_name_against_tax_id.py")
check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)


def pairs(*names):
    """The dump's (genus, epithet) pairs for one node."""
    from gutdb.transform import parse_scientific_name
    out = set()
    for name in names:
        genus, epithet = parse_scientific_name(check.normalize(name))
        out.add((genus.casefold(), epithet.casefold()))
    return out


class NameAgainstTaxIdTests(unittest.TestCase):
    def test_the_vibrio_fischeri_shape_is_caught(self):
        # 29494 is Vibrio furnissii. A row named fischeri holding it is a wrong
        # tax ID, and nothing else in the project notices unless two rows
        # collide on the ID.
        self.assertEqual(
            check.classify("Vibrio fischeri", "vibrio", "fischeri", "species",
                           pairs("Vibrio furnissii", "Vibrio furnissii NCTC 11218")),
            "DISAGREES")

    def test_a_synonym_of_the_node_agrees(self):
        # The trap: NCBI's scientific name for 668 is Aliivibrio fischeri, and
        # comparing against that alone would flag this correct row.
        self.assertEqual(
            check.classify("Vibrio fischeri", "vibrio", "fischeri", "species",
                           pairs("Aliivibrio fischeri", "Vibrio fischeri")),
            "agrees_exactly")

    def test_a_strain_designation_after_the_epithet_agrees(self):
        self.assertEqual(
            check.classify("Clostridium difficile", "clostridium", "difficile",
                           "species", pairs("Clostridium difficile CD196")),
            "agrees_on_the_epithet")

    def test_a_drifted_genus_with_the_same_epithet_agrees(self):
        # merge_duplicate_taxa.py canonicalizes this; it is not a wrong tax ID.
        self.assertEqual(
            check.classify("Prevotella ruminicola", "prevotella", "ruminicola",
                           "species", pairs("Xylanibacter ruminicola")),
            "agrees_on_the_epithet")

    def test_latin_gender_agrees(self):
        self.assertEqual(
            check.classify("Candida glabrata", "candida", "glabrata", "species",
                           pairs("Nakaseomyces glabratus")),
            "agrees_on_the_epithet")

    def test_quotes_and_parentheses_do_not_cause_a_mismatch(self):
        self.assertEqual(
            check.classify("Candidatus Fonsibacter ubiquis", "candidatus fonsibacter",
                           "ubiquis", "species",
                           pairs('"Candidatus Fonsibacter ubiquis" Henson et al. 2018')),
            # The authority trails the epithet, so this agrees on the epithet
            # rather than on the whole string.
            "agrees_on_the_epithet")

    def test_underscores_and_brackets_do_not_cause_a_mismatch(self):
        self.assertEqual(
            check.classify("Bacteroides sp. 1 1 14", "bacteroides", "sp. 1 1 14",
                           "species", pairs("Bacteroides sp. 1_1_14")),
            "agrees_exactly")
        self.assertEqual(
            check.classify("Clostridium leptum", "clostridium", "leptum", "species",
                           pairs("[Clostridium] leptum")),
            "agrees_exactly")

    def test_a_bare_genus_agrees_on_the_genus_alone(self):
        self.assertEqual(
            check.classify("Blautia", "blautia", "", "genus", pairs("Blautia")),
            "agrees_exactly")

    def test_a_node_with_only_type_material_is_not_judged(self):
        self.assertEqual(
            check.classify("Streptomyces albosporeus", "streptomyces",
                           "albosporeus", "species", pairs("DSM 41672")),
            "no_binomial_at_this_node")

    def test_a_tax_id_absent_from_the_dump_is_not_judged(self):
        self.assertEqual(
            check.classify("Mammaliicoccus sciuri", "mammaliicoccus", "sciuri",
                           "species", None),
            "tax_id_not_in_dump")


if __name__ == "__main__":
    unittest.main()
