"""A source row whose tax ID the table already holds must update, not duplicate.

MiMeDB offers `Bacteroides vulgatus` for tax 821, which this database holds as
`Phocaeicola vulgatus`. Inserting beside it splits one organism across two rows
and two ML features. The guards are what make the redirect safe, so most of
these tests are about when it must NOT happen.
"""
import unittest

from gutdb.pipeline import _species_taxid, _taxid_holder


class Dump:
    """Stands in for the taxdump: strain -> species, unknown -> None."""

    def __init__(self, mapping):
        self.mapping = mapping

    def species_node(self, taxid):
        return self.mapping.get(taxid)


class Connection:
    """Mock holding one taxa table's worth of (id, species, rank) by tax ID."""

    def __init__(self, by_taxid):
        self.by_taxid = by_taxid
        self.pending = []

    def cursor(self, **kwargs):
        return self

    def execute(self, sql, params):
        self.pending = self.by_taxid.get(params[0], [])

    def fetchall(self):
        return self.pending

    def close(self):
        pass


class SpeciesTaxidTests(unittest.TestCase):
    def test_strain_resolves_to_its_species(self):
        # 435590 is "Phocaeicola vulgatus ATCC 8482"; the species is 821.
        dump = Dump({435590: 821})
        self.assertEqual(_species_taxid("435590", dump), 821)

    def test_species_is_left_alone(self):
        self.assertEqual(_species_taxid(821, Dump({821: 821})), 821)

    def test_unknown_tax_id_is_kept_as_given(self):
        # NCBI deletes nodes (1232666, 1328388). Guessing is worse than keeping
        # what the source said.
        self.assertEqual(_species_taxid(1232666, Dump({})), 1232666)

    def test_without_a_dump_nothing_changes(self):
        self.assertEqual(_species_taxid(435590, None), 435590)

    def test_blank_survives(self):
        self.assertEqual(_species_taxid("", Dump({})), "")


class TaxidHolderTests(unittest.TestCase):
    def test_redirects_onto_the_row_holding_the_tax_id(self):
        conn = Connection({821: [(1569, "vulgatus", "species")]})
        self.assertEqual(_taxid_holder(conn, 821, "species", "vulgatus"), 1569)

    def test_latin_gender_still_counts_as_the_same_organism(self):
        # Candida glabrata -> Nakaseomyces glabratus is one organism renamed.
        conn = Connection({5478: [(11554, "glabratus", "species")]})
        self.assertEqual(_taxid_holder(conn, 5478, "species", "glabrata"), 11554)

    def test_a_different_epithet_is_refused(self):
        # tax 411483 is "F. prausnitzii A2-165", renamed F. duncaniae. A source
        # row carrying it must not relabel F. prausnitzii, so this inserts.
        conn = Connection({411483: [(17, "prausnitzii", "species")]})
        self.assertIsNone(_taxid_holder(conn, 411483, "species", "duncaniae"))

    def test_a_different_rank_is_refused(self):
        # Merging across ranks broke 663 rank groups once.
        conn = Connection({1239: [(900, "", "genus")]})
        self.assertIsNone(_taxid_holder(conn, 1239, "species", "firmus"))

    def test_two_holders_are_left_to_merge_duplicate_taxa(self):
        conn = Connection({821: [(1569, "vulgatus", "species"),
                                 (2152, "vulgatus", "species")]})
        self.assertIsNone(_taxid_holder(conn, 821, "species", "vulgatus"))

    def test_no_holder_means_insert(self):
        self.assertIsNone(_taxid_holder(Connection({}), 821, "species", "vulgatus"))

    def test_no_tax_id_means_insert(self):
        conn = Connection({821: [(1569, "vulgatus", "species")]})
        self.assertIsNone(_taxid_holder(conn, None, "species", "vulgatus"))
        self.assertIsNone(_taxid_holder(conn, "", "species", "vulgatus"))

    def test_a_genus_row_is_never_redirected_by_epithet(self):
        # A bare genus has no epithet to compare, so there is nothing to make
        # the match safe; it inserts.
        conn = Connection({1239: [(900, "", "genus")]})
        self.assertIsNone(_taxid_holder(conn, 1239, "genus", ""))


if __name__ == "__main__":
    unittest.main()
