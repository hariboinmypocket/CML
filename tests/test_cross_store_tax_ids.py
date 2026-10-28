"""The two stores must not disagree about which organism a row is.

This project's recurring defect is a repair that reached only one store, and
tax IDs are where it does the most damage: the ID is what joins a row to NCBI,
so a row can carry the right name, the right lineage and the right abundance
while pointing at a different organism entirely.

It was found here rather than reasoned about. `taxa` held `Oxalobacter
formigenes` with ncbi_tax_id 556268 in MySQL and 847 in the mirror. 847 is
`Oxalobacter formigenes`; 556268 is `Oxalobacter paraformigenes`, a species
split out of it in 2022 that lists only the strain `Oxalobacter formigenes
HOxBLS` under `includes`. So MySQL named one species and pointed at another
that still exists separately, with 686 abundance rows and 5 disease
associations hanging off it -- and the mirror was the store that had it right.

No counter could see this. `validate_sqlite_mirror` runs the same check SQL
against each store separately and prefixes the results `sqlite_`, so a value
that is individually legal in both stores but different between them passes
every check. Comparing the stores needs both connections at once, which is a
test rather than a counter, on the same grounds as the association drift check.

Rows are matched on (genus_key, species_key), never on id: the mirror keeps its
own taxa.id space, so an id-keyed comparison would compare two unrelated
organisms and could report agreement that does not exist.
"""
import sys
import unittest
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

# (genus_key, species_key) -> why this one still differs, and the decision it
# waits on. Every entry is asserted to be STILL divergent, so correcting one
# forces its removal here instead of leaving a stale exemption behind.
KNOWN_DIVERGENT = {
    ("phocaeicola", "coprocola"):
        "MySQL carries the strain id 470145, the mirror the species 310298. "
        "Retagging MySQL onto 310298 would collide with Bacteroides coprocola "
        "(id 781) in the same species layer -- the pre-rename row for the same "
        "organism, holding 1 association, so not the evidence-free synonym "
        "fix_misassigned_tax_ids.py is willing to drop. No sample holds both "
        "rows, so folding them is additive-safe, but moving an association "
        "between taxa is a data decision, and refusing beats guessing.",
}

TAXA_SQL = """SELECT genus_key, species_key, ncbi_tax_id, ncbi_rank, scientific_name
              FROM taxa WHERE ncbi_tax_id IS NOT NULL"""


def mysql_connection():
    try:
        from gutdb.config import Settings
        from gutdb.db import connect
        return connect(Settings.from_env())
    except Exception:
        return None


def mirror_connection():
    import sqlite3
    path = ROOT / "gut_microbiome.sqlite"
    if not path.exists():
        return None
    try:
        from gutdb.config import require_local
        return sqlite3.connect(require_local(path), timeout=60)
    except Exception:
        return None


def taxa_by_key(rows):
    """key -> {field: set of values}, so a duplicated key cannot hide a value."""
    out = defaultdict(lambda: {"tax_id": set(), "rank": set(), "name": set()})
    for genus_key, species_key, tax_id, rank, name in rows:
        entry = out[(genus_key, species_key)]
        entry["tax_id"].add(int(tax_id))
        entry["rank"].add(rank)
        entry["name"].add(name)
    return out


class CrossStoreTaxIdTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.connection = mysql_connection()
        if cls.connection is None:
            raise unittest.SkipTest("no database reachable")
        cls.mirror = mirror_connection()
        if cls.mirror is None:
            raise unittest.SkipTest("no readable SQLite mirror")
        cursor = cls.connection.cursor()
        try:
            cursor.execute(TAXA_SQL)
            cls.mysql_taxa = taxa_by_key(cursor.fetchall())
        finally:
            cursor.close()
        cls.mirror_taxa = taxa_by_key(cls.mirror.execute(TAXA_SQL).fetchall())
        # The mirror is a strict subset, so only shared keys can disagree.
        cls.shared = sorted(set(cls.mysql_taxa) & set(cls.mirror_taxa))

    @classmethod
    def tearDownClass(cls):
        for handle in ("connection", "mirror"):
            if getattr(cls, handle, None) is not None:
                getattr(cls, handle).close()

    def test_there_are_shared_taxa_to_compare(self):
        """Without this the two tests below pass by comparing nothing."""
        self.assertGreater(len(self.shared), 1000)

    def test_no_shared_taxon_points_at_a_different_organism(self):
        disagreements = []
        for key in self.shared:
            if key in KNOWN_DIVERGENT:
                continue
            mine = self.mysql_taxa[key]["tax_id"]
            theirs = self.mirror_taxa[key]["tax_id"]
            if mine != theirs:
                name = sorted(self.mysql_taxa[key]["name"])[0]
                disagreements.append(
                    "%s (%s/%s): MySQL %s vs mirror %s"
                    % (name, key[0], key[1], sorted(mine), sorted(theirs)))
        self.assertEqual(
            [], disagreements,
            "%d taxa carry a different NCBI tax ID in each store:\n  %s"
            % (len(disagreements), "\n  ".join(disagreements[:20])))

    def test_no_shared_taxon_disagrees_about_its_ncbi_rank(self):
        """A retag leaves ncbi_rank describing the id the row used to carry.

        That is the half-finished repair this project has already hit: MySQL
        kept 'strain' after Vibrio furnissii moved off a strain id.
        """
        disagreements = []
        for key in self.shared:
            if self.mysql_taxa[key]["tax_id"] != self.mirror_taxa[key]["tax_id"]:
                continue
            mine = self.mysql_taxa[key]["rank"]
            theirs = self.mirror_taxa[key]["rank"]
            if mine != theirs:
                name = sorted(self.mysql_taxa[key]["name"])[0]
                disagreements.append(
                    "%s: MySQL %s vs mirror %s" % (name, sorted(mine), sorted(theirs)))
        self.assertEqual(
            [], disagreements,
            "%d taxa disagree about ncbi_rank between the stores:\n  %s"
            % (len(disagreements), "\n  ".join(disagreements[:20])))


    def test_each_documented_exemption_is_still_divergent(self):
        """An exemption that has been fixed must be deleted, not left lying."""
        for key, why in KNOWN_DIVERGENT.items():
            with self.subTest(taxon="/".join(key)):
                self.assertIn(key, self.shared,
                              "%s is no longer in both stores; drop the "
                              "exemption" % (key,))
                self.assertNotEqual(
                    self.mysql_taxa[key]["tax_id"], self.mirror_taxa[key]["tax_id"],
                    "%s now agrees across the stores -- remove it from "
                    "KNOWN_DIVERGENT. Recorded reason: %s" % (key, why))


class OneSidedTaxIdTests(unittest.TestCase):
    """A row joinable to NCBI in one store and not the other.

    The blind spot in the class above, and in the first version of this file:
    its query selects `WHERE ncbi_tax_id IS NOT NULL`, so a key whose ID is
    missing on one side never enters both maps and is never compared. 120 rows
    were in that state -- 89 with the ID only in MySQL, 31 only in the mirror.
    114 are filled now; the rest are refused for reasons this asserts exist.
    """

    @classmethod
    def setUpClass(cls):
        root = ROOT / "taxdump"
        for required in ("names.dmp", "nodes.dmp"):
            if not (root / required).exists():
                raise unittest.SkipTest(f"{required} not present (gitignored)")
        cls.connection = mysql_connection()
        if cls.connection is None:
            raise unittest.SkipTest("no database reachable")
        cls.mirror = mirror_connection()
        if cls.mirror is None:
            raise unittest.SkipTest("no readable SQLite mirror")
        import importlib.util
        path = ROOT / "scripts" / "backfill_one_sided_taxids.py"
        spec = importlib.util.spec_from_file_location("backfill_one_sided_taxids", path)
        cls.mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.mod)
        # Decided by the script's own function, so the test cannot drift from it.
        cls.fills, cls.refusals, cls.total = cls.mod.analyse(
            cls.connection, cls.mirror, root)

    @classmethod
    def tearDownClass(cls):
        for handle in ("connection", "mirror"):
            if getattr(cls, handle, None) is not None:
                getattr(cls, handle).close()

    def test_nothing_is_merely_waiting_to_be_filled(self):
        pending = ["%s: %s <- %s" % (target, row["name"], tax)
                   for target, _key, row, tax, _rank in self.fills]
        self.assertEqual(
            [], pending,
            "%d row(s) could take a tax ID from the other store. Run "
            "scripts/backfill_one_sided_taxids.py --apply:\n  %s"
            % (len(pending), "\n  ".join(pending)))

    def test_every_remaining_gap_has_a_stated_reason(self):
        self.assertEqual(
            self.total, len(self.fills) + len(self.refusals),
            "a one-sided tax ID was neither filled nor refused")
        for key, message in self.refusals:
            with self.subTest(taxon="/".join(key)):
                self.assertTrue(message.strip(), "refused with an empty reason")


class KeySpellingTests(unittest.TestCase):
    """The two stores must not key one organism two different ways.

    MySQL held `Bacteroides sp. 1 1 14` where the mirror held `sp. 1_1_14`,
    which is NCBI's spelling. 24 keys differed that way, and because every
    two-store repair here matches on (genus_key, species_key), those rows were
    invisible to all of them -- including to the class above, which compares
    only keys present in both stores.
    """

    @classmethod
    def setUpClass(cls):
        cls.connection = mysql_connection()
        if cls.connection is None:
            raise unittest.SkipTest("no database reachable")
        cls.mirror = mirror_connection()
        if cls.mirror is None:
            raise unittest.SkipTest("no readable SQLite mirror")
        cursor = cls.connection.cursor()
        try:
            cursor.execute("SELECT genus_key, species_key FROM taxa")
            cls.mysql_keys = {(a, b) for a, b in cursor.fetchall()}
        finally:
            cursor.close()
        cls.mirror_keys = {(a, b) for a, b in cls.mirror.execute(
            "SELECT genus_key, species_key FROM taxa").fetchall()}

    @classmethod
    def tearDownClass(cls):
        for handle in ("connection", "mirror"):
            if getattr(cls, handle, None) is not None:
                getattr(cls, handle).close()

    @staticmethod
    def _squash(key):
        import re
        return (key[0], re.sub(r"[ _]+", " ", key[1]))

    def test_no_key_is_a_punctuation_variant_of_the_other_stores(self):
        mysql_squashed = {self._squash(k) for k in self.mysql_keys}
        mirror_squashed = {self._squash(k) for k in self.mirror_keys}
        variants = [("mirror", k) for k in self.mirror_keys - self.mysql_keys
                    if self._squash(k) in mysql_squashed]
        variants += [("mysql", k) for k in self.mysql_keys - self.mirror_keys
                     if self._squash(k) in mirror_squashed]
        self.assertEqual(
            [], variants,
            "%d key(s) differ between the stores only by spaces or underscores, "
            "so the same organism is keyed two ways:\n  %s"
            % (len(variants), "\n  ".join(repr(v) for v in variants[:20])))


if __name__ == "__main__":
    unittest.main()
