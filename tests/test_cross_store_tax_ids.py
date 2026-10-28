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


if __name__ == "__main__":
    unittest.main()
