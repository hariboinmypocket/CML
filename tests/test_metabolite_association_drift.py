"""The stored succinate associations must still be what the data produces.

taxon_metabolite_associations is DERIVED and nothing re-derives it, which is how
taxon_disease_adjudication lost a verdict: a taxa merge renamed a genus, the
foreign key cascade removed the row keyed on the old taxon_id, and no counter
noticed because every counter asks whether the rows present are consistent, not
whether a row that should exist does.

This table has the same `fk_tma_taxon ... ON DELETE CASCADE`, and the same
blindness was measured rather than assumed: deleting Bacteroides' 11 association
rows -- exactly what merging that genus would do -- moved none of the 76
counters. Bacteroides is also the genus most central to succinate, so the most
consequential rows were the unguarded ones.

A counter cannot close this. Whether a taxon SHOULD have a row depends on the
prevalence floor inside clr_matrix and on the support floors here, so any SQL
approximation would flag legitimately-excluded taxa as missing. Re-deriving is
the only exact check, and it costs about 3s for all eleven studies, so it is
affordable as a test rather than something to remember to run.
"""
import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

SUCCINATE_KEGG = "C00042"


def load_deriver():
    path = ROOT / "scripts" / "derive_taxon_metabolite_associations.py"
    spec = importlib.util.spec_from_file_location("derive_tma_drift", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def connection():
    try:
        from gutdb.config import Settings
        from gutdb.db import connect
        return connect(Settings.from_env())
    except Exception:
        return None


class MetaboliteAssociationDriftTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.connection = connection()
        if cls.connection is None:
            raise unittest.SkipTest("no database reachable")
        cls.deriver = load_deriver()
        cursor = cls.connection.cursor()
        try:
            cursor.execute("SELECT id FROM metabolites WHERE kegg_id = %s",
                           (SUCCINATE_KEGG,))
            found = cursor.fetchall()
        finally:
            cursor.close()
        if len(found) != 1:
            raise unittest.SkipTest(
                f"kegg_id {SUCCINATE_KEGG} is held by {len(found)} rows")
        cls.metabolite_id = found[0][0]
        # Once per class, not per test: the derivation is about 3s and all four
        # tests read the same two snapshots.
        cls.stored_rows = cls.stored(cls)
        cls.fresh_rows = cls.fresh(cls)

    @classmethod
    def tearDownClass(cls):
        if getattr(cls, "connection", None) is not None:
            cls.connection.close()

    def stored(self):
        cursor = self.connection.cursor()
        try:
            cursor.execute(
                """SELECT taxon_id, study_id, rho, n_samples, n_nonzero,
                          high_confidence_annotation
                   FROM taxon_metabolite_associations WHERE metabolite_id = %s""",
                (self.metabolite_id,))
            return {(r[0], r[1]): r for r in cursor.fetchall()}
        finally:
            cursor.close()

    def fresh(self):
        levels, confidence = self.deriver.metabolite_by_study(
            self.connection, self.metabolite_id)
        return {(r[0], r[1]): r
                for r in self.deriver.correlate(self.connection, levels,
                                                confidence, "genus")}


    def test_no_derived_pair_is_missing_from_the_table(self):
        """The cascade's failure mode: rows the data supports are simply gone."""
        missing = sorted(set(self.fresh_rows) - set(self.stored_rows))
        self.assertEqual(missing, [], f"{len(missing)} derived pairs absent from "
                                      "taxon_metabolite_associations")

    def test_no_stored_pair_has_lost_its_evidence(self):
        extra = sorted(set(self.stored_rows) - set(self.fresh_rows))
        self.assertEqual(extra, [], f"{len(extra)} stored pairs the data no "
                                    "longer produces")

    def test_the_stored_values_still_match_a_fresh_derivation(self):
        drifted = []
        for key in sorted(set(self.fresh_rows) & set(self.stored_rows)):
            fresh, stored = self.fresh_rows[key], self.stored_rows[key]
            if (round(float(fresh[2]), 4) != round(float(stored[2]), 4)
                    or fresh[5] != stored[3]
                    or fresh[6] != stored[4]
                    or bool(fresh[7]) != bool(stored[5])):
                drifted.append(key)
        self.assertEqual(drifted, [], f"{len(drifted)} pairs differ from a fresh run")

    def test_every_stored_taxon_still_exists(self):
        """Cheap, and it would catch a cascade that somehow left orphans."""
        cursor = self.connection.cursor()
        try:
            cursor.execute("""SELECT COUNT(*) FROM taxon_metabolite_associations a
                              WHERE NOT EXISTS (SELECT 1 FROM taxa t
                                                WHERE t.id = a.taxon_id)""")
            self.assertEqual(int(cursor.fetchone()[0]), 0)
        finally:
            cursor.close()


if __name__ == "__main__":
    unittest.main()
