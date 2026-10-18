"""Bounded LDA scaling, automatic refresh, and preservation of source effects."""
from pathlib import Path
import sqlite3
import unittest

from gutdb.lda import install_normalized_lda, normalize_lda_score


class LDANormalizationTests(unittest.TestCase):
    def setUp(self):
        self.connection = sqlite3.connect(":memory:")
        schema = Path(__file__).parents[1] / "gutdb" / "schema_sqlite.sql"
        self.connection.executescript(schema.read_text())
        self.connection.execute("INSERT INTO studies (id, project_accession) VALUES (1, 'TEST')")
        self.connection.execute("INSERT INTO diseases (id, name, name_key) VALUES (1, 'Test', 'test')")
        self.connection.execute("""INSERT INTO phenotype_comparisons
            (id, study_id, phenotype_a_id, phenotype_b_id, comparison_key)
            VALUES (1, 1, 1, 1, 'test')""")

    def tearDown(self):
        self.connection.close()

    def insert(self, row_id, score, kind="LDA"):
        self.connection.execute("""INSERT INTO taxa
            (id, genus, genus_key, scientific_name) VALUES (?, ?, ?, ?)""",
            (row_id, f'Taxon{row_id}', f'taxon{row_id}', f'Taxon{row_id}'))
        self.connection.execute("""INSERT INTO taxon_disease_associations
            (id, taxon_id, disease_id, comparison_id, direction, effect_type, effect_size)
            VALUES (?, ?, 1, 1, 'marker', ?, ?)""", (row_id, row_id, kind, score))

    def test_bounded_symmetric_monotonic_and_zero(self):
        scores = [-100, -5.62, -4, -1, 0, 1, 4, 5.73, 100]
        for i, score in enumerate(scores, 1):
            self.insert(i, score)
        actual = [r[0] for r in self.connection.execute(
            "SELECT lda_score_normalized FROM taxon_disease_associations ORDER BY id")]
        self.assertEqual(actual, sorted(actual))
        self.assertTrue(all(-1 <= n <= 1 for n in actual))
        self.assertEqual(actual[2:7], [-0.8, -0.5, 0, 0.5, 0.8])
        for score, scaled in zip(scores, actual):
            self.assertEqual(scaled, float(normalize_lda_score(score)))
            self.assertEqual(normalize_lda_score(-score), -normalize_lda_score(score))

    def test_missing_and_non_lda_scores_stay_null(self):
        self.insert(1, None)
        self.insert(2, -10.4403242349, "log2_fold_change")
        self.assertEqual(self.connection.execute("""SELECT effect_size, lda_score_normalized
            FROM taxon_disease_associations ORDER BY id""").fetchall(),
            [(None, None), (-10.4403242349, None)])

    def test_new_scores_and_type_changes_refresh_automatically(self):
        self.insert(1, 4)
        self.connection.execute("UPDATE taxon_disease_associations SET effect_size = -4")
        self.assertEqual(self.connection.execute("""SELECT lda_score_normalized
            FROM taxon_disease_associations""").fetchone()[0], -0.8)
        self.connection.execute("UPDATE taxon_disease_associations SET effect_type = 'reported'")
        self.assertIsNone(self.connection.execute("""SELECT lda_score_normalized
            FROM taxon_disease_associations""").fetchone()[0])

    def test_migration_and_repeat_preserve_scores(self):
        self.connection.execute("DROP VIEW v_taxon_disease_summary")
        self.connection.execute("ALTER TABLE taxon_disease_associations DROP COLUMN lda_score_normalized")
        self.insert(1, 4.31)
        self.assertTrue(install_normalized_lda(self.connection, sqlite=True))
        self.assertFalse(install_normalized_lda(self.connection, sqlite=True))
        self.assertEqual(self.connection.execute("""SELECT effect_size, lda_score_normalized
            FROM taxon_disease_associations""").fetchone(), (4.31, 0.812))
        fields = [r[0] for r in self.connection.execute("SELECT * FROM v_taxon_disease_summary").description]
        self.assertIn("lda_score_normalized", fields)


if __name__ == "__main__":
    unittest.main()
