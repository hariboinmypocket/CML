"""The README's current-state figures must still be true.

They rotted silently. `13,328` was stated as the modelling-ready cohort size in
three places; the CKD withdrawal had already taken it to 12,896 and the
replicate-representative fix then took it to 13,008, so the project's own record
was wrong through two changes and nothing noticed. One figure did not even
reconcile with its own sentence -- "389 subjects contribute 1,340 ... so 971 are
repeat visits", where 1,340 - 389 is 951.

No expected value is written here. Each figure is LOCATED in the README by its
surrounding wording and compared against the database, so there is nothing to
keep in step. If a rewording breaks the pattern the test fails loudly, which is
correct: it means the anchor needs re-pointing rather than silently passing.

Skipped when no database is reachable, so the suite still runs anywhere.
"""
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

# (label, pattern with one capturing group, metric key). The metrics are
# computed in four queries rather than one per figure, because
# v_ml_ready_samples costs about four seconds a scan -- almost all of it a
# pre-existing rollup over 3M abundance rows -- and scanning it once per figure
# made this the slowest thing in the suite by an order of magnitude.
FIGURES = [
    ("modelling-ready cohort size",
     r"the modelling-ready cohort is ([\d,]+) samples", "deduplicated"),
    ("cohort before the subject collapse",
     r"taking the cohort from ([\d,]+) to [\d,]+\.", "deduplicated"),
    ("cohort after the subject collapse",
     r"taking the cohort from [\d,]+ to ([\d,]+)\.", "one_per_subject"),
    ("rows with no subject_id in the cohort",
     r"([\d,]+) of the [\d,]+\s+deduplicated rows", "no_subject"),
    ("rows the subject flag constrains",
     r"genuinely constrains the ([\d,]+) MullerMM", "with_subject"),
    ("repeat visits in the cohort",
     r"so \*\*([\d,]+) are repeat visits", "repeat_visits"),
    ("subjects with more than one row",
     r"([\d,]+) subjects contribute [\d,]+ of the cohort's rows", "multi_subjects"),
    ("rows those subjects contribute",
     r"[\d,]+ subjects contribute ([\d,]+) of the cohort's rows", "multi_subject_rows"),
    ("metabolites", r"took metabolites to \*\*([\d,]+)\*\*", "metabolites"),
    ("metabolite levels", r"levels to\s+\*\*([\d,]+)\*\*", "metabolite_levels"),
]

# The mirror coverage table. Both columns are checked, because the point of the
# table is the RATIO between them, and a figure that drifted on either side
# would misstate it.
MIRROR_TABLE = [
    ("samples", r"\| samples \| ([\d,]+) \| [\d,]+ \|", "mysql_samples"),
    ("samples, mirror", r"\| samples \| [\d,]+ \| ([\d,]+) \|", "mirror_samples"),
    ("abundance rows", r"\| abundance rows \| ([\d,]+) \| [\d,]+ \|", "mysql_abundance"),
    ("abundance rows, mirror", r"\| abundance rows \| [\d,]+ \| ([\d,]+) \|", "mirror_abundance"),
    ("taxa", r"\| taxa \| ([\d,]+) \| [\d,]+ \|", "mysql_taxa"),
    ("taxa, mirror", r"\| taxa \| [\d,]+ \| ([\d,]+) \|", "mirror_taxa"),
    ("genus layer", r"genus layer \| ([\d,]+) \| [\d,]+ \|", "mysql_genus"),
    ("genus layer, mirror", r"genus layer \| [\d,]+ \| ([\d,]+) \|", "mirror_genus"),
    ("species layer", r"species layer \| ([\d,]+) \| [\d,]+ \|", "mysql_species"),
    ("species layer, mirror", r"species layer \| [\d,]+ \| ([\d,]+) \|", "mirror_species"),
]

# Both layers in one scan. Counting them separately meant two passes over 3M
# abundance rows per store, which doubled the cost of this file for nothing.
LAYERS = """
    SELECT COUNT(DISTINCT CASE WHEN t.taxonomic_rank = 'genus'
                               THEN a.sample_id END),
           COUNT(DISTINCT CASE WHEN t.taxonomic_rank = 'species'
                               THEN a.sample_id END)
    FROM sample_taxon_abundances a
    JOIN taxa t ON t.id = a.taxon_id"""


def coverage(mysql, mirror) -> dict[str, int]:
    cursor = mysql.cursor()
    try:
        def one(sql):
            cursor.execute(sql)
            return int(cursor.fetchone()[0])
        measured = dict(
            mysql_samples=one("SELECT COUNT(*) FROM samples"),
            mysql_abundance=one("SELECT COUNT(*) FROM sample_taxon_abundances"),
            mysql_taxa=one("SELECT COUNT(*) FROM taxa"))
        cursor.execute(LAYERS)
        genus, species = cursor.fetchone()
        measured.update(mysql_genus=int(genus), mysql_species=int(species))
    finally:
        cursor.close()

    def o(sql):
        return int(mirror.execute(sql).fetchone()[0])
    measured.update(
        mirror_samples=o("SELECT COUNT(*) FROM samples"),
        mirror_abundance=o("SELECT COUNT(*) FROM sample_taxon_abundances"),
        mirror_taxa=o("SELECT COUNT(*) FROM taxa"))
    genus, species = mirror.execute(LAYERS).fetchone()
    measured.update(mirror_genus=int(genus), mirror_species=int(species))
    return measured


def metrics(connection) -> dict[str, int]:
    cursor = connection.cursor()
    try:
        cursor.execute("""
            SELECT SUM(ml_ready_deduplicated = 1),
                   SUM(ml_ready_one_per_subject = 1),
                   SUM(ml_ready_deduplicated = 1 AND subject_id IS NULL),
                   SUM(ml_ready_deduplicated = 1 AND subject_id IS NOT NULL)
            FROM v_ml_ready_samples""")
        deduplicated, one_per_subject, no_subject, with_subject = cursor.fetchone()
        cursor.execute("""
            SELECT COUNT(*), COALESCE(SUM(c), 0), COALESCE(SUM(c - 1), 0) FROM (
                SELECT COUNT(*) c FROM samples s
                JOIN v_ml_ready_samples v ON v.sample_id = s.id
                WHERE v.ml_ready_deduplicated = 1 AND s.subject_id IS NOT NULL
                GROUP BY s.study_id, s.subject_id HAVING c > 1) d""")
        multi_subjects, multi_subject_rows, repeat_visits = cursor.fetchone()
        cursor.execute("SELECT COUNT(*) FROM metabolites")
        metabolites = cursor.fetchone()[0]
        cursor.execute("SELECT COUNT(*) FROM sample_metabolite_levels")
        metabolite_levels = cursor.fetchone()[0]
    finally:
        cursor.close()
    return {k: int(v) for k, v in dict(
        deduplicated=deduplicated, one_per_subject=one_per_subject,
        no_subject=no_subject, with_subject=with_subject,
        multi_subjects=multi_subjects, multi_subject_rows=multi_subject_rows,
        repeat_visits=repeat_visits, metabolites=metabolites,
        metabolite_levels=metabolite_levels).items()}


def mirror_connection():
    import sqlite3
    path = ROOT / "gut_microbiome.sqlite"
    if not path.exists():
        return None
    try:
        return sqlite3.connect(path, timeout=60)
    except Exception:
        return None


def connection():
    try:
        from gutdb.config import Settings
        from gutdb.db import connect
        return connect(Settings.from_env())
    except Exception:
        return None


class DocumentedFigureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = (ROOT / "README.md").read_text()
        cls.connection = connection()
        if cls.connection is None:
            raise unittest.SkipTest("no database reachable")
        cls.mirror = mirror_connection()

    @classmethod
    def tearDownClass(cls):
        if getattr(cls, "connection", None) is not None:
            cls.connection.close()
        if getattr(cls, "mirror", None) is not None:
            cls.mirror.close()

    def test_the_mirror_coverage_table_matches_both_stores(self):
        """The table's point is the ratio, so both columns are checked."""
        if self.mirror is None:
            self.skipTest("no mirror file")
        measured = coverage(self.connection, self.mirror)
        for label, pattern, key in MIRROR_TABLE:
            with self.subTest(figure=label):
                found = re.search(pattern, self.text)
                self.assertIsNotNone(found, f"no anchor for {label!r}")
                documented = int(found.group(1).replace(",", ""))
                self.assertEqual(documented, measured[key],
                                 f"README says {documented:,} for {label}, "
                                 f"stores say {measured[key]:,}")

    def test_the_mirror_is_a_strict_subset_of_mysql(self):
        """It may hold less. It must never hold a run MySQL does not."""
        if self.mirror is None:
            self.skipTest("no mirror file")
        cursor = self.connection.cursor()
        try:
            cursor.execute("SELECT run_accession FROM samples")
            live = {r[0] for r in cursor.fetchall()}
        finally:
            cursor.close()
        held = {r[0] for r in self.mirror.execute("SELECT run_accession FROM samples")}
        self.assertEqual(sorted(held - live)[:10], [],
                         "mirror holds run accessions absent from MySQL")

    def test_every_anchor_still_matches_exactly_once(self):
        """A reworded sentence must fail here, not quietly stop being checked."""
        for label, pattern, _ in FIGURES + MIRROR_TABLE:
            with self.subTest(figure=label):
                self.assertEqual(len(re.findall(pattern, self.text)), 1,
                                 f"anchor for {label!r} no longer matches exactly once")

    def test_documented_figures_match_the_database(self):
        measured = metrics(self.connection)
        for label, pattern, key in FIGURES:
            with self.subTest(figure=label):
                found = re.search(pattern, self.text)
                self.assertIsNotNone(found, f"no anchor for {label!r}")
                documented = int(found.group(1).replace(",", ""))
                self.assertEqual(documented, measured[key],
                                 f"README says {documented:,} for {label}, "
                                 f"database says {measured[key]:,}")


if __name__ == "__main__":
    unittest.main()
