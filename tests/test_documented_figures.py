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

# The sample_attributes section. Its 362 appears twice in different sentences,
# so each is anchored on its own wording rather than on the number.
ATTRIBUTES = [
    ("attribute keys", r"keeps ([\d,]+) keys exactly as the source wrote them", "keys"),
    ("attribute keys, second mention", r"of the ([\d,]+) keys are unmapped", "keys"),
    ("attribute keys, third mention", r"looks_like_an_interval` for all ([\d,]+) keys", "keys"),
    ("unmapped keys", r"\*\*([\d,]+) of the [\d,]+ keys are unmapped", "unmapped"),
    ("keys mixing numbers and words", r"([\d,]+) keys mix numbers with words", "mixing"),
    ("PhysicalExercise numeric rows", r"given by ([\d,]+) of [\d,]+", "exercise_numeric"),
    ("PhysicalExercise total", r"given by [\d,]+ of ([\d,]+)", "exercise_total"),
]


# Figures written into gutdb/schema.sql comments. The same wrong count for
# mixed-value keys was in both files, which is the reason to check both: a
# figure repeated in two places goes stale in two places.
SCHEMA_FIGURES = [
    ("schema: attribute keys", r"-- ([\d,]+) keys, and [\d,]+ of them mix", "keys"),
    ("schema: keys mixing", r"-- [\d,]+ keys, and ([\d,]+) of them mix", "mixing"),
    ("schema: numeric Stage rows", r"mean of the ([\d,]+) rows that happen to hold", "stage_numeric"),
    ("schema: Stage total", r"stage out of ([\d,]+) whose values run", "stage_total"),
    ("schema: rank disagreements", r"-- ([\d,]+) rows disagree\.", "rank_disagree"),
    ("schema: genus layer coarser", r"rows disagree\. ([\d,]+) sit in the genus layer",
     "genus_coarser"),
    ("schema: families in the genus layer", r"coarser -- ([\d,]+) families", "genus_family"),
    ("schema: species layer finer", r"and ([\d,]+) sit\n\s*-- in the species layer",
     "species_finer"),
    ("schema: qc_status empty", r"which is why ([\d,]+) rows are\n\s*-- empty", "qc_null"),
    ("schema: negative metabolite levels",
     r"z-transformed and hold ([\d,]+) negative levels", "levels_negative"),
    ("schema: negative rows for log1p",
     r"NaN on the ([\d,]+) negative rows", "levels_negative"),
    ("schema: non-negative rows in centred cohorts",
     r"log of a log on the other ([\d,]+) in those", "centred_nonnegative"),
    ("schema: levels in the centred cohorts",
     r"([\d,]+) levels between them", "centred_total"),
]


def attribute_metrics(connection) -> dict[str, int]:
    cursor = connection.cursor()
    try:
        def one(sql, args=()):
            cursor.execute(sql, args)
            return int(cursor.fetchone()[0])
        return dict(
            keys=one("SELECT COUNT(*) FROM v_sample_attribute_keys"),
            # The project's own flag, not a re-derivation of it.
            mixing=one("SELECT COUNT(*) FROM v_sample_attribute_keys "
                       "WHERE mixes_numbers_and_words = 1"),
            unmapped=one("""SELECT COUNT(*) FROM (
                SELECT DISTINCT a.attribute FROM sample_attributes a
                WHERE NOT EXISTS (SELECT 1 FROM sample_attribute_aliases al
                                  WHERE al.attribute = a.attribute)) d"""),
            stage_total=one("SELECT COUNT(*) FROM sample_attributes WHERE attribute = 'Stage'"),
            stage_numeric=one("SELECT COUNT(*) FROM sample_attributes "
                              "WHERE attribute = 'Stage' AND value_numeric IS NOT NULL"),
            rank_disagree=one("""SELECT COUNT(*) FROM taxa
                WHERE ncbi_rank IS NOT NULL
                  AND ((taxonomic_rank = 'genus'
                        AND ncbi_rank NOT IN ('genus', 'subgenus', 'species',
                                              'subspecies', 'strain', 'no rank', 'clade'))
                    OR (taxonomic_rank = 'species'
                        AND ncbi_rank NOT IN ('species', 'no rank', 'clade')))"""),
            exercise_numeric=one("SELECT COUNT(*) FROM sample_attributes "
                                 "WHERE attribute = 'PhysicalExercise_LastWeek.h.' "
                                 "AND value_numeric IS NOT NULL"),
            exercise_total=one("SELECT COUNT(*) FROM sample_attributes "
                               "WHERE attribute = 'PhysicalExercise_LastWeek.h.'"),
            # The rank vocabularies are the ones v_abundance_genus and
            # v_abundance_species use, so these cannot pass by agreeing with a
            # definition the views do not apply.
            genus_coarser=one("""SELECT COUNT(*) FROM taxa
                WHERE taxonomic_rank = 'genus' AND ncbi_rank IS NOT NULL
                  AND ncbi_rank NOT IN ('genus', 'subgenus', 'species',
                                        'subspecies', 'strain', 'no rank', 'clade')"""),
            genus_family=one("SELECT COUNT(*) FROM taxa WHERE taxonomic_rank = 'genus' "
                             "AND ncbi_rank = 'family'"),
            species_finer=one("""SELECT COUNT(*) FROM taxa
                WHERE taxonomic_rank = 'species' AND ncbi_rank IS NOT NULL
                  AND ncbi_rank NOT IN ('species', 'no rank', 'clade')"""),
            qc_null=one("SELECT COUNT(*) FROM samples WHERE qc_status IS NULL"),
            levels_negative=one("SELECT COUNT(*) FROM sample_metabolite_levels "
                                "WHERE level < 0"),
            centred_total=one("""SELECT COUNT(*) FROM sample_metabolite_levels m
                JOIN samples s ON s.id = m.sample_id
                JOIN studies st ON st.id = s.study_id
                WHERE st.metabolite_scale = 'centred'"""),
            centred_nonnegative=one("""SELECT COUNT(*) FROM sample_metabolite_levels m
                JOIN samples s ON s.id = m.sample_id
                JOIN studies st ON st.id = s.study_id
                WHERE st.metabolite_scale = 'centred' AND m.level >= 0"""))
    finally:
        cursor.close()


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
        # exists() is true for a file iCloud has evicted, and connect()
        # succeeds lazily on one, so the hang lands on the first query -- which
        # is how this suite once stalled past ten minutes. require_local turns
        # that into the skip below.
        from gutdb.config import require_local
        return sqlite3.connect(require_local(path), timeout=60)
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
        cls.schema = (ROOT / "gutdb" / "schema.sql").read_text()
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

    def test_the_attribute_figures_match_the_database(self):
        measured = attribute_metrics(self.connection)
        for label, pattern, key in ATTRIBUTES:
            with self.subTest(figure=label):
                found = re.search(pattern, self.text)
                self.assertIsNotNone(found, f"no anchor for {label!r}")
                documented = int(found.group(1).replace(",", ""))
                self.assertEqual(documented, measured[key],
                                 f"README says {documented:,} for {label}, "
                                 f"database says {measured[key]:,}")

    def test_the_schema_comment_figures_match_the_database(self):
        measured = attribute_metrics(self.connection)
        for label, pattern, key in SCHEMA_FIGURES:
            with self.subTest(figure=label):
                found = re.search(pattern, self.schema)
                self.assertIsNotNone(found, f"no anchor for {label!r}")
                documented = int(found.group(1).replace(",", ""))
                self.assertEqual(documented, measured[key],
                                 f"schema.sql says {documented:,} for {label}, "
                                 f"database says {measured[key]:,}")

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
        for label, pattern, _ in FIGURES + MIRROR_TABLE + ATTRIBUTES:
            with self.subTest(figure=label):
                self.assertEqual(len(re.findall(pattern, self.text)), 1,
                                 f"anchor for {label!r} no longer matches exactly once")
        for label, pattern, _ in SCHEMA_FIGURES:
            with self.subTest(figure=label):
                self.assertEqual(len(re.findall(pattern, self.schema)), 1,
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
