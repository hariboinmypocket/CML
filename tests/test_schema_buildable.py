"""The schemas must build on a machine that has never seen this database.

This is the one class of defect the validate counters cannot reach. They run
against an already-built database, where every view exists, so a view defined
before the view it selects from passes every check and still makes
`gutdb-cli init` fail on a clean machine. That is exactly what happened:
v_attribute_concept_coverage was defined ~300 lines above v_ml_ready_samples
and MySQL refused it with "Table '<db>.v_ml_ready_samples' doesn't exist".

The MySQL schema is checked by reading it rather than by applying it, so this
needs no server and no credentials. The SQLite schema is actually built,
because it can be.
"""
import re
import sqlite3
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFINITION = re.compile(r"CREATE (?:OR REPLACE )?VIEW (?:IF NOT EXISTS )?(\w+)")


def view_positions(text):
    return {m.group(1): m.start() for m in DEFINITION.finditer(text)}


def forward_references(text):
    """Views that select from a view defined later in the same file."""
    views = view_positions(text)
    found = []
    for name, start in views.items():
        following = [p for p in views.values() if p > start]
        body = text[start:min(following) if following else len(text)]
        for other, other_start in views.items():
            if other == name or other_start <= start:
                continue
            if re.search(r"\b(?:FROM|JOIN)\s+" + re.escape(other) + r"\b", body):
                line = text[:start].count("\n") + 1
                other_line = text[:other_start].count("\n") + 1
                found.append(f"{name} (line {line}) selects from "
                             f"{other} (defined later, line {other_line})")
    return found


class SchemaOrderingTests(unittest.TestCase):
    def test_no_view_selects_from_a_view_defined_later(self):
        for name in ("schema.sql", "schema_sqlite.sql"):
            text = (ROOT / "gutdb" / name).read_text()
            with self.subTest(schema=name):
                self.assertEqual(forward_references(text), [])

    def test_both_schemas_actually_define_views(self):
        """A regex that stopped matching would make the check above vacuous."""
        for name, least in (("schema.sql", 10), ("schema_sqlite.sql", 6)):
            text = (ROOT / "gutdb" / name).read_text()
            with self.subTest(schema=name):
                self.assertGreaterEqual(len(view_positions(text)), least)

    def test_the_check_catches_a_forward_reference(self):
        """Pin the detector itself, so a silent regex failure is not a pass."""
        bad = ("CREATE OR REPLACE VIEW v_first AS SELECT * FROM v_second;\n"
               "CREATE OR REPLACE VIEW v_second AS SELECT 1;\n")
        self.assertEqual(len(forward_references(bad)), 1)
        good = ("CREATE OR REPLACE VIEW v_second AS SELECT 1;\n"
                "CREATE OR REPLACE VIEW v_first AS SELECT * FROM v_second;\n")
        self.assertEqual(forward_references(good), [])


class SqliteFreshBuildTests(unittest.TestCase):
    def test_a_fresh_mirror_builds_and_every_view_resolves(self):
        schema = (ROOT / "gutdb" / "schema_sqlite.sql").read_text()
        with tempfile.TemporaryDirectory() as directory:
            connection = sqlite3.connect(Path(directory) / "fresh.sqlite")
            try:
                connection.executescript(schema)
                views = [r[0] for r in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type='view'")]
                self.assertTrue(views)
                for view in views:
                    # SQLite accepts a view over a missing table at creation
                    # time and only fails when it is queried.
                    with self.subTest(view=view):
                        connection.execute(f"SELECT * FROM {view} LIMIT 0")
            finally:
                connection.close()


if __name__ == "__main__":
    unittest.main()
