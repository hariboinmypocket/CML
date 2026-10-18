from decimal import Decimal
import unittest

from gutdb.transform import reconcile_demographic_cells, verified_demographic_link
from scripts.repair_ckd_demographics import assert_after, apply_sql, fetch_targets, plan
import sqlite3


class CKDRepairTests(unittest.TestCase):
    def test_all_three_verification_checks_required(self):
        record = dict(published_values_match_audit=True, ena_library_to_run_verified=True,
                      exact_source_identifier_in_ena=True)
        self.assertTrue(verified_demographic_link(record))
        for key in record:
            self.assertFalse(verified_demographic_link({**record, key: False}))
            self.assertFalse(verified_demographic_link({**record, key: "true"}))
        self.assertFalse(verified_demographic_link({"previous_join_basis": "documented"}))

    def test_clear_only_matching_imported_cells(self):
        source = dict(sex="female", age_years="40", bmi="22.5")
        current = dict(sex="Female", age_years=Decimal("40"), bmi=25)
        edits, conflicts = reconcile_demographic_cells(current, source, False)
        self.assertEqual(edits, dict(sex=None, age_years=None))
        self.assertEqual(conflicts, ["bmi"])
        current.update(edits)
        self.assertEqual(reconcile_demographic_cells(current, source, False)[0], {})

    def test_verified_fill_does_not_invent_bmi_or_overwrite(self):
        edits, conflicts = reconcile_demographic_cells(
            dict(sex="male", age_years=None, bmi=None),
            dict(sex="female", age_years="52", bmi=""), True)
        self.assertEqual(edits, {"age_years": "52"})
        self.assertEqual(conflicts, ["sex"])

    def test_both_store_comparison_accepts_numeric_sql_types(self):
        before = {"R": dict(id=1, sex=None, age_years=None, bmi=None)}
        after = {"R": dict(id=1, sex="male", age_years=Decimal("30"), bmi=None)}
        assert_after(before, after, {"R": dict(sex="male", age_years="30")})
        after["R"]["id"] = 99
        with self.assertRaises(RuntimeError):
            assert_after(before, after, {"R": dict(sex="male", age_years="30")})

    def test_sqlite_preserves_samples_and_local_identity(self):
        db = sqlite3.connect(":memory:")
        db.row_factory = sqlite3.Row
        db.executescript("""
            CREATE TABLE studies (id INTEGER PRIMARY KEY, project_accession TEXT);
            CREATE TABLE samples (id INTEGER PRIMARY KEY, study_id INTEGER,
                run_accession TEXT UNIQUE, sex TEXT, age_years REAL, bmi REAL);
            INSERT INTO studies VALUES(88, 'PRJNA562327');
            INSERT INTO samples VALUES(1904, 88, 'bad', 'female', 40, 22);
            INSERT INTO samples VALUES(987, 88, 'good', NULL, NULL, NULL);
        """)
        sources = {run: dict(sex="female", age_years="40", bmi="22") for run in ("bad", "good")}
        evidence = {run: dict(published_values_match_audit=True, ena_library_to_run_verified=True,
                             exact_source_identifier_in_ena=(run == "good")) for run in sources}
        before = fetch_targets(db, False, sources)
        edits, conflicts, counts = plan(before, sources, evidence)
        self.assertFalse(conflicts)
        self.assertEqual(sum(counts.values()), 6)
        apply_sql(db, False, edits)
        after = fetch_targets(db, False, sources)
        assert_after(before, after, edits)
        self.assertEqual(after["bad"]["id"], 1904)
        self.assertEqual(after["good"]["age_years"], 40)
        self.assertEqual(plan(after, sources, evidence)[0], {})
        self.assertEqual(db.execute("SELECT COUNT(*) FROM samples").fetchone()[0], 2)
        db.close()


if __name__ == "__main__":
    unittest.main()
