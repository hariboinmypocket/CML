"""Demographic reconciliation must preserve sample identity and missingness."""
import unittest

from gutdb.pipeline import GMREPO_SAMPLE_FIELDS, prepare_gmrepo_sample_export


class Connection:
    """Two queries now run here, so the mock has to tell them apart."""

    def __init__(self, rows=(), withdrawals=()):
        self.rows = list(rows)
        self.withdrawals = [{"run_accession": run, "field": field}
                            for run, field in withdrawals]
        self.pending = None

    def cursor(self, **kwargs):
        return self

    def execute(self, sql, *args):
        self.pending = (self.withdrawals if "demographic_withdrawals" in sql
                        else self.rows)

    def fetchall(self):
        return self.pending

    def close(self):
        pass


def raw(**updates):
    row = dict.fromkeys(GMREPO_SAMPLE_FIELDS, "")
    row.update(project_id="PROJECT", run_id="RUN", sample_id="SAMPLE", qc_status="0")
    row.update(updates)
    return row


def stored(**updates):
    row = dict(run_accession="RUN", project_accession="PROJECT", gmrepo_sample_id="SAMPLE",
               sex=None, age_years=None, bmi=None, disease_name=None, mesh_id=None)
    row.update(updates)
    return row


class SampleExportTests(unittest.TestCase):
    def test_fills_demographics_from_exact_record_and_keeps_qc(self):
        connection = Connection([stored(sex="female", age_years=42.5, bmi=23.1)])
        rows, audit = prepare_gmrepo_sample_export(connection, [raw()])
        self.assertEqual((rows[0]["sex"], rows[0]["age_years"], rows[0]["bmi"]),
                         ("female", "42.5", "23.1"))
        self.assertEqual(rows[0]["qc_status"], "0")
        self.assertEqual({a["field"] for a in audit}, {"sex", "age_years", "bmi"})
        self.assertEqual(prepare_gmrepo_sample_export(connection, rows)[1], [])

    def test_withdrawn_cells_are_not_refilled_from_the_source(self):
        """A cleared cell and a never-measured cell look identical here.

        Without the ledger the fallback to the raw source value re-imports
        exactly what a withdrawal removed, and the audit cannot tell: the value
        matches the source, so nothing reads as a change.
        """
        connection = Connection(
            [stored()], withdrawals=[("RUN", "sex"), ("RUN", "age_years"),
                                     ("RUN", "bmi")])
        rows, audit = prepare_gmrepo_sample_export(
            connection, [raw(sex="female", age_years="54", bmi="23.4")])
        self.assertEqual([rows[0][f] for f in ("sex", "age_years", "bmi")], ["", "", ""])
        self.assertEqual({a["source"] for a in audit},
                         {"withheld_by_demographic_withdrawals"})

    def test_withdrawal_is_per_field_and_per_run(self):
        connection = Connection([stored(), stored(run_accession="OTHER")],
                                withdrawals=[("RUN", "bmi")])
        rows, _ = prepare_gmrepo_sample_export(
            connection, [raw(sex="female", age_years="54", bmi="23.4"),
                         raw(run_id="OTHER", sample_id="SAMPLE", bmi="31.0")])
        self.assertEqual([rows[0][f] for f in ("sex", "age_years", "bmi")],
                         ["female", "54.0", ""])
        self.assertEqual(rows[1]["bmi"], "31.0")

    def test_a_withdrawn_cell_still_stored_in_the_database_is_not_recovered(self):
        connection = Connection([stored(sex="female", age_years=54.0)],
                                withdrawals=[("RUN", "sex"), ("RUN", "age_years")])
        rows, _ = prepare_gmrepo_sample_export(connection, [raw()])
        self.assertEqual([rows[0][f] for f in ("sex", "age_years")], ["", ""])

    def test_identity_conflicts_abort(self):
        for change in ({"project_accession": "OTHER"}, {"gmrepo_sample_id": "OTHER"}):
            with self.subTest(change=change), self.assertRaisesRegex(ValueError, "identity mismatch"):
                prepare_gmrepo_sample_export(Connection([stored(**change)]), [raw()])

    def test_unreported_values_stay_blank_and_invalid_sentinels_are_removed(self):
        rows, _ = prepare_gmrepo_sample_export(Connection(), [raw(age_years="-88", bmi="0")])
        self.assertEqual([rows[0][f] for f in ("sex", "age_years", "bmi")], ["", "", ""])
        rows, _ = prepare_gmrepo_sample_export(Connection(), [raw(age_years="0")])
        self.assertEqual(rows[0]["age_years"], "0.0")

    def test_flattened_phenotype_survives_without_raw_columns(self):
        row = raw(phenotypes='[{"term":"Health","disease":"D006262"}]',
                  nr_phenotypes="1", longitude="1", latitude="2", data_quality="curated")
        rows, _ = prepare_gmrepo_sample_export(Connection(), [row])
        self.assertEqual((rows[0]["disease_name"], rows[0]["mesh_id"]), ("Health", "D006262"))
        self.assertEqual(set(rows[0]), set(GMREPO_SAMPLE_FIELDS))
        self.assertFalse(set(rows[0]) & {"phenotypes", "nr_phenotypes", "longitude", "latitude", "data_quality"})

    def test_ambiguous_phenotypes_are_not_guessed(self):
        row = raw(phenotypes='[{"term":"Health","disease":"D006262"},'
                             '{"term":"IBS","disease":"D043183"}]')
        rows, _ = prepare_gmrepo_sample_export(Connection(), [row])
        self.assertEqual((rows[0]["disease_name"], rows[0]["mesh_id"]), ("", ""))

    def test_curated_disease_and_mesh_pair_is_preserved(self):
        connection = Connection([stored(disease_name="Health", mesh_id="D006262")])
        rows, audit = prepare_gmrepo_sample_export(connection, [raw(disease_name="IBS", mesh_id="D043183")])
        self.assertEqual((rows[0]["disease_name"], rows[0]["mesh_id"]), ("Health", "D006262"))
        self.assertEqual({a["field"] for a in audit}, {"disease_name", "mesh_id"})


if __name__ == "__main__":
    unittest.main()
