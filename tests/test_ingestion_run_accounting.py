"""Every row a loader reads must be accounted for in the run's own row.

ingestion_runs is the only permanent record of what a load did, and it was the
last table validate did not mention. Its four counters are read, inserted,
updated and skipped, and the subtraction an auditor performs on them --
rows_read - (inserted + updated + skipped) -- returned a non-zero residual for
14 of the 92 runs.

enrich_taxonomy was one cause and is fixed: it counted its three drop paths in a
local dict that the CLI printed once and discarded, so a row NCBI could not
resolve was recorded nowhere at all. These tests pin the fixed accounting, path
by path, because the identity is only worth anything if every path maintains it.

The other cause is NOT a defect and must stay unasserted; the last test pins
that decision.
"""
import unittest

from gutdb.pipeline import (
    _MYSQL_ONLY_CHECKS,
    _validate_checks,
    enrich_taxonomy,
)


class Client:
    """Stands in for TaxonomyClient: a name map and a tax ID -> lineage map."""

    def __init__(self, names=None, lineages=None):
        self.names = names or {}
        self.lineages = lineages or {}

    def resolve_name(self, name, domain):
        return self.names.get(name)

    def fetch_lineages(self, batch):
        return {taxid: self.lineages[taxid]
                for taxid in batch if taxid in self.lineages}


class Connection:
    """One taxa SELECT, then UPDATEs, plus the two ingestion_runs writes.

    `recorded` holds what _finish_run put in the run's row, which is the thing
    under test: not what the loader counted in memory, but what survives.
    """

    def __init__(self, pending):
        self.pending = pending
        self.recorded = None
        self.updates = 0

    def cursor(self, **kwargs):
        return self

    @property
    def lastrowid(self):
        return 1

    def execute(self, sql, params=None):
        statement = " ".join(sql.split())
        if statement.startswith("UPDATE ingestion_runs"):
            _completed, status, read, inserted, updated, skipped, _err, _id = params
            self.recorded = {
                "status": status, "read": read, "inserted": inserted,
                "updated": updated, "skipped": skipped,
            }
        elif statement.startswith("UPDATE taxa"):
            self.updates += 1

    def fetchall(self):
        return self.pending

    def commit(self):
        pass

    def rollback(self):
        pass

    def close(self):
        pass


def run(pending, client, resolve_missing_ids=True):
    connection = Connection(pending)
    stats, counters = enrich_taxonomy(
        connection, client, limit=None, batch_size=150, sleep_seconds=0,
        resolve_missing_ids=resolve_missing_ids,
    )
    return connection.recorded, stats, counters


class RowAccountingTests(unittest.TestCase):
    """rows_read == inserted + updated + skipped, by every path that drops a row."""

    def assert_balances(self, recorded):
        self.assertEqual(
            recorded["read"],
            recorded["inserted"] + recorded["updated"] + recorded["skipped"],
            f"run row does not balance: {recorded}",
        )

    def test_every_row_enriched(self):
        recorded, _stats, _c = run(
            [(10, "Bacteroides fragilis", 817, "Bacteria")],
            Client(lineages={"817": {"superkingdom": "Bacteria",
                                     "phylum": "Bacteroidetes"}}),
        )
        self.assertEqual(recorded["read"], 1)
        self.assertEqual(recorded["updated"], 1)
        self.assertEqual(recorded["skipped"], 0)
        self.assert_balances(recorded)

    def test_name_that_cannot_be_resolved_is_skipped(self):
        """The 443-row path: read, looked up, no match, counted nowhere."""
        recorded, _stats, counters = run(
            [(11, "Unknown organism", None, "Bacteria")],
            Client(names={}),
        )
        self.assertEqual(counters["names_unresolved"], 1)
        self.assertEqual(recorded["read"], 1)
        self.assertEqual(recorded["updated"], 0)
        self.assertEqual(recorded["skipped"], 1)
        self.assert_balances(recorded)

    def test_row_with_no_taxid_and_no_resolution_is_skipped(self):
        """--no-name-resolution: this path incremented NO counter at all."""
        recorded, _stats, counters = run(
            [(12, "Unknown organism", None, "Bacteria")],
            Client(),
            resolve_missing_ids=False,
        )
        self.assertEqual(counters["no_taxid_unresolved"], 1)
        self.assertEqual(recorded["skipped"], 1)
        self.assert_balances(recorded)

    def test_missing_lineage_is_counted_per_row_not_per_taxid(self):
        """Two taxa rows share one tax ID, and NCBI returns nothing for it.

        The drop is counted once per tax ID in `counters` and must be counted
        twice in rows_skipped, because rows_read counted rows. Writing
        `stats.skipped += 1` here leaves the row short by one and this test is
        what notices.
        """
        recorded, _stats, counters = run(
            [(13, "Blautia wexlerae", 1,
              "Bacteria"), (14, "Blautia luti", 1, "Bacteria")],
            Client(lineages={}),
        )
        self.assertEqual(counters["no_lineage_returned"], 1)
        self.assertEqual(recorded["read"], 2)
        self.assertEqual(recorded["skipped"], 2)
        self.assert_balances(recorded)

    def test_all_paths_at_once(self):
        """One enriched, one unresolvable, two sharing a tax ID with no lineage."""
        recorded, _stats, _c = run(
            [
                (20, "Bacteroides fragilis", 817, "Bacteria"),
                (21, "Unknown organism", None, "Bacteria"),
                (22, "Blautia wexlerae", 1, "Bacteria"),
                (23, "Blautia luti", 1, "Bacteria"),
            ],
            Client(lineages={"817": {"superkingdom": "Bacteria",
                                     "phylum": "Bacteroidetes"}}),
        )
        self.assertEqual(recorded["read"], 4)
        self.assertEqual(recorded["updated"], 1)
        self.assertEqual(recorded["skipped"], 3)
        self.assert_balances(recorded)

    def test_a_failed_run_is_not_held_to_the_identity(self):
        """A run that dies mid-way records a snapshot, and that is correct.

        Run 16, NCBI_BioSample_demographics, is 150 rows short for exactly this
        reason. The counters assert nothing about failed runs.
        """
        class Exploding(Client):
            def fetch_lineages(self, batch):
                raise RuntimeError("NCBI unreachable")

        with self.assertRaises(RuntimeError):
            run([(30, "Bacteroides fragilis", 817, "Bacteria")], Exploding())


class RunCounterTests(unittest.TestCase):
    """The two counters validate now carries for this table."""

    def setUp(self):
        self.checks = _validate_checks()

    def test_both_counters_exist_and_are_portable(self):
        """The mirror has ingestion_runs too, so neither is MySQL-only."""
        for name in ("run_accounted_for_more_rows_than_it_read",
                     "run_abandoned_in_the_running_state"):
            self.assertIn(name, self.checks)
            self.assertNotIn(name, _MYSQL_ONLY_CHECKS)

    def test_abandoned_run_counter_cannot_fire_on_a_live_load(self):
        """It compares against the table's own clock, never NOW().

        An age cutoff would have had to be guessed, and a counter that fires
        while a load is legitimately in flight blocks a push for no reason.
        """
        sql = self.checks["run_abandoned_in_the_running_state"]
        self.assertIn("MAX(completed_at)", sql)
        self.assertNotIn("NOW()", sql.upper())
        self.assertNotIn("INTERVAL", sql.upper())

    def test_the_row_identity_is_not_asserted(self):
        """Deliberate, and this test is why it stays that way.

        rows_read = inserted + updated + skipped fails in 14 of 92 runs. Five
        are pre-fix enrich_taxonomy runs whose stored rows are history, one is a
        failed run's snapshot, and seven are MiMeDB, which counts rows_read in
        SOURCE ROWS while its outcomes are TAXA: 2,584 records collapse onto
        1,679 organisms and the residual 905 is what the collapse absorbed. A
        counter on the identity would fire 14 times forever and gate nothing.
        """
        for name, sql in self.checks.items():
            if not isinstance(name, tuple) and "ingestion_runs" in sql:
                flattened = " ".join(sql.split())
                self.assertNotIn("<> rows_read", flattened)
                self.assertNotIn("!= rows_read", flattened)
                self.assertNotIn("< rows_read", flattened)


if __name__ == "__main__":
    unittest.main()
