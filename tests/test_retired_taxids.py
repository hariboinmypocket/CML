"""No row may point at a retired tax ID that could simply be repointed.

38 taxa carried a tax ID absent from nodes.dmp, with 30,950 abundance rows and
186 associations hanging off them. A retired ID joins to nothing through NCBI,
and nothing noticed: `ncbi_rank` goes NULL for such a row, because the rank
backfill cannot answer for an ID the dump has never heard of, and a NULL there
looked like missing data rather than a broken join.

`repoint_merged_taxids.py` fixed the 25 whose replacement was unoccupied. The
13 that remain are organism merges, not ID corrections -- Allisonella merges
onto Dialister, Anaeroglobus and Colibacter both onto Megasphaera, Raoultella
onto Klebsiella -- and the database holds both names as separate rows, so their
per-(sample, rank) shares would have to be summed or the sum to 1 breaks.

So this does not assert a count, and keeps no list of exceptions to go stale.
It asserts the shape: every row still pointing at a retired ID must be one the
repointer refused for a stated reason. A newly retired ID with a free target
fails here, which is the case a person can fix by running the script.

This cannot be a validate counter. Deciding it needs merged.dmp, and no SQL
reaches a file. Skipped when taxdump/ is absent, which it is on a fresh clone --
the dumps are gitignored.
"""
import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

TAXDUMP = ROOT / "taxdump"
MERGED = TAXDUMP / "merged.dmp"
NODES = TAXDUMP / "nodes.dmp"


def load_repointer():
    path = ROOT / "scripts" / "repoint_merged_taxids.py"
    spec = importlib.util.spec_from_file_location("repoint_merged_taxids", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class RetiredTaxIdTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for required in (MERGED, NODES):
            if not required.exists():
                raise unittest.SkipTest(f"{required.name} not present (gitignored)")
        cls.mod = load_repointer()
        try:
            from gutdb.config import Settings
            from gutdb.db import connect
            cls.mysql = connect(Settings.from_env())
        except Exception:
            raise unittest.SkipTest("no database reachable")
        import sqlite3
        cls.lite = None
        mirror = ROOT / "gut_microbiome.sqlite"
        if mirror.exists():
            try:
                from gutdb.config import require_local
                cls.lite = sqlite3.connect(require_local(mirror), timeout=60)
            except Exception:
                cls.lite = None
        # Read once: nodes.dmp is 232MB and both tests work off the same plan.
        cls.merged = cls.mod.read_merged(MERGED)
        cls.ranks = cls.mod.read_ranks(NODES)
        cls.plans = {}
        cls.plans["mysql"] = cls.mod.plan(
            cls.mod.Store("mysql", cls.mysql, "%s"), cls.merged, cls.ranks)
        if cls.lite is not None:
            cls.plans["sqlite"] = cls.mod.plan(
                cls.mod.Store("sqlite", cls.lite, "?"), cls.merged, cls.ranks)

    @classmethod
    def tearDownClass(cls):
        for handle in ("mysql", "lite"):
            if getattr(cls, handle, None) is not None:
                getattr(cls, handle).close()

    def test_the_dumps_were_actually_read(self):
        """Without this the two tests below pass by planning against nothing."""
        self.assertGreater(len(self.merged), 10_000)
        self.assertGreater(len(self.ranks), 1_000_000)

    def test_nothing_is_merely_waiting_to_be_repointed(self):
        # The same cross-store veto the script applies, from the script, so this
        # cannot drift from it: a key refused by any store is refused for all.
        blocked = self.mod.blocked_keys(self.plans)
        pending = []
        for label, (repoints, _) in self.plans.items():
            for key, row_id, name, old, new, rank in repoints:
                if key in blocked:
                    continue
                pending.append("%s: %s (%s -> %s)" % (label, name, old, new))
        self.assertEqual(
            [], pending,
            "%d row(s) point at a retired tax ID whose replacement is free. Run "
            "scripts/repoint_merged_taxids.py --apply:\n  %s"
            % (len(pending), "\n  ".join(pending)))

    def test_every_remaining_retired_id_has_a_stated_reason(self):
        """A refusal must explain itself, so none is a silent leftover."""
        for label, (_, refusals) in self.plans.items():
            dead = self.mod.Store(
                label, self.mysql if label == "mysql" else self.lite,
                "%s" if label == "mysql" else "?").rows(self.mod.SELECT_DEAD)
            explained = {key for key, _ in refusals} | set(
                self.mod.blocked_keys(self.plans))
            for row in dead:
                key = (row[4], row[5])
                with self.subTest(store=label, taxon=row[1]):
                    self.assertIn(
                        key, explained,
                        "%r still points at retired %s with no reason recorded"
                        % (row[1], row[2]))
                    reasons = [m for k, m in refusals if k == key]
                    reasons += [self.mod.blocked_keys(self.plans).get(key, "")]
                    self.assertTrue(
                        any(reasons), "%r is refused with an empty reason" % (row[1],))


if __name__ == "__main__":
    unittest.main()
