"""A row with a tax ID must not leave a lineage cell the dump could fill.

Only emptiness is asserted. Stored lineage names that disagree with the dump are
a separate and already-decided matter: scripts/reconcile_lineage_ranks.py
separates a rename (the stored name is an older spelling of the same node --
`Actinobacteria` and `Actinomycetes` are both class 1760) from a
reclassification, and deliberately leaves the renames alone, because adopting
them would rewrite hundreds of rows to a name the rest of the table does not
use. Nothing here reopens that. A cell with no value has no competing name, so
filling it decides nothing.

100 cells were empty when this was written -- 49 in MySQL, 51 in the mirror --
and some only became fillable once the 114 one-sided tax IDs were filled, since
a row with no tax ID has no lineage to look up.

This cannot be a validate counter: deciding it needs the dump, and no SQL
reaches a file. Skipped when taxdump/ is absent, which it is on a fresh clone.
"""
import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

TAXDUMP = ROOT / "taxdump"


def load_filler():
    path = ROOT / "scripts" / "fill_empty_lineage_from_dump.py"
    spec = importlib.util.spec_from_file_location("fill_empty_lineage_from_dump", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class LineageCompletenessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for required in ("names.dmp", "nodes.dmp"):
            if not (TAXDUMP / required).exists():
                raise unittest.SkipTest(f"{required} not present (gitignored)")
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
        cls.mod = load_filler()
        # The dump is read once: nodes.dmp alone is 232MB.
        cls.dump = cls.mod.load_taxdump(str(TAXDUMP))
        cls.plans = {"mysql": cls.mod.plan(cls.mysql, False, cls.dump)}
        if cls.lite is not None:
            cls.plans["sqlite"] = cls.mod.plan(cls.lite, True, cls.dump)

    @classmethod
    def tearDownClass(cls):
        for handle in ("mysql", "lite"):
            if getattr(cls, handle, None) is not None:
                getattr(cls, handle).close()

    def test_the_dump_was_actually_read(self):
        """Without this the test below passes by planning against nothing."""
        self.assertGreater(len(self.dump.rank_of), 1_000_000)

    def test_no_empty_lineage_cell_is_waiting_on_the_dump(self):
        pending = []
        for label, (fills, _refusals) in self.plans.items():
            pending.extend("%s: %s.%s" % (label, name, column)
                           for _row_id, name, column, _value in fills)
        self.assertEqual(
            [], pending,
            "%d lineage cell(s) are empty where the dump has an answer. Run "
            "scripts/fill_empty_lineage_from_dump.py --apply:\n  %s"
            % (len(pending), "\n  ".join(pending[:20])))

    def test_a_superkingdom_that_is_not_top_level_is_refused(self):
        """`Taxdump.lineage` resolves superkingdom through kingdom as a
        fallback, so a virus row can be handed `Heunggongvirae` for a cell whose
        correct value is `Viruses`. The filler must refuse that, not write it.
        """
        for label, (fills, _refusals) in self.plans.items():
            for _row_id, name, column, value in fills:
                if column == "superkingdom":
                    with self.subTest(store=label, taxon=name):
                        self.assertIn(value.strip().casefold(), self.mod.TOP_LEVEL)


if __name__ == "__main__":
    unittest.main()
