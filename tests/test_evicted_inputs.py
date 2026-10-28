"""An input file that is not on this disk must fail, not hang.

The Desktop is an iCloud-synced folder. When the disk filled to 98%, macOS
reclaimed space by evicting this project's large inputs to iCloud: 71 files and
1.6GB, among them the SQLite mirror, both NCBI taxdump files, every GMrepo
abundance CSV and every muller2022 table. An evicted file keeps all of its
metadata, so exists(), st_size and mtime stay correct and open() still
succeeds. Only the first read blocks, and with no room to land the download it
blocks forever.

That is the least visible failure this project has had. It has no exception, no
log line and no counter -- the process simply stops, which read as failing
hardware for a whole session. A rebuild of the mirror from the data/ CSVs would
have stalled the same way, because those CSVs were evicted too.

The checks here are on the flag, never on st_blocks: APFS compression also
reports zero allocated blocks, and most of this repo's own source is compressed
and perfectly local, so st_blocks would condemn healthy files.
"""
import importlib.util
import os
import sys
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gutdb import config, taxdump
from gutdb.config import SF_DATALESS, EvictedFileError, require_local
from gutdb.pipeline import open_dict_rows

UF_COMPRESSED = 0x00000020

_spec = importlib.util.spec_from_file_location(
    "build_gutdb_sqlite", ROOT / "scripts" / "build_gutdb_sqlite.py"
)
build_mirror = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(build_mirror)

_cleanup_spec = importlib.util.spec_from_file_location(
    "cleanup_gutdb", ROOT / "scripts" / "cleanup_gutdb.py"
)
cleanup_gutdb = importlib.util.module_from_spec(_cleanup_spec)
_cleanup_spec.loader.exec_module(cleanup_gutdb)


def first_evicted_input():
    """A genuinely dataless file under data/, if this machine has one."""
    for path in (ROOT / "data").rglob("*"):
        if not path.is_file():
            continue
        flags = config._file_flags(path)
        if flags is not None and flags & SF_DATALESS:
            return path
    return None


class RequireLocalTests(unittest.TestCase):
    def test_a_dataless_file_is_refused(self):
        target = ROOT / "README.md"
        with mock.patch.object(config, "_file_flags", return_value=SF_DATALESS):
            with self.assertRaises(EvictedFileError):
                require_local(target)

    def test_a_local_file_passes_through(self):
        target = ROOT / "README.md"
        self.assertEqual(require_local(target), target)

    def test_compression_alone_is_not_eviction(self):
        """The confound that hid this for a session.

        gutdb/pipeline.py is stored compressed and reports zero allocated
        blocks while being entirely local, so a blocks-based check would
        refuse to read the pipeline's own source.
        """
        with mock.patch.object(config, "_file_flags", return_value=UF_COMPRESSED):
            target = ROOT / "gutdb" / "pipeline.py"
            self.assertEqual(require_local(target), target)

    def test_a_platform_without_the_field_is_not_second_guessed(self):
        with mock.patch.object(config, "_file_flags", return_value=None):
            target = ROOT / "README.md"
            self.assertEqual(require_local(target), target)

    def test_the_message_says_which_file_and_how_to_recover_it(self):
        target = ROOT / "data" / "gmrepo_species_abundances.csv"
        with mock.patch.object(config, "_file_flags", return_value=SF_DATALESS):
            with self.assertRaises(EvictedFileError) as caught:
                require_local(target)
        message = str(caught.exception)
        self.assertIn("gmrepo_species_abundances.csv", message)
        self.assertIn("brctl download", message)


class EveryInputReaderIsGuardedTests(unittest.TestCase):
    """Each place that opens a data file for reading must check first.

    Asserted by calling the readers, not by reading their source, so moving a
    read to a new helper cannot quietly drop the guard.
    """

    def setUp(self):
        patcher = mock.patch.object(config, "_file_flags", return_value=SF_DATALESS)
        patcher.start()
        self.addCleanup(patcher.stop)
        self.target = ROOT / "README.md"

    def test_open_dict_rows_is_guarded(self):
        # A generator, so it must be iterated before its body runs.
        with self.assertRaises(EvictedFileError):
            next(iter(open_dict_rows(self.target)))

    def test_the_taxdump_node_reader_is_guarded(self):
        with self.assertRaises(EvictedFileError):
            taxdump.SpeciesNodes(self.target)

    def test_the_mirror_builders_reader_is_guarded(self):
        with self.assertRaises(EvictedFileError):
            next(iter(build_mirror.read_rows(self.target)))

    def test_the_mirror_only_cleanups_trait_sources_are_guarded(self):
        """cleanup_gutdb.py reads three files the build does not.

        Two of them -- mimedb_microbes_v2.csv and
        gut_microbiome_925_oxygen_pH_spore.csv -- are evicted on this machine
        right now, so running it would have hung rather than failed. It is the
        script that back-fills the taxa traits, and it reaches only the mirror,
        so it is the one most likely to be run to close that gap.

        build_gutdb_sqlite.py's schema read is guarded the same way but has no
        test here on purpose: it lives in main(), which deletes the mirror
        before reading, so calling it from a test would destroy the database.
        """
        # Each of the three reads is checked on its own. Mocking every file
        # dataless would not do it: the first guard to fire ends the function,
        # so two of the three could lose their guard and the test would still
        # pass -- which it did, until this was written this way. Nor can
        # _file_flags be mocked per path, because mimedb_microbes_v2.csv really
        # is evicted here, and reporting it local would hang the read it is
        # meant to prevent. So the guard is spied on and the reads stubbed.
        seen = []

        def spy(path):
            seen.append(str(path))
            return path

        empty = pd.DataFrame()
        with mock.patch.object(cleanup_gutdb, "require_local", side_effect=spy), \
                mock.patch.object(cleanup_gutdb.pd, "read_csv", return_value=empty):
            try:
                cleanup_gutdb.load_local_trait_sources()
            except Exception:
                # An empty frame will not satisfy the body; only the guard
                # calls made before that point matter here.
                pass
        for source in ("final_taxa.tsv", "mimedb_microbes_v2.csv",
                       "gut_microbiome_925_oxygen_pH_spore.csv"):
            with self.subTest(source=source):
                self.assertTrue(
                    any(source in path for path in seen),
                    "%s is read without going through require_local" % source)


class RealEvictedFileTests(unittest.TestCase):
    def test_a_genuinely_evicted_input_is_refused(self):
        """Against a real dataless file, with nothing mocked.

        Skipped once iCloud has materialised everything again, or once
        "Optimize Mac Storage" is off -- which is the state to aim for.
        """
        evicted = first_evicted_input()
        if evicted is None:
            self.skipTest("no evicted file under data/ on this machine")
        self.assertGreater(evicted.stat().st_size, 0)
        self.assertTrue(evicted.exists(), "an evicted file still looks present")
        with self.assertRaises(EvictedFileError):
            next(iter(open_dict_rows(evicted)))


if __name__ == "__main__":
    unittest.main()
