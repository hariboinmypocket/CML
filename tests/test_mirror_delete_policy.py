"""A repair script must not turn on SQLite's foreign keys before deleting taxa.

This looks backwards, so the reasoning is recorded with the numbers it was
measured on. `sample_taxon_abundances.taxon_id` and
`taxon_disease_associations.taxon_id` both declare ON DELETE CASCADE in
gutdb/schema_sqlite.sql, and SQLite enforces foreign keys only when a
connection asks -- the default is off.

Deleting one taxon that held 16,815 abundance rows and 111 associations:

    foreign_keys = OFF   the delete succeeds and leaves 16,815 orphan
                         abundances and 111 orphan associations, so
                         orphan_abundances and orphan_associations both fire and
                         the pre-push hook refuses the push.

    foreign_keys = ON    the delete succeeds, cascades, and removes all 16,926
                         rows. Both counters read zero. Validate is clean and
                         the data is gone.

So enabling enforcement in these scripts would convert a loud, gated failure
into a silent one. It is the same mechanism the README records for MySQL, where
foreign keys are always enforced and `fk_tma_taxon ... ON DELETE CASCADE` is how
the adjudication lost a verdict with nothing noticing.

The safety model of a repair script is: prove the row carries no evidence, then
delete it. If that proof is ever wrong, orphans are what make the mistake
visible. Cascade would destroy the evidence of the mistake along with the data.

This is not a blanket rule. build_gutdb_sqlite.py prunes `samples` and
cleanup_gutdb.py de-duplicates `diseases`, both with enforcement on and both
intending the cascade -- cleanup_gutdb even deletes the child associations
itself first. The rule is specific to deleting `taxa` from the mirror in a
repair.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCHEMA = ROOT / "gutdb" / "schema_sqlite.sql"
CASCADING_CHILDREN = ("sample_taxon_abundances", "taxon_disease_associations")


def scripts_deleting_taxa():
    out = []
    for path in sorted((ROOT / "scripts").glob("*.py")):
        text = path.read_text()
        if "sqlite3.connect" not in text:
            continue
        if re.search(r"DELETE\s+FROM\s+`?taxa`?", text, re.I):
            out.append(path)
    return out


class MirrorDeletePolicyTests(unittest.TestCase):
    def test_the_premise_still_holds_cascade_is_declared(self):
        """If the schema ever stops cascading, this policy needs rethinking."""
        text = SCHEMA.read_text()
        for child in CASCADING_CHILDREN:
            with self.subTest(table=child):
                block = re.search(
                    rf"CREATE TABLE(?: IF NOT EXISTS)? {child}\b(.*?);",
                    text, re.S | re.I)
                self.assertIsNotNone(block, f"{child} not found in the schema")
                self.assertRegex(
                    block.group(1),
                    r"taxon_id[^,]*REFERENCES\s+taxa\(id\)\s+ON DELETE CASCADE",
                    f"{child}.taxon_id no longer cascades; re-read this test's "
                    f"docstring before changing the policy")

    def test_the_orphan_counters_that_make_the_policy_work_exist(self):
        from gutdb.pipeline import _validate_checks
        checks = _validate_checks()
        for name in ("orphan_abundances", "orphan_associations"):
            with self.subTest(counter=name):
                self.assertIn(name, checks)

    def test_no_repair_script_enables_foreign_keys_before_deleting_taxa(self):
        offenders = []
        for path in scripts_deleting_taxa():
            if re.search(r"PRAGMA\s+foreign_keys\s*=\s*ON", path.read_text(), re.I):
                offenders.append(path.name)
        self.assertEqual(
            [], offenders,
            "these delete taxa from the mirror with foreign keys enforced, so a "
            "wrong delete cascades silently instead of leaving orphans the "
            "counters catch: %s" % ", ".join(offenders))

    def test_there_are_scripts_this_applies_to(self):
        """Without this the rule above could pass by covering nothing."""
        self.assertGreaterEqual(len(scripts_deleting_taxa()), 4)


if __name__ == "__main__":
    unittest.main()
