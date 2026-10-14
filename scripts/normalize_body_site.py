#!/usr/bin/env python3
"""Give body_site one vocabulary in both stores, and record where it came from.

The column held three vocabularies for one concept and the two stores disagreed
about the same rows:

    MySQL   GMrepo    NULL                    33,630
            MullerMM  NULL                     2,464
            PubMed    'feces'                    479
            MGnify    'human gastric mucosa'     120
    SQLite  GMrepo    'gut'                   23,494

scripts/cleanup_gutdb.py has filled the mirror's NULLs with 'gut' since it was
written; nothing ever did the same to MySQL. So `SELECT ... WHERE body_site =
'feces'` returns 479 rows against MySQL and none against the mirror, and
`body_site IS NULL` returns 36,094 against one and 0 against the other. That is
the same defect shape as the (SeqCode) annotations: a repair applied to one store
only.

'gut' is also the wrong word for the job. The stomach is part of the gut, so
'gut' does not separate the 120 gastric biopsies from the 36,000 stool samples --
and that separation is the only thing this column is for here. The vocabulary is
therefore the sampled material:

    feces            lumen content, voided stool
    gastric mucosa   endoscopic biopsy of the stomach wall

WHY THIS MATTERS ENOUGH TO FIX. MGYS00001188 contributes 120 gastric mucosa
samples and NO control arm of its own, so Gastritis, Metaplasia and Stomach
Neoplasms are each modelled by comparing stomach biopsies against pooled stool
controls. All three score an AUC of exactly 1.000 on CLR genus features -- and
body site alone, a single indicator, also scores exactly 1.000. The top genera do
not overlap at all: Haematomicrobium, Helicobacter, Tsukamurella, Gordonia and
Thermus in the mucosa against Bacteroides, Phocaeicola, Faecalibacterium,
Prevotella and Bifidobacterium in stool. These are not disease models.

NOTHING HERE IS PER-SAMPLE METADATA, and the script checks that rather than
assuming it: body_site is constant within every study in the database, because
each source sets it for the whole cohort. The values are filled from the source's
documented scope, which is evidence of a different kind from a measurement and is
labelled as such in the audit:

  GMrepo     curates human GUT metagenomes; its sample export carries no
             body-site column at all, so the constant is the repository's scope.
  MullerMM   every one of the 13 cohorts' own publication titles, in
             data/muller2022/*/metadata.tsv, says faecal, fecal, stool or gut.
  PubMed     the single study is the ME/CFS stool cohort; 'feces' was written by
             scripts/load_biomapai_study.py as a constant, not read from source.
  MGnify     'human gastric mucosa' arrived through the MGnify sync, whose title
             is "Gastric Mucosa-Associated Microbial Signatures of Early Gastric
             Cancer".

Writes data/body_site_audit.csv. Idempotent. Dry run by default; --apply writes.
"""
from __future__ import annotations

import csv
import os
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect

DRY_RUN = "--apply" not in sys.argv
SQLITE = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")
AUDIT = Path("data/body_site_audit.csv")
AUDIT_FIELDS = ("store", "source_database", "stored_body_site", "normalized",
                "basis", "n_samples")

# Existing values that mean the same material under a different word. 'gut' is
# folded rather than kept because it cannot separate stool from a gastric biopsy.
SYNONYMS = {
    "gut": "feces",
    "stool": "feces",
    "faeces": "feces",
    "feces": "feces",
    "human gastric mucosa": "gastric mucosa",
    "gastric mucosa": "gastric mucosa",
}

# What to fill a NULL with, per source database, and the basis for whatever the
# column ends up holding for that source -- whether this script filled it or
# found it already there.
#
# The basis is NOT "was the cell already populated". Both values that were
# already in place came from constants in this project's own code, not from the
# sources: the mirror's 'gut' was written by scripts/cleanup_gutdb.py and
# PubMed's 'feces' is a literal in scripts/load_biomapai_study.py. Deriving the
# basis from whether a cell was empty would have labelled both of those
# `reported`, which is the mislabel this column exists to avoid.
BY_SOURCE = {
    "GMrepo": ("feces", "source_scope"),
    "MullerMM": ("feces", "source_scope"),
    "PubMed": ("feces", "source_scope"),
    "MicrobiomeHD": ("feces", "source_scope"),
    # The only source that carried a body site through its own sync, per study.
    "MGnify": (None, "reported"),
}


def normalize(value):
    if value is None or not str(value).strip():
        return None
    return SYNONYMS.get(str(value).strip().casefold())


def audit_store(label, cur, placeholder, rows_out) -> dict[str, int]:
    """Fold synonyms, then fill NULLs from the source's scope. Returns counters."""
    counts = {"folded": 0, "filled": 0, "unmapped": 0}

    cur.execute("""SELECT st.source_database, s.body_site, COUNT(*)
                   FROM samples s JOIN studies st ON st.id = s.study_id
                   GROUP BY st.source_database, s.body_site""")
    for source, stored, n in cur.fetchall():
        if stored is not None and str(stored).strip():
            tidy = normalize(stored)
            if tidy is None:
                print(f"   {label}: {source} {stored!r} is not in the vocabulary "
                      f"({n} samples) -- left alone, add it to SYNONYMS")
                counts["unmapped"] += n
                rows_out.append({"store": label, "source_database": source,
                                 "stored_body_site": stored, "normalized": "",
                                 "basis": "unmapped", "n_samples": n})
                continue
            rows_out.append({"store": label, "source_database": source,
                             "stored_body_site": stored, "normalized": tidy,
                             "basis": BY_SOURCE.get(source,
                                                    (None, "unknown_source"))[1],
                             "n_samples": n})
            if tidy != stored:
                print(f"   {label}: {source} {stored!r} -> {tidy!r}  {n} samples")
                counts["folded"] += n
                if not DRY_RUN:
                    cur.execute(
                        f"""UPDATE samples s JOIN studies st ON st.id = s.study_id
                            SET s.body_site = {placeholder}
                            WHERE st.source_database = {placeholder}
                              AND s.body_site = {placeholder}"""
                        if placeholder == "%s" else
                        f"""UPDATE samples SET body_site = {placeholder}
                            WHERE body_site = {placeholder} AND study_id IN (
                              SELECT id FROM studies WHERE source_database = {placeholder})""",
                        (tidy, source, stored) if placeholder == "%s"
                        else (tidy, stored, source))
            continue

        fill, basis = BY_SOURCE.get(source, (None, "unknown_source"))
        if fill is None:
            print(f"   {label}: {source} has {n} samples with no body_site and no "
                  f"documented scope -- left NULL")
            rows_out.append({"store": label, "source_database": source,
                             "stored_body_site": "", "normalized": "",
                             "basis": basis, "n_samples": n})
            continue
        print(f"   {label}: {source} NULL -> {fill!r} ({basis})  {n} samples")
        counts["filled"] += n
        rows_out.append({"store": label, "source_database": source,
                         "stored_body_site": "", "normalized": fill,
                         "basis": basis, "n_samples": n})
        if not DRY_RUN:
            if placeholder == "%s":
                cur.execute("""UPDATE samples s JOIN studies st ON st.id = s.study_id
                               SET s.body_site = %s
                               WHERE st.source_database = %s
                                 AND (s.body_site IS NULL OR TRIM(s.body_site) = '')""",
                            (fill, source))
            else:
                cur.execute("""UPDATE samples SET body_site = ?
                               WHERE (body_site IS NULL OR TRIM(body_site) = '')
                                 AND study_id IN (SELECT id FROM studies
                                                  WHERE source_database = ?)""",
                            (fill, source))
    return counts


def check_study_constant(label, cur) -> int:
    """body_site must stay a per-study constant; a study holding two is a defect.

    Not an assumption this script relies on blindly: if a source ever starts
    reporting it per sample, filling from the source's scope would overwrite real
    metadata, and this is what would catch that.
    """
    cur.execute("""SELECT COUNT(*) FROM (
                     SELECT study_id FROM samples
                     GROUP BY study_id
                     HAVING COUNT(DISTINCT COALESCE(body_site, '~')) > 1) x""")
    n = cur.fetchone()[0]
    print(f"   {label}: studies holding more than one body_site: {n}")
    return n


def main() -> int:
    rows: list[dict] = []
    mixed = 0

    print("mysql:")
    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    mixed += check_study_constant("mysql", cur)
    counts = audit_store("mysql", cur, "%s", rows)
    if not DRY_RUN:
        conn.commit()
        print(f"   mysql: folded {counts['folded']}, filled {counts['filled']}")
    conn.close()

    if Path(SQLITE).exists():
        print("\nsqlite mirror:")
        sq = sqlite3.connect(SQLITE)
        sc = sq.cursor()
        mixed += check_study_constant("sqlite", sc)
        sq_counts = audit_store("sqlite", sc, "?", rows)
        if not DRY_RUN:
            sq.commit()
            print(f"   sqlite: folded {sq_counts['folded']}, "
                  f"filled {sq_counts['filled']}")
        sq.close()
    else:
        print(f"\nsqlite mirror {SQLITE} not present; skipped")

    AUDIT.parent.mkdir(parents=True, exist_ok=True)
    with AUDIT.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=AUDIT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {AUDIT} ({len(rows)} rows)")

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    return 1 if mixed else 0


if __name__ == "__main__":
    raise SystemExit(main())
