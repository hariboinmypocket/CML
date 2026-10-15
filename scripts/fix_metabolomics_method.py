#!/usr/bin/env python3
"""Make studies.metabolomics_method answer the question its name asks.

All 13 metabolomics studies held the literal string 'see mtb.map.tsv'. That file
exists, but it maps compound names to KEGG and HMDB identifiers -- its columns are
Compound, KEGG, HMDB, High.Confidence.Annotation -- and it contains no method, no
instrument and no units. So the column named for the method held a filename, and
the file it named held something else.

That matters because this column carries a standing instruction. Metabolite
levels are the one measurement in this database that cannot be pooled across
studies: instruments, targeted versus untargeted designs, extraction protocols
and units all differ between cohorts, and nothing sums to 1, so the rank
invariant has no analogue. The documented rule is to check metabolomics_method
before combining anything -- and that rule was uncheckable, because every study
gave the same answer and the answer was a pointer to the wrong file.

The method itself is NOT in the Muller collection; it redistributes the matrices,
not the protocols. What the collection does carry, in every cohort's
metadata.tsv, is the DOI of the paper that describes the protocol. That is the
right thing to store here: it identifies the protocol uniquely, it makes the
comparability rule mechanical -- levels are comparable only where the DOI matches
-- and it does not claim to be a method when no method was published to us.

Values become 'doi:<doi>' so the column is self-describing and nobody mistakes it
for an instrument name.

Reads data/muller2022/*/metadata.tsv. Dry run by default; --apply writes.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect

DRY_RUN = "--apply" not in sys.argv
COLLECTION = Path("data/muller2022")
SOURCE_DB = "MullerMM"


def doi_by_cohort() -> dict[str, str]:
    """One DOI per cohort directory, taken from its own metadata.tsv.

    Read from every row rather than just the first: a cohort whose rows disagreed
    about the DOI would mean two protocols pooled under one study accession, and
    that is worth refusing rather than averaging.
    """
    out: dict[str, str] = {}
    for directory in sorted(COLLECTION.glob("*/")):
        meta = directory / "metadata.tsv"
        if not meta.exists():
            continue
        with meta.open(newline="", encoding="utf-8-sig") as fh:
            dois = {(row.get("DOI") or "").strip()
                    for row in csv.DictReader(fh, delimiter="\t")}
        dois.discard("")
        if len(dois) > 1:
            print(f"   {directory.name}: {len(dois)} different DOIs in one cohort "
                  f"{sorted(dois)} -- skipped, this needs a human")
            continue
        if dois:
            out[directory.name] = dois.pop()
    return out


def main() -> int:
    dois = doi_by_cohort()
    print(f"cohorts with a DOI in their metadata: {len(dois)}")

    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    cur.execute("""SELECT project_accession, metabolomics_method
                   FROM studies
                   WHERE source_database = %s AND metabolomics_method IS NOT NULL""",
                (SOURCE_DB,))
    stored = cur.fetchall()
    print(f"studies carrying a metabolomics_method: {len(stored)}")

    changed = missing = 0
    for accession, current in stored:
        cohort = accession.removeprefix(f"{SOURCE_DB}_")
        doi = dois.get(cohort)
        if not doi:
            print(f"   {accession}: no DOI in the collection -- left as {current!r}")
            missing += 1
            continue
        value = f"doi:{doi}"
        if value == current:
            continue
        print(f"   {accession}: {current!r} -> {value!r}")
        changed += 1
        if not DRY_RUN:
            cur.execute("UPDATE studies SET metabolomics_method = %s "
                        "WHERE project_accession = %s", (value, accession))

    if not DRY_RUN:
        conn.commit()
        print(f"\nupdated {changed} studies ({missing} left unchanged)")

    # The point of the column is that a DOI tells two protocols apart, so a
    # duplicate DOI across accessions is the one thing worth reporting back.
    cur.execute("""SELECT metabolomics_method, COUNT(*) n,
                          GROUP_CONCAT(project_accession) accessions
                   FROM studies
                   WHERE metabolomics_method LIKE 'doi:%%'
                   GROUP BY metabolomics_method HAVING n > 1""")
    shared = cur.fetchall()
    if shared:
        print("\ntwo studies sharing one protocol DOI (their levels ARE comparable):")
        for method, n, accessions in shared:
            print(f"   {method}  x{n}  {accessions}")
    else:
        print("\nevery metabolomics study has its own protocol DOI; no two "
              "cohorts' levels are comparable")
    conn.close()

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
