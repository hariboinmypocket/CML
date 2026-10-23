#!/usr/bin/env python3
"""Seed demographic_withdrawals from the CKD link verification.

A cleared cell and a never-measured cell look identical in the data, so the
withdrawal has to be recorded somewhere the pipeline reads. Derived from
results/ckd_link_verification/verification.json rather than from a hand-written
list: a run in PRJNA562327 keeps its demographics only where the
clinical-to-run link is verified against public ENA aliases, which is the same
predicate the original repair applied. Regenerating this is therefore a
re-derivation, not a transcription of whatever the table happens to hold.

Preview by default; --apply writes to MySQL and refreshes the CSV the mirror
build reads.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.transform import verified_demographic_link

PROJECT_ACCESSION = "PRJNA562327"
FIELDS = ("sex", "age_years", "bmi")
REASON = "ena_link_unverified_PRJNA562327"
VERIFICATION = Path("results/ckd_link_verification/verification.json")
LEDGER = Path("data/ckd_demographic_withdrawals.csv")


def derive(connection, verification: Path) -> tuple[list[tuple[str, str, str]], int]:
    records = json.loads(verification.read_text())["records"]
    verified = {r["run_accession"] for r in records if verified_demographic_link(r)}
    cursor = connection.cursor(dictionary=True)
    try:
        cursor.execute(
            "SELECT s.run_accession, s.sex, s.age_years, s.bmi FROM samples s "
            "JOIN studies st ON st.id = s.study_id "
            "WHERE st.project_accession = %s", (PROJECT_ACCESSION,))
        samples = cursor.fetchall()
    finally:
        cursor.close()
    # A stored value outside the verified set means the withdrawal has already
    # been undone; seeding over it would hide that rather than report it.
    stray = [(r["run_accession"], f) for r in samples for f in FIELDS
             if r[f] is not None and r["run_accession"] not in verified]
    return ([(r["run_accession"], f, REASON) for r in samples for f in FIELDS
             if r["run_accession"] not in verified], len(stray))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env-file", default=".env")
    parser.add_argument("--verification", type=Path, default=VERIFICATION)
    parser.add_argument("--ledger", type=Path, default=LEDGER)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    connection = connect(Settings.from_env(args.env_file))
    try:
        rows, stray = derive(connection, args.verification)
        runs = len({r[0] for r in rows})
        print(f"{PROJECT_ACCESSION}: withholding {len(rows):,} cells "
              f"across {runs:,} unverified runs")
        if stray:
            print(f"WARNING: {stray} stored value(s) sit outside the verified set; "
                  "the withdrawal has already been undone somewhere")
        if not args.apply:
            print("Preview only. Use --apply to seed the ledger.")
            return
        cursor = connection.cursor()
        try:
            cursor.executemany(
                "INSERT INTO demographic_withdrawals (run_accession, field, reason) "
                "VALUES (%s, %s, %s) ON DUPLICATE KEY UPDATE reason = VALUES(reason)",
                rows)
            connection.commit()
            cursor.execute("SELECT COUNT(*) FROM demographic_withdrawals")
            total = cursor.fetchone()[0]
        finally:
            cursor.close()
    finally:
        connection.close()

    args.ledger.parent.mkdir(parents=True, exist_ok=True)
    with args.ledger.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["run_accession", "field", "reason"])
        writer.writerows(sorted(rows))
    print(f"MySQL ledger rows: {total:,}")
    print(f"Ledger CSV for the mirror build: {args.ledger.resolve()}")


if __name__ == "__main__":
    main()
