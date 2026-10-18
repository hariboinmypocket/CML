#!/usr/bin/env python3
"""Withdraw unverified Ren 2020 clinical links; synchronize verified cells.

Dry run by default. --apply backs up both active CSVs, provenance, a consistent
SQLite snapshot and targeted MySQL records before any demographic update.
For a repeat run, --audit-source accepts the original audit in that backup.
No sample, abundance, disease or taxon row is deleted or repointed.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import shutil
import sqlite3
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from gutdb.config import Settings
from gutdb.db import connect
from gutdb.transform import clean_str, reconcile_demographic_cells, same_demographic_value, verified_demographic_link

PROJECT = "PRJNA562327"
AUDIT = ROOT / "data/ckd_study_demographics_audit.csv"
SAMPLES = ROOT / "data/gmrepo_project_samples.csv"
MIRROR = ROOT / "gut_microbiome.sqlite"
VERIFICATION = ROOT / "results/ckd_link_verification/verification.json"
REPORT = ROOT / "results/ckd_link_verification/cleanup.json"
FIELDS = ("sex", "age_years", "bmi")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


def csv_bytes(fields, rows):
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue().encode("utf-8")


def census(connection, mysql):
    cur = connection.cursor()
    result = {}
    for table in ("samples", "sample_taxon_abundances", "taxa", "taxon_disease_associations"):
        cur.execute(f"SELECT COUNT(*) FROM {table}")
        result[table] = cur.fetchone()[0]
    cur.close()
    return result


def fetch_targets(connection, mysql, runs, lock=False):
    placeholder = "%s" if mysql else "?"
    cur = connection.cursor(dictionary=True) if mysql else connection.cursor()
    query = """SELECT s.*, t.project_accession AS repair_project
               FROM samples s JOIN studies t ON t.id=s.study_id
               WHERE s.run_accession IN (""" + ",".join([placeholder] * len(runs)) + ")"
    if mysql and lock:
        query += " FOR UPDATE"
    cur.execute(query, tuple(runs))
    rows = [dict(row) for row in cur.fetchall()]
    cur.close()
    result = {row["run_accession"]: row for row in rows}
    if len(rows) != len(runs) or set(result) != set(runs):
        raise RuntimeError("Missing or duplicated target run in a database store")
    if any(row["repair_project"] != PROJECT for row in rows):
        raise RuntimeError("Unexpected project identity; refusing repair")
    return result


def plan(current, sources, evidence):
    edits, conflicts, totals = {}, {}, Counter()
    for run, source in sources.items():
        verified = verified_demographic_link(evidence[run])
        changes, differing = reconcile_demographic_cells(current[run], source, verified)
        if changes:
            edits[run] = changes
        if differing:
            conflicts[run] = differing
        for field in changes:
            totals[("filled_" if verified else "cleared_") + field] += 1
    return edits, conflicts, dict(totals)


def apply_sql(connection, mysql, edits):
    placeholder = "%s" if mysql else "?"
    cur = connection.cursor()
    for run, changes in edits.items():
        fields = list(changes)
        cur.execute("UPDATE samples SET " + ", ".join(f"{field}={placeholder}" for field in fields)
                    + f" WHERE run_accession={placeholder}",
                    tuple(changes[field] for field in fields) + (run,))
        if cur.rowcount != 1:
            raise RuntimeError("Target changed during repair")
    cur.close()


def assert_after(before, after, edits):
    for run, original in before.items():
        expected = dict(original)
        expected.update(edits.get(run, {}))
        actual = after[run]
        # SQLite REAL and MySQL DECIMAL differ, but each is compared in its own store.
        other_keys = set(expected) - set(FIELDS)
        same_cells = all(
            (not clean_str(expected[field]) and not clean_str(actual[field]))
            or same_demographic_value(field, expected[field], actual[field])
            for field in FIELDS
        )
        if not same_cells or any(expected[key] != actual[key] for key in other_keys):
            raise RuntimeError(f"Unexpected target-row change: {run}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--audit-source", type=Path, default=AUDIT)
    args = parser.parse_args()
    verification = json.loads(VERIFICATION.read_text())
    if digest(args.audit_source) != verification["input_audit_sha256"]:
        parser.error("Audit differs from verified snapshot; use --audit-source ORIGINAL_BACKUP_AUDIT")
    audit_fields, source_rows = read_csv(args.audit_source)
    sources = {row["run_accession"]: row for row in source_rows}
    evidence = {row["run_accession"]: row for row in verification["records"]}
    if len(sources) != len(source_rows) or set(sources) != set(evidence):
        raise RuntimeError("Verification/audit run coverage differs")
    for run, row in sources.items():
        if any(row[key] != evidence[run][key] for key in ("supplement_id", "library_name")):
            raise RuntimeError("Verification identifier differs from audit")
    kept = []
    for row in source_rows:
        if verified_demographic_link(evidence[row["run_accession"]]):
            kept.append({**row, "join_basis": "source_alias_matches_after_separator_normalization"})
    sample_fields, samples = read_csv(SAMPLES)
    sample_index = {row["run_id"]: row for row in samples if row["run_id"] in sources}
    if len(sample_index) != len(sources) or sum(row["run_id"] in sources for row in samples) != len(sources):
        raise RuntimeError("Missing or duplicated target run in sample CSV")
    if any(row["project_id"] != PROJECT for row in sample_index.values()):
        raise RuntimeError("Unexpected sample CSV project")
    pinned_hashes = {path: digest(path) for path in (AUDIT, SAMPLES)}
    csv_edits, csv_conflicts, csv_counts = plan(sample_index, sources, evidence)
    for run, changes in csv_edits.items():
        sample_index[run].update({key: "" if value is None else value for key, value in changes.items()})

    mysql = connect(Settings.from_env(ROOT / ".env"))
    mirror = sqlite3.connect(f"file:{MIRROR}?mode=rw", uri=True, timeout=60)
    mirror.row_factory = sqlite3.Row
    backup = None
    try:
        mysql.start_transaction(isolation_level="REPEATABLE READ")
        mirror.execute("BEGIN IMMEDIATE" if args.apply else "BEGIN")
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        before_counts = {"mysql": census(mysql, True), "sqlite": census(mirror, False)}
        originals = {"mysql": fetch_targets(mysql, True, sources, args.apply),
                     "sqlite": fetch_targets(mirror, False, sources)}
        plans = {name: plan(rows, sources, evidence) for name, rows in originals.items()}
        report = {"snapshot_utc": timestamp, "apply": args.apply, "project": PROJECT,
                  "audit_rows_before": len(source_rows), "audit_rows_kept": len(kept),
                  "audit_rows_removed": len(source_rows) - len(kept),
                  "census_before": before_counts,
                  "changes": {name: value[2] for name, value in plans.items()} | {"csv": csv_counts},
                  "preserved_conflicts": {name: value[1] for name, value in plans.items()} | {"csv": csv_conflicts}}
        if not args.apply:
            print(json.dumps(report, indent=2))
            return 0

        backup = ROOT / "data/archive" / ("ckd_verified_cleanup_" + timestamp)
        backup.mkdir(parents=True, exist_ok=False)
        for path in (AUDIT, SAMPLES, ROOT / "data/gmrepo_demographic_provenance.csv", VERIFICATION):
            if path.exists():
                shutil.copy2(path, backup / path.name)
        # A later repair's active audit is already pruned; retain the full
        # hash-verified input separately so every generation remains rerunnable.
        shutil.copy2(args.audit_source, backup / "verified_input_audit.csv")
        # Separate reader + SQLite backup API includes committed WAL pages.
        # The held BEGIN IMMEDIATE blocks concurrent writers during the snapshot.
        with sqlite3.connect(f"file:{MIRROR}?mode=ro", uri=True) as reader:
            with sqlite3.connect(backup / MIRROR.name) as destination:
                reader.backup(destination)
        (backup / "mysql_target_rows.json").write_text(json.dumps(originals["mysql"], default=str, indent=2))
        report["backup_directory"] = str(backup)
        report["source_audit_path"] = str(backup / "verified_input_audit.csv")
        report["commit_state"] = "backed_up_no_commits"
        (backup / "repair.json").write_text(json.dumps(report, indent=2))
        staged = {AUDIT: csv_bytes(audit_fields, kept), SAMPLES: csv_bytes(sample_fields, samples)}
        for name, connection, is_mysql in (("mysql", mysql, True), ("sqlite", mirror, False)):
            apply_sql(connection, is_mysql, plans[name][0])
            assert_after(originals[name], fetch_targets(connection, is_mysql, sources), plans[name][0])
        report["census_after"] = {"mysql": census(mysql, True), "sqlite": census(mirror, False)}
        if report["census_after"] != before_counts:
            raise RuntimeError("Unexpected census change during repair")
        for path, expected in pinned_hashes.items():
            if digest(path) != expected:
                raise RuntimeError("Active CSV changed concurrently; refusing file replacement")
        # Two stores/files cannot commit atomically. Persistent recovery report
        # records each step, and compare-and-clear/fill is safe to run again.
        mysql.commit()
        report["commit_state"] = "mysql_committed"
        (backup / "repair.json").write_text(json.dumps(report, indent=2))
        mirror.commit()
        report["commit_state"] = "both_databases_committed"
        (backup / "repair.json").write_text(json.dumps(report, indent=2))
        for path, data in staged.items():
            if digest(path) != pinned_hashes[path]:
                raise RuntimeError("CSV changed after SQL commit; recover with archived audit")
            temporary = path.with_suffix(path.suffix + ".ckd-repair-tmp")
            temporary.write_bytes(data)
            temporary.replace(path)
        report["commit_state"] = "complete"
        report["active_audit_sha256"] = digest(AUDIT)
        report["active_samples_sha256"] = digest(SAMPLES)
        (backup / "repair.json").write_text(json.dumps(report, indent=2))
        REPORT.write_text(json.dumps(report, indent=2))
        print(json.dumps(report, indent=2))
        return 0
    finally:
        mysql.rollback()
        mirror.rollback()
        mysql.close()
        mirror.close()


if __name__ == "__main__":
    raise SystemExit(main())
