#!/usr/bin/env python3
"""Refresh the GMrepo sample CSV from MySQL metadata and remove raw/helper fields.

The database is read only. Preview by default; --apply preserves an original
backup and a field-level audit before atomically replacing the CSV.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import os
from pathlib import Path
import shutil
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.pipeline import (GMREPO_SAMPLE_FIELDS, prepare_gmrepo_sample_export,
                            write_sample_export_audit)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--path", type=Path, default=Path("data/gmrepo_project_samples.csv"))
    parser.add_argument("--env-file", default=".env")
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    with args.path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    runs = [r["run_id"] for r in rows]
    if len(set(runs)) != len(runs):
        raise ValueError("Duplicate run IDs; inspect before cleanup")
    connection = connect(Settings.from_env(args.env_file))
    try:
        connection.start_transaction(readonly=True, consistent_snapshot=True)
        output, audit = prepare_gmrepo_sample_export(connection, rows)
    finally:
        connection.rollback()
        connection.close()
    assert [r["run_id"] for r in output] == runs
    assert [r["qc_status"] for r in output] == [r["qc_status"] for r in rows]
    for field in ("sex", "age_years", "bmi"):
        before = sum(bool(r.get(field)) for r in rows)
        after = sum(bool(r.get(field)) for r in output)
        filled = sum(a["field"] == field and not a["before"] and bool(a["after"]) for a in audit)
        cleared = sum(a["field"] == field and bool(a["before"]) and not a["after"] for a in audit)
        print(f"{field}: present {before:,} -> {after:,}; filled {filled:,}; "
              f"invalid values cleared {cleared:,}; remaining missing {len(output)-after:,}")
    print("Label changes:", dict(Counter(a["field"] for a in audit
                                       if a["field"] in {"disease_name", "mesh_id"})))
    if not args.apply:
        print("Preview only. Use --apply to save the cleaned export.")
        return
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    backup_dir = args.path.parent / "archive" / ("gmrepo_samples_cleanup_" + stamp)
    backup_dir.mkdir(parents=True, exist_ok=False)
    backup = backup_dir / args.path.name
    shutil.copy2(args.path, backup)
    audit_path = backup_dir / "field_changes.csv"
    write_sample_export_audit(audit_path, audit)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                         dir=args.path.parent, delete=False) as handle:
            temporary = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=GMREPO_SAMPLE_FIELDS)
            writer.writeheader()
            writer.writerows(output)
        with temporary.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            assert reader.fieldnames == GMREPO_SAMPLE_FIELDS
            assert list(reader) == [{k: str(v) if v is not None else "" for k, v in r.items()}
                                    for r in output]
        os.replace(temporary, args.path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    print(f"Saved {len(output):,} rows, {len(GMREPO_SAMPLE_FIELDS)} columns: {args.path.resolve()}")
    print(f"Original backup: {backup.resolve()}")
    print(f"Field-level audit: {audit_path.resolve()}")


if __name__ == "__main__":
    main()
