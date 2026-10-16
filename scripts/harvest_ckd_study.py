#!/usr/bin/env python3
"""Verification-gated CKD demographic synchronization (dry run by default).

The former numeric-prefix joins for HC and CKD were unverified and are retired.
Only the pinned verification manifest can authorize a clinical-to-run link.
Delegate to the shared two-store repair so this entry point cannot reintroduce
withdrawn values or silently update only MySQL. No network or prefix guessing.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

from repair_ckd_demographics import AUDIT, REPORT, digest, main


if __name__ == "__main__":
    # After pruning, use the original audit behind the cleanup report.
    # The delegated repair checks its hash against the verification manifest.
    if "--audit-source" not in sys.argv and REPORT.exists():
        report = json.loads(REPORT.read_text())
        if report.get("commit_state") != "complete":
            raise SystemExit("Incomplete cleanup; inspect the recovery report before importing")
        if digest(AUDIT) != report["active_audit_sha256"]:
            raise SystemExit("Active CKD audit changed since verification; reverify before importing")
        source = Path(report.get("source_audit_path") or
                      str(Path(report["backup_directory"]) / AUDIT.name))
        sys.argv.extend(["--audit-source", str(source)])
    raise SystemExit(main())
