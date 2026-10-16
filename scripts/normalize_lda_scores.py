#!/usr/bin/env python3
"""Add a generated, three-decimal softsign LDA score. Dry run unless --apply.

MySQL is authoritative. Use --sqlite to update the local mirror separately.
Raw association fields are fingerprinted before/after; every derived value is
checked against an independent Decimal calculation. No source score is edited.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import sqlite3
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.lda import install_normalized_lda, normalize_lda_score


def raw_fingerprint(connection) -> str:
    cursor = connection.cursor()
    try:
        cursor.execute("""SELECT id, taxon_id, disease_id, comparison_id,
            direction, effect_type, effect_size, p_value, q_value, source_database
            FROM taxon_disease_associations ORDER BY id""")
        digest = hashlib.sha256()
        for row in cursor:
            digest.update(repr(tuple(row)).encode("utf-8"))
            digest.update(b"\n")
        return digest.hexdigest()
    finally:
        cursor.close()


def verify(connection) -> None:
    cursor = connection.cursor()
    try:
        cursor.execute("""SELECT id, effect_type, effect_size, lda_score_normalized
            FROM taxon_disease_associations""")
        count = 0
        for row_id, effect_type, score, actual in cursor:
            expected = normalize_lda_score(score) if effect_type == "LDA" else None
            if ((expected is None) != (actual is None)
                    or (expected is not None and abs(float(expected) - float(actual)) > 1e-9)):
                raise ValueError(f"Normalized LDA mismatch at association {row_id}")
            count += expected is not None
        print(f"Verified {count:,} normalized LDA scores; other effect types remain NULL.")
        cursor.execute("""SELECT MIN(lda_score_normalized), MAX(lda_score_normalized)
            FROM taxon_disease_associations""")
        lo, hi = cursor.fetchone()
        print(f"Normalized range: {lo} to {hi}")
    finally:
        cursor.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--env-file", default=".env")
    parser.add_argument("--sqlite", type=Path, metavar="PATH")
    args = parser.parse_args()
    if args.sqlite:
        # mode=rw prevents a typo from silently creating an empty database.
        connection = sqlite3.connect(args.sqlite.resolve().as_uri() + "?mode=rw", uri=True)
    else:
        connection = connect(Settings.from_env(args.env_file))
    try:
        cursor = connection.cursor()
        cursor.execute("""SELECT COUNT(*), MIN(effect_size), MAX(effect_size)
            FROM taxon_disease_associations WHERE effect_type = 'LDA'""")
        print("LDA count and current range:", cursor.fetchone())
        cursor.close()
        if not args.apply:
            print("Dry run. Formula: ROUND(score / (1 + ABS(score)), 3). Use --apply to install.")
            return
        before = raw_fingerprint(connection)
        if args.sqlite:
            connection.execute("BEGIN IMMEDIATE")
        added = install_normalized_lda(connection, sqlite=bool(args.sqlite))
        verify(connection)
        if raw_fingerprint(connection) != before:
            raise ValueError("Raw association data changed during migration; inspect before continuing")
        connection.commit()
        print(f"Column {'added' if added else 'already present'}; summary view updated.")
        print("All original association fields unchanged (SHA-256 verified).")
    except Exception:
        connection.rollback()
        raise
    finally:
        connection.close()


if __name__ == "__main__":
    main()
