#!/usr/bin/env python3
"""Fill empty lineage cells from the phylum the row's own genus already has.

143 taxa carry no phylum. scripts/backfill_lineage.py resolves lineage from the
NCBI dump, which cannot help most of these: 91 of them have no tax ID at all.
But 62 do not need NCBI, because the answer is already in this database -- rows
sharing their genus_key carry a phylum and they do not. That is an internal
inconsistency, not missing knowledge:

    Waltera                            genus     7,944 abundance rows
    Merdimmobilis hominis              species   3,195
    Faecalibacterium butyricigenerans  species   3,169   9 faecalibacterium peers say Firmicutes
    Pseudoruminococcus                 genus     2,541
    Erysipelatoclostridium             genus     1,851

THE GUARD IS THE POINT. A genus_key is derived from a name's first token, so it
can group organisms that are not relatives, and filling from an ambiguous key
would write a confident wrong lineage onto a well-measured row. Two keys are
ambiguous today:

  `bacteroides` carries both Bacteroidetes (77 rows) and Firmicutes (1), and the
  Firmicutes one is RIGHT: NCBI writes it `[Bacteroides] pectinophilus`, the
  brackets meaning the name is misapplied, and the organism is a Clostridium
  relative. One name whose genus token contradicts its lineage.

  `unclassified` carries Proteobacteria and Firmicutes across three MiMeDB rows
  named "Unclassified <family> bacterium ...", where the first token is the word
  "Unclassified". Both lineages are correct for their families; the key is simply
  meaningless for those rows.

Neither currently has a phylum-less row waiting on it, so nothing is skipped
today. The guard exists so that when one arrives it is skipped rather than
guessed, which is the same rule normalize_reclassified_phyla.py applies to
reclassifications: a value only moves when the evidence is unambiguous.

Each rank is filled independently and only where the surviving peers agree on
exactly one value. Both stores. Idempotent. Dry run by default; --apply writes.
"""
from __future__ import annotations

import csv
import os
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect

DRY_RUN = "--apply" not in sys.argv
SQLITE = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")
AUDIT = Path("data/lineage_from_genus_audit.csv")
AUDIT_FIELDS = ("store", "taxon_id", "scientific_name", "taxonomic_rank", "column",
                "filled_with", "peer_rows", "abundance_rows", "note")
RANKS = ("superkingdom", "phylum", "class_name", "order_name", "family")


def fill_store(label: str, cur, placeholder: str, rows_out: list[dict]) -> Counter:
    counts: Counter = Counter()
    for column in RANKS:
        # The consensus value per genus_key, and how many rows back it.
        cur.execute(f"""SELECT genus_key, {column}, COUNT(*)
                        FROM taxa
                        WHERE {column} IS NOT NULL AND TRIM({column}) <> ''
                          AND genus_key IS NOT NULL AND TRIM(genus_key) <> ''
                        GROUP BY genus_key, {column}""")
        by_key: dict[str, dict[str, int]] = defaultdict(dict)
        for key, value, n in cur.fetchall():
            by_key[key][value] = n

        cur.execute(f"""SELECT id, scientific_name, taxonomic_rank, genus_key
                        FROM taxa
                        WHERE ({column} IS NULL OR TRIM({column}) = '')
                          AND genus_key IS NOT NULL AND TRIM(genus_key) <> ''""")
        targets = cur.fetchall()
        updates: list[tuple[str, int]] = []
        for taxon_id, name, rank, key in targets:
            options = by_key.get(key, {})
            if not options:
                counts[f"{column}:no_peer_has_it"] += 1
                continue
            if len(options) > 1:
                counts[f"{column}:peers_disagree"] += 1
                rows_out.append({
                    "store": label, "taxon_id": taxon_id, "scientific_name": name,
                    "taxonomic_rank": rank, "column": column, "filled_with": "",
                    "peer_rows": sum(options.values()), "abundance_rows": "",
                    "note": "skipped: genus_key maps to "
                            + " | ".join(sorted(options)),
                })
                continue
            value, n_peers = next(iter(options.items()))
            counts[f"{column}:filled"] += 1
            updates.append((value, taxon_id))
            rows_out.append({
                "store": label, "taxon_id": taxon_id, "scientific_name": name,
                "taxonomic_rank": rank, "column": column, "filled_with": value,
                "peer_rows": n_peers, "abundance_rows": "",
                "note": "from the consensus of rows sharing this genus_key",
            })
        if not DRY_RUN and updates:
            cur.executemany(
                f"UPDATE taxa SET {column} = {placeholder} WHERE id = {placeholder}",
                updates)
    return counts


def report(label: str, cur) -> None:
    cur.execute("""SELECT COUNT(*) FROM taxa
                   WHERE phylum IS NULL OR TRIM(phylum) = ''""")
    print(f"   {label}: taxa still with no phylum: {cur.fetchone()[0]}")


def main() -> int:
    rows: list[dict] = []

    print("mysql:")
    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    report("mysql", cur)
    counts = fill_store("mysql", cur, "%s", rows)
    if not DRY_RUN:
        conn.commit()
    for key in sorted(counts):
        print(f"      {key:34} {counts[key]:>5}")
    report("mysql", cur)
    conn.close()

    if Path(SQLITE).exists():
        print("\nsqlite mirror:")
        sq = sqlite3.connect(SQLITE)
        sc = sq.cursor()
        report("sqlite", sc)
        sq_counts = fill_store("sqlite", sc, "?", rows)
        if not DRY_RUN:
            sq.commit()
        for key in sorted(sq_counts):
            print(f"      {key:34} {sq_counts[key]:>5}")
        report("sqlite", sc)
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
