#!/usr/bin/env python3
"""Record what rank NCBI gives each tax ID, beside the layer the row sits in.

taxonomic_rank in this database is the abundance LAYER a row was loaded into --
GMrepo publishes a genus array and a species array, and the loaders keep that
split because the per-(sample, rank) sum to 1 is defined over it. It is not a
claim about what the organism is, and for 152 rows the two differ:

    stored genus    -> NCBI family        90    Muribaculaceae, Tannerellaceae
    stored genus    -> NCBI order         19    Candidatus Gastranaerophilales
    stored genus    -> NCBI class          7    Candidatus Sericytochromatia, Mammalia
    stored genus    -> NCBI phylum         5    Candidatus Latescibacterota
    stored genus    -> other               7    species 4, subgenus, infraclass, infraorder
    stored species  -> NCBI strain        22    Bifidobacterium longum (1322347)
    stored species  -> NCBI subspecies      1    Streptomyces lavendulae (58340)
    stored species  -> other                1    species subgroup

The strain row said 277 when this was written and says 22 now, and the drop is
NCBI's, not a repair: this table compares the layer against the CURRENT dump, so
it moves when a dump lands. The stale thing was ncbi_rank in MySQL, which still
held the older answer for 257 rows long after the mirror had been refreshed --
see tests/test_cross_store_tax_ids.py, which is what noticed.

THE 123 COARSE GENUS-LAYER ROWS ARE NOT AN ERROR TO MOVE. GMrepo's own export
labels them taxonomic_rank = genus, so the database is faithful, and they are the
unclassified residual within a family rather than a duplicate of it: mean shares
run 0.17% (Erysipelotrichaceae) to 0.6% (Tannerellaceae) while the real genera
beside them carry far more, which is what a 16S classifier produces when it can
place a read to family but no further. 3,443 samples hold such a row alongside
genuine genera of the same family, and because the parts are disjoint the sum to
1 is correct. Moving them out of the genus layer would break that invariant and
throw away real mass.

So this records the fact and changes nothing else. What was missing was never the
rows; it was that nothing SAID they were coarser, so a model treated
Muribaculaceae as a peer of Bacteroides. v_abundance_genus now exposes
coarser_than_genus.

THE 283 SPECIES-LAYER ROWS ARE LEFT TO THE SESSIONS ALREADY ON THEM. A
species-named row carrying a strain tax ID is the `retaxid` case in
merge_duplicate_taxa.py's OVERRIDES, and two concurrent branches are working
through taxon identity right now -- the concatenated-genus repair and the gram
stain corrections both came out of it. Repointing tax IDs from here would collide
with that. They are recorded and reported, not touched.

Reads taxdump/nodes.dmp (override with GUTDB_TAXDUMP).
Idempotent. Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import csv
import os
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect

DRY_RUN = "--apply" not in sys.argv
TAXDUMP = Path(os.environ.get("GUTDB_TAXDUMP", "taxdump"))
AUDIT = Path("data/ncbi_rank_audit.csv")
# A dry run must not overwrite the record of what was last applied. Writing it
# unconditionally destroys the baseline the dry run exists to be compared
# against, which is the third place in this repo that has done it --
# adjudicate_conflicts.py and one concurrent session's script were the others.
# A tool that reports on state does not write the file that records the state.
DRY_RUN_AUDIT = AUDIT.with_name(AUDIT.stem + "_dry_run.csv")
AUDIT_FIELDS = ("taxon_id", "scientific_name", "ncbi_tax_id", "stored_rank",
                "ncbi_rank", "abundance_rows", "associations", "note")
# Ranks that are at or below genus, so a genus-layer row carrying one is not
# coarser than its layer. 'no rank' and 'clade' say nothing either way.
NOT_COARSER = {"genus", "subgenus", "species", "subspecies", "strain",
               "species group", "species subgroup", "isolate", "serotype",
               "no rank", "clade"}


def ranks_from_dump() -> dict[int, str]:
    nodes = TAXDUMP / "nodes.dmp"
    if not nodes.exists():
        print(f"   {nodes} not present; nothing to backfill")
        return {}
    out: dict[int, str] = {}
    with nodes.open(encoding="utf-8") as fh:
        for line in fh:
            parts = line.split("\t|\t")
            if len(parts) > 2:
                out[int(parts[0])] = parts[2].strip()
    return out


def main() -> int:
    ranks = ranks_from_dump()
    print(f"ranks read from the dump: {len(ranks)}")
    if not ranks:
        return 1

    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    cur.execute("""SELECT t.id, t.scientific_name, t.ncbi_tax_id, t.taxonomic_rank,
                          t.ncbi_rank,
                          (SELECT COUNT(*) FROM sample_taxon_abundances a
                           WHERE a.taxon_id = t.id),
                          (SELECT COUNT(*) FROM taxon_disease_associations x
                           WHERE x.taxon_id = t.id)
                   FROM taxa t WHERE t.ncbi_tax_id IS NOT NULL""")
    rows = cur.fetchall()

    updates: list[tuple[str, int]] = []
    disagree: list[dict] = []
    counts: Counter = Counter()
    for taxon_id, name, taxid, stored_rank, current, n_ab, n_ev in rows:
        actual = ranks.get(taxid)
        if not actual:
            counts["taxid_not_in_dump"] += 1
            continue
        if actual != current:
            updates.append((actual, taxon_id))
        if actual == stored_rank or actual in ("no rank", "clade"):
            counts["agrees"] += 1
            continue
        coarser = (stored_rank == "genus" and actual not in NOT_COARSER)
        counts[f"{stored_rank} -> {actual}"] += 1
        disagree.append({
            "taxon_id": taxon_id, "scientific_name": name, "ncbi_tax_id": taxid,
            "stored_rank": stored_rank, "ncbi_rank": actual,
            "abundance_rows": n_ab, "associations": n_ev,
            "note": ("coarser than its layer: unclassified residual, left in place"
                     if coarser else
                     "finer than its layer: tax ID work, left to the taxon-identity"
                     " sessions"),
        })

    print(f"\nrows with a tax ID: {len(rows)}; "
          f"not in this dump: {counts['taxid_not_in_dump']}; "
          f"agreeing: {counts['agrees']}")
    print("\ndisagreements:")
    for key, n in sorted(((k, v) for k, v in counts.items() if " -> " in k),
                         key=lambda kv: -kv[1]):
        print(f"   {key:36} {n:>5}")

    coarse = [d for d in disagree if d["note"].startswith("coarser")]
    fine = [d for d in disagree if d["note"].startswith("finer")]
    print(f"\ncoarser than their layer: {len(coarse)} rows, "
          f"{sum(d['abundance_rows'] for d in coarse):,} abundance rows, "
          f"{sum(d['associations'] for d in coarse)} associations")
    print(f"finer than their layer:   {len(fine)} rows, "
          f"{sum(d['abundance_rows'] for d in fine):,} abundance rows "
          f"(left alone)")

    if not DRY_RUN and updates:
        cur.executemany("UPDATE taxa SET ncbi_rank = %s WHERE id = %s", updates)
        conn.commit()
        print(f"\nncbi_rank written for {len(updates)} rows")

    destination = DRY_RUN_AUDIT if DRY_RUN else AUDIT
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=AUDIT_FIELDS)
        writer.writeheader()
        writer.writerows(sorted(disagree, key=lambda d: -d["abundance_rows"]))
    print(f"wrote {destination} ({len(disagree)} rows)")
    conn.close()

    backfill_mirror(ranks)

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    return 0


def backfill_mirror(ranks: dict[int, str]) -> int:
    """Fill ncbi_rank in the SQLite mirror too, adding the column if absent.

    The mirror had no ncbi_rank column at all until v_abundance_genus's rank
    flags needed it, which left the genus_rank_flags_inconsistent counter unable
    to run against this store -- reported as FAILED_TO_RUN rather than passing
    quietly, which is the one thing that design was for. Same one-store gap this
    project keeps producing, this time in a column rather than a repair.

    Only the fill is done here, not the disagreement analysis: the audit above
    reports MySQL's 411 rank disagreements and the mirror's would be the same
    organisms under different row ids, which is noise rather than information.
    """
    import sqlite3

    path = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")
    if not Path(path).exists():
        print(f"\nsqlite mirror {path} not present; skipped")
        return 0

    print("\nsqlite mirror:")
    conn = sqlite3.connect(path)
    cur = conn.cursor()
    cur.execute("PRAGMA table_info(taxa)")
    has_column = "ncbi_rank" in {row[1] for row in cur.fetchall()}
    if not has_column:
        if DRY_RUN:
            print("   would add the ncbi_rank column and fill it")
            conn.close()
            return 0
        cur.execute("ALTER TABLE taxa ADD COLUMN ncbi_rank TEXT")
        print("   added the ncbi_rank column")

    cur.execute("""SELECT id, ncbi_tax_id, ncbi_rank FROM taxa
                   WHERE ncbi_tax_id IS NOT NULL""")
    updates = []
    not_in_dump = 0
    for taxon_id, taxid, current in cur.fetchall():
        actual = ranks.get(int(taxid))
        if not actual:
            not_in_dump += 1
            continue
        if actual != current:
            updates.append((actual, taxon_id))
    print(f"   {len(updates)} rows to fill, {not_in_dump} tax IDs not in the dump")
    if not DRY_RUN and updates:
        cur.executemany("UPDATE taxa SET ncbi_rank = ? WHERE id = ?", updates)
        conn.commit()
        print(f"   ncbi_rank written for {len(updates)} rows")
    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
