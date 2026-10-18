#!/usr/bin/env python3
"""Take the misattributed bracket mass back out of MySQL's genus layer.

scripts/fix_genus_layer_purity.py repointed GMrepo's bracketed genus-array
entries onto their leading name token, so three organisms' genus-layer mass was
added to a genus they do not belong to. NCBI's brackets mean the name is
MISAPPLIED and assigns these three no genus at all -- they are incertae sedis,
and all three are Bacillota:

    84026   [Clostridium] methylpentosum   family Oscillospiraceae, no genus
    384638  [Bacteroides] pectinophilus    no family, no genus
    51123   [Eubacterium] saphenum         family Anaerovoracaceae, no genus

So a Firmicute's mass went under the Bacteroidetes genus Bacteroides. Measured
here: 218 samples averaging 1.5% of the genus composition and reaching 11.6%, in
the most-measured genus in this database.

THE AMOUNT IS MEASURED, NOT ASSUMED, and that is what makes this safe. The
repoint left no audit, so the script reconciles every affected sample against the
GMrepo source before touching it: MySQL's stored genus value must equal
source(parent) + source(bracket) to within 1e-6. All 1,215 affected samples
reconcile -- 218 Bacteroides, 953 Clostridium, 44 Eubacterium -- and any sample
that does not is skipped and reported rather than adjusted by a guessed amount.

The source files disagree for 682 of those samples: gmrepo_reload.csv holds
exactly twice the value in gmrepo_run_abundances.csv, a post-duplication
artifact. Rather than pick a file, the reconciliation tries every candidate value
and accepts the one that balances against what MySQL actually stores. The
database decides, not the CSV.

WHAT IT WRITES. For each affected sample: subtract the bracket value from the
parent genus row, and insert a genus-rank row for the bracket organism carrying
it. Mass is conserved exactly, so no rank group changes and no renormalisation is
needed -- the sum-to-1 counter is the end-to-end check. A parent row left at or
below zero is deleted rather than kept as a zero, which happens where the
bracket organism was the only contributor.

The new taxa mirror what the SQLite store already holds after
fix_mirror_bracket_names.py: scientific_name and genus both the full bracketed
string, rank genus, the NCBI tax ID, and the NCBI lineage above genus. `genus`
carries the bracketed string rather than being left empty because
validate's taxa_missing_genus counter requires every taxon to have one, and the
bracketed string is the honest name of this genus-layer feature.

Writes data/bracket_repoint_unwind_audit.csv, which doubles as the recovery
record the original repoint never left. Idempotent. Dry run by default.
"""
from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect

DRY_RUN = "--apply" not in sys.argv
AUDIT = Path("data/bracket_repoint_unwind_audit.csv")
AUDIT_FIELDS = ("run_accession", "sample_id", "bracket_name", "parent_genus",
                "stored_before", "bracket_value", "parent_after", "action")
SOURCES = ("data/gmrepo_all_abundances.csv", "data/gmrepo_run_abundances.csv",
           "data/gmrepo_reload.csv")
# bracket organism -> (parent token it was folded into, taxid, phylum, family)
BRACKETS = {
    "[Bacteroides] pectinophilus": ("Bacteroides", 384638, "Firmicutes", ""),
    "[Clostridium] methylpentosum": ("Clostridium", 84026, "Firmicutes",
                                     "Oscillospiraceae"),
    "[Eubacterium] saphenum": ("Eubacterium", 51123, "Firmicutes",
                               "Anaerovoracaceae"),
}
TOL = 1e-6


def load_source() -> dict[str, dict[str, set]]:
    """run -> name -> {candidate values}, for the brackets and their parents."""
    wanted = set(BRACKETS) | {p for p, *_ in BRACKETS.values()}
    out: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))
    for path in SOURCES:
        try:
            fh = open(path, newline="")
        except FileNotFoundError:
            continue
        with fh:
            for row in csv.DictReader(fh):
                name = (row.get("scientific_name") or "").strip()
                if name in wanted and row.get("taxonomic_rank") == "genus":
                    run = row.get("run_id") or row.get("run_accession")
                    out[run][name].add(float(row["relative_abundance"]))
    return out


def already_unwound(cur, sample_id: int, bracket: str) -> bool:
    """True when this sample already carries a row for the bracket organism."""
    cur.execute("""SELECT 1 FROM sample_taxon_abundances a
                   JOIN taxa t ON t.id = a.taxon_id
                   WHERE a.sample_id = %s AND t.scientific_name = %s
                     AND t.taxonomic_rank = 'genus' LIMIT 1""",
                (sample_id, bracket))
    return cur.fetchone() is not None


def ensure_bracket_taxon(cur, name: str) -> int:
    parent, taxid, phylum, family = BRACKETS[name]
    cur.execute("""SELECT id FROM taxa
                   WHERE scientific_name = %s AND taxonomic_rank = 'genus'""",
                (name,))
    row = cur.fetchone()
    if row:
        return row[0]
    cur.execute(
        """INSERT INTO taxa (superkingdom, phylum, family, genus, species,
                             genus_key, species_key, scientific_name,
                             taxonomic_rank, ncbi_tax_id, ncbi_rank,
                             source_database)
           VALUES ('Bacteria', %s, %s, %s, '', %s, '', %s, 'genus', %s,
                   'species', 'GMrepo')""",
        (phylum, family or None, name, name.casefold(), name, taxid))
    return cur.lastrowid


def main() -> int:
    source = load_source()
    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()

    cur.execute("""SELECT s.run_accession, s.id, t.scientific_name, t.id,
                          a.relative_abundance
                   FROM sample_taxon_abundances a
                   JOIN taxa t ON t.id = a.taxon_id AND t.taxonomic_rank = 'genus'
                   JOIN samples s ON s.id = a.sample_id
                   WHERE t.scientific_name IN (%s, %s, %s)""",
                tuple(p for p, *_ in BRACKETS.values()))
    stored: dict[tuple[str, str], tuple[int, int, float]] = {}
    for run, sample_id, parent, taxon_id, value in cur.fetchall():
        stored[(run, parent)] = (sample_id, taxon_id, float(value))

    rows_out: list[dict] = []
    planned = skipped = done = 0
    for bracket, (parent, *_rest) in BRACKETS.items():
        for run, names in source.items():
            if bracket not in names:
                continue
            key = (run, parent)
            if key not in stored:
                continue
            sample_id, parent_taxon, got = stored[key]
            # Reconcile: which candidate bracket value balances the books?
            match = None
            for p in names[bracket]:
                for b in (names.get(parent) or {0.0}):
                    if abs(got - (b + p)) < TOL:
                        match = p
                        break
                if match is not None:
                    break
            if match is None:
                # Already unwound is the common reason the books no longer
                # balance, and it must not read as a failure: once the mass has
                # moved, stored(parent) is deliberately less than
                # source(parent) + source(bracket). Check for the bracket row
                # before calling it a skip.
                if already_unwound(cur, sample_id, bracket):
                    done += 1
                    continue
                skipped += 1
                rows_out.append({
                    "run_accession": run, "sample_id": sample_id,
                    "bracket_name": bracket, "parent_genus": parent,
                    "stored_before": got, "bracket_value": "",
                    "parent_after": "", "action": "skipped: does not reconcile",
                })
                continue
            after = got - match
            action = "subtract_and_split" if after > TOL else "move_whole_row"
            planned += 1
            rows_out.append({
                "run_accession": run, "sample_id": sample_id,
                "bracket_name": bracket, "parent_genus": parent,
                "stored_before": got, "bracket_value": match,
                "parent_after": round(after, 10), "action": action,
            })
            if DRY_RUN:
                continue
            bracket_taxon = ensure_bracket_taxon(cur, bracket)
            if after > TOL:
                cur.execute("""UPDATE sample_taxon_abundances
                               SET relative_abundance = %s
                               WHERE sample_id = %s AND taxon_id = %s""",
                            (after, sample_id, parent_taxon))
            else:
                cur.execute("""DELETE FROM sample_taxon_abundances
                               WHERE sample_id = %s AND taxon_id = %s""",
                            (sample_id, parent_taxon))
            cur.execute("""INSERT INTO sample_taxon_abundances
                             (sample_id, taxon_id, relative_abundance)
                           VALUES (%s, %s, %s)
                           ON DUPLICATE KEY UPDATE
                             relative_abundance = relative_abundance + VALUES(relative_abundance)""",
                        (sample_id, bracket_taxon, match))

    if not DRY_RUN:
        conn.commit()

    by_action: dict[str, int] = defaultdict(int)
    moved = 0.0
    for r in rows_out:
        by_action[r["action"]] += 1
        if r["bracket_value"] != "":
            moved += float(r["bracket_value"])
    print(f"affected samples: {planned}   already unwound: {done}   "
          f"skipped: {skipped}")
    for action, n in sorted(by_action.items()):
        print(f"   {action:34} {n}")
    print(f"mass moved out of the parent genera: {moved:.4f}")

    cur.execute("""SELECT COUNT(*) FROM (
                     SELECT a.sample_id, t.taxonomic_rank,
                            SUM(a.relative_abundance) total
                     FROM sample_taxon_abundances a
                     JOIN taxa t ON t.id = a.taxon_id
                     GROUP BY a.sample_id, t.taxonomic_rank
                     HAVING ABS(total - 1.0) > 0.01) bad""")
    print(f"rank groups not summing to 1: {cur.fetchone()[0]} "
          f"(must be 0 -- mass is conserved)")
    conn.close()

    # NEVER overwrite a non-empty audit with nothing. This file records what was
    # DONE -- it is the recovery record the original repoint never left -- so a
    # clean re-run finding nothing to do must leave it alone. The opposite rule
    # applies to data/phylum_reclassification_audit.csv, which records what is
    # OUTSTANDING and is truncated when a run comes back clean; a re-run of this
    # script wiped 1,215 rows before the distinction was drawn.
    AUDIT.parent.mkdir(parents=True, exist_ok=True)
    if rows_out:
        with AUDIT.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=AUDIT_FIELDS)
            writer.writeheader()
            writer.writerows(rows_out)
        print(f"wrote {AUDIT} ({len(rows_out)} rows) -- the recovery record the "
              f"original repoint never left")
    elif AUDIT.exists():
        print(f"nothing to record; {AUDIT} left as it is")

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    return 1 if skipped else 0


if __name__ == "__main__":
    raise SystemExit(main())
