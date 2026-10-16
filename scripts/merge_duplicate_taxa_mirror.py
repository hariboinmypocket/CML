#!/usr/bin/env python3
"""Merge the mirror's duplicate tax IDs, reusing the MySQL script's decisions.

scripts/merge_duplicate_taxa.py fixed this for MySQL and is MySQL-only, so the
mirror still holds the pre-merge state: 104 tax IDs spread over 208 rows. They
are the species-level form of the Firmicutes/Bacillota split already repaired at
phylum, class, order and family -- same tax ID, two names, one organism:

    33038   Ruminococcus gnavus   || Mediterraneibacter gnavus   8,363 abundances
    40520   Blautia obeum         || Ruminococcus obeum          8,029
    39485   Eubacterium eligen    || Lachnospira eligens         7,448
    165179  Prevotella copri      || Segatella copri             4,688

AND 582 (sample, tax ID) CELLS COUNT ONE ORGANISM TWICE, across 12,392 of the
mirror's 23,494 samples. The per-rank sums are still exactly 1.0, because
renormalisation happened after the duplication -- so the sum-to-1 counter cannot
see it, and none of validate's 47 counters did. The duplicated organism holds
double share while every other taxon in those samples is squeezed.

THE DECISION LOGIC IS IMPORTED, NOT REIMPLEMENTED. norm(), load_groups(), plan()
and OVERRIDES come from the MySQL script, which is the whole point: a second
canonicalisation rule would be a second answer to "which name is this organism",
and the two stores would disagree about it the way they disagreed about
everything else this sweep has touched. Only the write phase is rewritten, for
SQLite's placeholders and its own id space.

NCBI's current names come from the local taxdump rather than the API the MySQL
script calls, so this needs no network and no NCBI_EMAIL. Same data, same field.

Anything plan() declines -- an unreviewed epithet conflict, or a tax ID whose
OVERRIDES entry is not a plain merge -- is reported and left alone. The OVERRIDES
actions that name MySQL row ids ('retaxid', 'delete_and_retaxid') cannot be
applied here, because those ids mean a different organism in this store; they are
skipped explicitly rather than translated by guesswork.

Writes data/taxa_merge_mirror_audit.csv. Idempotent. Dry run by default.
"""
from __future__ import annotations

import csv
import importlib.util
import os
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

_spec = importlib.util.spec_from_file_location(
    "merge_duplicate_taxa", ROOT / "scripts" / "merge_duplicate_taxa.py")
_mysql_merge = importlib.util.module_from_spec(_spec)
_argv = sys.argv
sys.argv = [_argv[0]]              # the module reads sys.argv for its own DRY_RUN
_spec.loader.exec_module(_mysql_merge)
sys.argv = _argv

norm = _mysql_merge.norm
load_groups = _mysql_merge.load_groups
plan = _mysql_merge.plan
OVERRIDES = _mysql_merge.OVERRIDES
TRAIT_COLUMNS = getattr(_mysql_merge, "TRAIT_COLUMNS", ())
NA_COLUMNS = getattr(_mysql_merge, "NA_COLUMNS", ())

DRY_RUN = "--apply" not in sys.argv
SQLITE = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")
TAXDUMP = Path(os.environ.get("GUTDB_TAXDUMP", "taxdump"))
AUDIT = Path("data/taxa_merge_mirror_audit.csv")
AUDIT_FIELDS = ("taxid", "canonical_id", "canonical_name", "loser_ids",
                "loser_names", "abundances_moved", "abundances_summed",
                "associations_moved", "rename_to", "note")


def ncbi_names_from_dump(taxids: set[int]) -> dict[str, dict[str, str]]:
    """{taxid: {'_scientific_name': name}} -- the one field plan() reads."""
    names_file = TAXDUMP / "names.dmp"
    if not names_file.exists():
        raise SystemExit(f"taxdump not found at {TAXDUMP}/; set GUTDB_TAXDUMP")
    wanted = {str(t) for t in taxids}
    out: dict[str, dict[str, str]] = {}
    with names_file.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.split("\t|")
            if len(parts) < 4 or parts[3].strip() != "scientific name":
                continue
            taxid = parts[0].strip()
            if taxid in wanted:
                out[taxid] = {"_scientific_name": parts[1].strip()}
    return out


def touched_samples_for(cur, canonical_id: int, losers: list[dict]) -> set[int]:
    """Samples whose rank group will lose mass, collected BEFORE the delete."""
    loser_ids = [r["id"] for r in losers]
    marks = ",".join("?" * len(loser_ids))
    cur.execute(f"""SELECT DISTINCT l.sample_id FROM sample_taxon_abundances l
                    WHERE l.taxon_id IN ({marks})
                      AND EXISTS (SELECT 1 FROM sample_taxon_abundances c
                                  WHERE c.taxon_id = ? AND c.sample_id = l.sample_id)""",
                (*loser_ids, canonical_id))
    return {r[0] for r in cur.fetchall()}


def renormalize(cur, sample_ids: set[int]) -> int:
    """Rescale the rank groups that lost a duplicated copy back to 1.0.

    Removing a double-counted row drops its group below 1.0, and that is correct
    arithmetic rather than damage: the stored values were normalised WHILE the
    duplicate was present, so the duplicated organism held double share and every
    other taxon in the sample was squeezed. Rescaling restores the true
    proportions. The MySQL script does exactly this; omitting it here left 49
    groups summing to about 0.99 and is why this function exists.

    Only touched groups are rescaled, so clean samples are never perturbed by
    float noise.
    """
    rescaled = 0
    for sample_id in sorted(sample_ids):
        for rank in ("genus", "species"):
            cur.execute("""SELECT SUM(a.relative_abundance)
                           FROM sample_taxon_abundances a
                           JOIN taxa t ON t.id = a.taxon_id
                           WHERE a.sample_id = ? AND t.taxonomic_rank = ?""",
                        (sample_id, rank))
            total = cur.fetchone()[0]
            if not total or abs(float(total) - 1.0) < 1e-9:
                continue
            cur.execute("""UPDATE sample_taxon_abundances
                           SET relative_abundance = relative_abundance / ?
                           WHERE sample_id = ?
                             AND taxon_id IN (SELECT id FROM taxa
                                              WHERE taxonomic_rank = ?)""",
                        (float(total), sample_id, rank))
            rescaled += 1
    return rescaled


def merge_one(cur, canonical_id: int, losers: list[dict]) -> tuple[int, int, int]:
    """Repoint a group's abundances and associations onto the canonical row."""
    loser_ids = [r["id"] for r in losers]
    marks = ",".join("?" * len(loser_ids))

    # Abundances: where the canonical row is ALREADY in the sample, the two rows
    # are one measurement stored twice, so the loser's row is dropped rather than
    # added. Adding would double the organism's share; the MySQL script verified
    # every overlapping pair holds an identical value.
    cur.execute(f"""SELECT COUNT(*) FROM sample_taxon_abundances l
                    WHERE l.taxon_id IN ({marks})
                      AND EXISTS (SELECT 1 FROM sample_taxon_abundances c
                                  WHERE c.taxon_id = ? AND c.sample_id = l.sample_id)""",
                (*loser_ids, canonical_id))
    overlapping = cur.fetchone()[0]
    cur.execute(f"""DELETE FROM sample_taxon_abundances
                    WHERE taxon_id IN ({marks})
                      AND sample_id IN (SELECT sample_id FROM sample_taxon_abundances
                                        WHERE taxon_id = ?)""",
                (*loser_ids, canonical_id))
    cur.execute(f"""UPDATE sample_taxon_abundances SET taxon_id = ?
                    WHERE taxon_id IN ({marks})""", (canonical_id, *loser_ids))
    moved = cur.rowcount

    cur.execute(f"""DELETE FROM taxon_disease_associations
                    WHERE taxon_id IN ({marks})
                      AND (disease_id, comparison_id) IN (
                          SELECT disease_id, comparison_id
                          FROM taxon_disease_associations WHERE taxon_id = ?)""",
                (*loser_ids, canonical_id))
    cur.execute(f"""UPDATE taxon_disease_associations SET taxon_id = ?
                    WHERE taxon_id IN ({marks})""", (canonical_id, *loser_ids))
    assoc_moved = cur.rowcount

    cur.execute(f"DELETE FROM taxa WHERE id IN ({marks})", loser_ids)
    return moved, overlapping, assoc_moved


def main() -> int:
    if not Path(SQLITE).exists():
        print(f"sqlite mirror {SQLITE} not present; nothing to do")
        return 0

    conn = sqlite3.connect(SQLITE)
    cur = conn.cursor()
    groups = load_groups(cur)
    if not groups:
        print("no duplicate tax IDs in the mirror; nothing to do")
        conn.close()
        return 0
    print(f"duplicate tax IDs in the mirror: {len(groups)} "
          f"over {sum(len(v) for v in groups.values())} rows")

    ncbi = ncbi_names_from_dump(set(groups))
    missing = [t for t in groups if str(t) not in ncbi]
    if missing:
        print(f"   {len(missing)} tax IDs are not in the dump and are skipped: "
              f"{missing[:5]}")
    merges, skipped = plan(groups, ncbi)

    # A GROUP SPANNING TWO RANKS IS NEVER MERGEABLE, whatever the tax ID says.
    # The abundance layer is defined per rank and the sum-to-1 invariant holds
    # per (sample, rank), so repointing a genus-layer row onto a species-layer
    # canonical moves mass from one composition into another. Three groups here
    # are the bracket-notation artifacts -- 'Eubacterium saphenum' beside
    # 'Eubacterium] saphenum', the opening bracket lost -- with the mangled row
    # sitting at genus rank and carrying genus-layer abundance.
    #
    # Merging them created one-row species groups in 663 samples, each summing
    # to about 1e-04 instead of 1.0, which is how this guard came to be written.
    # MySQL does not have these rows: fix_genus_layer_purity.py repointed its
    # genus-layer binomials to their parent genus, and the mirror never got that
    # repair, so the script that fixed MySQL never had to refuse them.
    mixed_rank = []
    for entry in list(merges):
        ranks = {entry["canonical"]["rank"]} | {r["rank"] for r in entry["losers"]}
        if len(ranks) > 1:
            merges.remove(entry)
            mixed_rank.append(entry)
            skipped.append({
                "taxid": entry["taxid"], "action": "spans_two_ranks",
                "note": " | ".join(f"{r['sci']} ({r['rank']})"
                                   for r in [entry["canonical"], *entry["losers"]]),
            })
    print(f"   plan: {len(merges)} merges, {len(skipped)} declined "
          f"({len(mixed_rank)} of them for spanning two ranks)")

    rows_out: list[dict] = []
    touched: set[int] = set()
    applied = 0
    for entry in merges:
        taxid = entry["taxid"]
        if str(taxid) not in ncbi:
            continue
        canon, losers = entry["canonical"], entry["losers"]
        moved = overlapping = assoc_moved = 0
        if not DRY_RUN:
            touched |= touched_samples_for(cur, canon["id"], losers)
            moved, overlapping, assoc_moved = merge_one(cur, canon["id"], losers)
            if entry["rename_to"]:
                cur.execute("UPDATE taxa SET scientific_name = ? WHERE id = ?",
                            (entry["rename_to"], canon["id"]))
        applied += 1
        rows_out.append({
            "taxid": taxid, "canonical_id": canon["id"],
            "canonical_name": canon["sci"],
            "loser_ids": " ".join(str(r["id"]) for r in losers),
            "loser_names": " || ".join(r["sci"] for r in losers),
            "abundances_moved": moved, "abundances_summed": overlapping,
            "associations_moved": assoc_moved,
            "rename_to": entry["rename_to"], "note": entry["note"],
        })

    for entry in skipped:
        action = entry.get("action", "")
        rows_out.append({
            "taxid": entry["taxid"], "canonical_id": "", "canonical_name": "",
            "loser_ids": "", "loser_names": entry.get("note", ""),
            "abundances_moved": "", "abundances_summed": "",
            "associations_moved": "", "rename_to": "",
            "note": f"declined: {action}" + (
                " (names a MySQL row id, meaningless here)"
                if action in ("retaxid", "delete_and_retaxid") else ""),
        })
        print(f"   declined {entry['taxid']}: {action}")

    if not DRY_RUN:
        rescaled = renormalize(cur, touched)
        conn.commit()
        print(f"\n{len(touched)} samples lost a duplicated copy; "
              f"{rescaled} rank groups rescaled back to 1.0")

    cur.execute("""SELECT COUNT(*) FROM (SELECT ncbi_tax_id FROM taxa
                   WHERE ncbi_tax_id IS NOT NULL
                   GROUP BY ncbi_tax_id HAVING COUNT(*) > 1) d""")
    print(f"\nduplicate tax IDs remaining: {cur.fetchone()[0]}")
    cur.execute("""SELECT COUNT(*) FROM (
                     SELECT a.sample_id, t.ncbi_tax_id
                     FROM sample_taxon_abundances a JOIN taxa t ON t.id = a.taxon_id
                     WHERE t.ncbi_tax_id IS NOT NULL
                     GROUP BY a.sample_id, t.ncbi_tax_id
                     HAVING COUNT(DISTINCT a.taxon_id) > 1) d""")
    print(f"(sample, tax ID) cells counted twice: {cur.fetchone()[0]}")
    conn.close()

    AUDIT.parent.mkdir(parents=True, exist_ok=True)
    with AUDIT.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=AUDIT_FIELDS)
        writer.writeheader()
        writer.writerows(rows_out)
    print(f"wrote {AUDIT} ({len(rows_out)} rows)")

    if DRY_RUN:
        print(f"\n*** DRY RUN - nothing written. {applied} merges planned. "
              f"Re-run with --apply ***")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
