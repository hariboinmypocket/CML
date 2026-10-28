#!/usr/bin/env python3
"""Reconcile stored class, order and family against the NCBI dump, in both stores.

`normalize_reclassified_phyla.py` reconciles the phylum and stops there, so the
three ranks between phylum and genus have never been checked against the dump.
Measured against each row's own `ncbi_tax_id`: 336 rows disagree on class, 533 on
order, 184 on family.

Almost none of that is corruption. It is the dump moving on, in two kinds that
must not be treated alike -- the distinction `normalize_reclassified_phyla.py`
already draws for the phylum, where a shared tax ID makes a rename that belongs
in PHYLUM_SYNONYMS and a different tax ID makes a reclassification that needs a
human in APPROVED:

  rename   The stored name is a synonym of the node NCBI now names differently:
           `Actinobacteria` and `Actinomycetes` are both class 1760 (255 rows),
           `Bacillales` and `Caryophanales` both order 1385 (108 rows). Nothing
           moved; one spelling is older. Folding these is safe in principle, but
           it rewrites hundreds of rows to adopt a name the rest of the table
           does not use, which is a vocabulary decision rather than a repair --
           so they are reported and left alone.

  moved    The stored name resolves to a different node of the same rank, i.e.
           the organism was reclassified: `Saccharomycetes` (4891) to
           `Pichiomycetes` (3316584), `Deltaproteobacteria` (28221) to
           `Desulfovibrionia` (2893115). Each is a taxonomic decision about one
           group, so it is surfaced per row and written only when listed in
           APPROVED_ROWS below.

WHAT THIS WRITES TODAY

Two rows, both found by the gram-stain check in the session that repaired the
concatenated-genus names. `Veillonella sp. 3 1 44` (1895) and `Veillonella sp.
6 1 27` (1897) store class Clostridia and order Clostridiales. Veillonella is
Negativicutes, and three things say so independently: the dump, for each row's
own tax ID (457416, 450749); every other Veillonella row in the table, all 12 of
them including the genus-rank row; and the rows' own `family` cell, which already
reads Veillonellaceae -- a family that cannot sit inside Clostridia. So these
rows are internally inconsistent rather than merely dated, which is what sets
them apart from the other 334 class disagreements.

Their gram_stain of 'negative' is correct and is left alone. Veillonella is
genuinely gram-negative; the stored class was the error, and the gram check only
surfaced it because the two happened to coincide. The 18 other gram-negative
Clostridia rows are real staining observations on rows NCBI does place in
Clostridia -- Roseburia, Lachnospira eligens, Clostridium symbiosum -- several
carrying thousands of abundance rows, and nothing here touches them.

Both rows read order `Eubacteriales` rather than `Clostridiales`: the synonym
fold in `fold_phylum_synonyms.py` moved them there while this entry was being
written, and the guard below refused to act until the expectation was re-measured
-- which is the guard doing its job, on a database that changes hourly.

Keyed on `(genus_key, species_key)` with the stored values to expect, so a row
reclassified or renamed since this list was written is skipped rather than
rewritten on stale grounds. Re-running is a no-op.

Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import csv
import os
import sqlite3
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import require_local
from gutdb.config import Settings
from gutdb.db import connect
from gutdb.transform import clean_str

DRY_RUN = "--apply" not in sys.argv
TAXDUMP = os.environ.get("GUTDB_TAXDUMP", "taxdump")
SQLITE_PATH = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")
AUDIT = Path("data/lineage_rank_audit.csv")

# taxa column -> NCBI rank
RANKS = (("class_name", "class"), ("order_name", "order"), ("family", "family"))

# Rows this script may rewrite, keyed by what identifies the organism rather
# than by taxa.id -- the SQLite mirror has its own id space, and a key cannot
# silently land on a different row.
APPROVED_ROWS: dict[tuple[str, str], dict] = {
    ("veillonella", "sp. 3 1 44"): {
        "expect": {"class_name": "Clostridia", "order_name": "Eubacteriales",
                   "family": "Veillonellaceae"},
        "set": {"class_name": "Negativicutes", "order_name": "Veillonellales"},
        "note": "Veillonella is Negativicutes: the dump for tax 457416, the "
                "other 12 Veillonella rows, and this row's own family cell all "
                "say so; family Veillonellaceae cannot sit inside Clostridia",
    },
    ("veillonella", "sp. 6 1 27"): {
        "expect": {"class_name": "Clostridia", "order_name": "Eubacteriales",
                   "family": "Veillonellaceae"},
        "set": {"class_name": "Negativicutes", "order_name": "Veillonellales"},
        "note": "same as Veillonella sp. 3 1 44, for tax 450749",
    },
}


def load_taxdump() -> tuple[dict[str, int], dict[int, str], dict[int, dict[str, str]]]:
    """Return (name -> taxid, taxid -> rank, taxid -> {rank: ancestor name}).

    Parsed the way normalize_reclassified_phyla.py does it, one rank wider:
    names.dmp for identity including synonyms, nodes.dmp for the parent chain,
    then names.dmp again for the ancestor names at the three ranks of interest.
    """
    names_path, nodes_path = f"{TAXDUMP}/names.dmp", f"{TAXDUMP}/nodes.dmp"
    if not (os.path.exists(names_path) and os.path.exists(nodes_path)):
        raise SystemExit(
            f"taxdump not found at {TAXDUMP}/. Unzip taxdmp.zip there or set GUTDB_TAXDUMP.")

    # exists() passes for a file iCloud has evicted, whose first read then
    # blocks instead of failing; these are the repo's two largest files.
    names_path = str(require_local(names_path))
    nodes_path = str(require_local(nodes_path))

    wanted = {rank for _, rank in RANKS}
    parent: dict[int, int] = {}
    rank_of: dict[int, str] = {}
    with open(nodes_path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split("\t|")
            taxid, par, rank = int(parts[0]), int(parts[1]), parts[2].strip()
            parent[taxid] = par
            if rank in wanted:
                rank_of[taxid] = rank

    # Names of nodes AT these ranks, indexed by rank so a homonym cannot answer
    # for the wrong layer: the bare string "Actinobacteria" is a synonym of the
    # *phylum* 201174, and matching it without the rank is how 255 rows would be
    # read as reclassified when the class simply has two names.
    #
    # Indexed on the first token, because a name above genus is a single word
    # and the dump carries the superseded ones with their authority attached --
    # class 1760 is "Actinomycetes" with "Actinobacteria Stackebrandt et al.
    # 1997" among its synonyms, so a whole-string match finds neither.
    node_name: dict[int, str] = {}
    rank_names: dict[tuple[str, str], set[int]] = {}
    with open(names_path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split("\t|")
            if len(parts) < 4:
                continue
            taxid = int(parts[0])
            rank = rank_of.get(taxid)
            if not rank:
                continue
            name, name_class = parts[1].strip(), parts[3].strip()
            if name_class == "scientific name":
                node_name[taxid] = name
            if name_class in ("scientific name", "synonym", "equivalent name",
                              "includes", "authority"):
                first = name.split()[0].casefold() if name.split() else ""
                if first:
                    rank_names.setdefault((rank, first), set()).add(taxid)

    lineage: dict[int, dict[str, str]] = {}

    def resolve(taxid: int) -> dict[str, str]:
        if taxid in lineage:
            return lineage[taxid]
        found: dict[str, str] = {}
        node, guard = taxid, 0
        while node and node != 1 and guard < 80:
            rank = rank_of.get(node)
            # The deepest node of a rank wins, so only fill an empty slot.
            if rank and rank not in found:
                found[rank] = node_name.get(node, "")
            node = parent.get(node, 0)
            guard += 1
        lineage[taxid] = found
        return found

    class Lineages:
        def get(self, taxid: int) -> dict[str, str]:
            return resolve(taxid) if taxid else {}

    return rank_names, rank_of, Lineages()


def classify(stored: str, ncbi: str, rank: str,
             rank_names: dict[tuple[str, str], set[int]]) -> str:
    """Why a stored rank name differs from the dump's, by the tax ID identity test.

    Two names resolving to the same node are one taxon under two spellings; two
    nodes are a reclassification. `rename` is expected here and is never written:
    fold_phylum_synonyms.py deliberately folds these onto the MAJORITY STORED
    name rather than NCBI's scientific name, because NCBI prefers names this
    database does not use (Kitasatosporales over Streptomycetales) and some carry
    a nomenclatural code the schema strips. So a stored name that differs from
    the dump while naming the same node is upstream policy, not a defect.
    """
    stored_nodes = rank_names.get((rank, stored.split()[0].casefold()), set())
    ncbi_nodes = rank_names.get((rank, ncbi.split()[0].casefold()), set())
    if not stored_nodes:
        return "unresolvable"
    if stored_nodes & ncbi_nodes:
        return "rename"
    return "moved"


def run_store(store: str, cursor, rank_names, rank_of, lineages,
              audit: list[dict]) -> dict[str, int]:
    placeholder = "%s" if store == "mysql" else "?"
    cursor.execute(
        "SELECT id, scientific_name, genus_key, species_key, ncbi_tax_id, "
        "class_name, order_name, family FROM taxa")
    rows = cursor.fetchall()

    stats = Counter()
    pairs: dict[tuple[str, str, str], int] = Counter()
    written = 0

    for (taxon_id, name, genus_key, species_key, taxid,
         class_name, order_name, family) in rows:
        stored = {"class_name": class_name, "order_name": order_name, "family": family}
        if not taxid:
            stats["no_tax_id"] += 1
            continue
        lineage = lineages.get(int(taxid))
        if not lineage:
            stats["unresolvable_tax_id"] += 1
            continue

        disagreements = {}
        for column, rank in RANKS:
            have, want = clean_str(stored[column]), clean_str(lineage.get(rank, ""))
            if not have or not want or have == want:
                stats["agrees" if have and want else "not_typed"] += 1
                continue
            kind = classify(have, want, rank, rank_names)
            stats[f"{rank}_{kind}"] += 1
            pairs[(rank, have, want)] += 1
            disagreements[column] = (have, want, kind)

        approved = APPROVED_ROWS.get((genus_key, species_key))
        if not disagreements or not approved:
            continue

        mismatched = [c for c, v in approved["expect"].items()
                      if clean_str(stored[c]) != v]
        if mismatched:
            print(f"    SKIPPED {taxon_id} {name!r}: expected "
                  + ", ".join(f"{c}={approved['expect'][c]!r}" for c in mismatched)
                  + ", found "
                  + ", ".join(f"{c}={clean_str(stored[c])!r}" for c in mismatched)
                  + " -- the row has changed since this entry was written")
            stats["approved_skipped"] += 1
            continue

        print(f"    {taxon_id:>8} {name[:34]:<34} "
              + "  ".join(f"{c}: {clean_str(stored[c])} -> {v}"
                          for c, v in approved["set"].items()))
        if not DRY_RUN:
            assignments = ", ".join(f"{c}={placeholder}" for c in approved["set"])
            cursor.execute(
                f"UPDATE taxa SET {assignments} WHERE id={placeholder}",
                tuple(approved["set"].values()) + (taxon_id,))
        written += 1
        audit.append({
            "store": store, "taxon_id": taxon_id, "scientific_name": name,
            "ncbi_tax_id": taxid,
            "changes": "; ".join(f"{c}: {clean_str(stored[c])} -> {v}"
                                 for c, v in approved["set"].items()),
            "note": approved["note"],
        })

    stats["approved_written"] = written
    for (rank, have, want), count in sorted(
            pairs.items(), key=lambda item: -item[1])[:10]:
        kind = classify(have, want, rank, rank_names)
        print(f"      {rank:<7} {have:<24} -> {want:<24} {count:>4}  {kind}")
    return stats


def main() -> int:
    print(f"loading taxdump from {TAXDUMP}/ ...")
    rank_names, rank_of, lineages = load_taxdump()
    print(f"  nodes at class/order/family {len(rank_of):,}   names indexed {len(rank_names):,}")

    audit: list[dict] = []
    totals = {}

    print("\nmysql")
    conn = connect(Settings.from_env(".env"))
    cursor = conn.cursor()
    totals["mysql"] = run_store("mysql", cursor, rank_names, rank_of, lineages, audit)
    if not DRY_RUN:
        conn.commit()
    cursor.close()
    conn.close()

    print("\nsqlite")
    if Path(SQLITE_PATH).exists():
        lite = sqlite3.connect(SQLITE_PATH)
        cursor = lite.cursor()
        totals["sqlite"] = run_store("sqlite", cursor, rank_names, rank_of,
                                     lineages, audit)
        if not DRY_RUN:
            lite.commit()
        cursor.close()
        lite.close()
    else:
        print(f"  (no SQLite mirror at {SQLITE_PATH}; skipped)")

    if audit:
        AUDIT.parent.mkdir(parents=True, exist_ok=True)
        with AUDIT.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(audit[0].keys()))
            writer.writeheader()
            writer.writerows(audit)
        print(f"\naudit written to {AUDIT}")

    for store, stats in totals.items():
        print(f"\n{store}")
        for key in sorted(stats):
            if key in ("agrees", "not_typed"):
                continue
            print(f"  {key:30} {stats[key]}")

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
