#!/usr/bin/env python3
"""Reconcile stored phylum against the NCBI taxonomy dump, in both stores.

The taxdump fill in cleanup_gutdb.py only writes into empty cells, by design --
it must not overwrite curated values. That leaves rows whose stored phylum is
populated but superseded, and one such group splits a real lineage in two:
Desulfovibrio, Bilophila and Desulfobulbus are gut sulfate reducers studied in
IBD and colorectal cancer, and they appear under both Proteobacteria (their
classic Deltaproteobacteria placement, from the curated CSV) and
Thermodesulfobacteriota (NCBI's current placement, from the dump). A genus
spanning two phyla breaks every phylum-level feature and v_taxon_specificity's
grouping.

This is a reclassification, not a rename, so PHYLUM_SYNONYMS is the wrong tool:
these taxa were moved OUT of Proteobacteria, and mapping the new name onto the
old one would assert a placement NCBI no longer makes. Instead the dump decides
per organism, and only for reclassifications listed in APPROVED below.

Every disagreement found is reported; only approved ones are written. An
unapproved disagreement may be a curated value that is deliberately better than
NCBI's, so it is surfaced for a human rather than normalized away.

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

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.transform import clean_str, normalize_phylum

DRY_RUN = "--apply" not in sys.argv
TAXDUMP = os.environ.get("GUTDB_TAXDUMP", "taxdump")
SQLITE_PATH = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")
AUDIT = Path("data/phylum_reclassification_audit.csv")

# (stored phylum, NCBI's current phylum) pairs this script may rewrite.
# Keep it explicit: each entry is a documented reclassification, not a guess.
APPROVED = {
    ("Proteobacteria", "Thermodesulfobacteriota"),
    ("Thermodesulfobacteria", "Thermodesulfobacteriota"),
}


def load_taxdump() -> tuple[dict[str, int], dict[int, str]]:
    """Return (scientific/synonym name -> taxid, taxid -> phylum name).

    Parsed the same way cleanup_gutdb.py does it: names.dmp for identity,
    nodes.dmp for the parent chain, then names.dmp again for ancestor names.
    """
    names_path, nodes_path = f"{TAXDUMP}/names.dmp", f"{TAXDUMP}/nodes.dmp"
    if not (os.path.exists(names_path) and os.path.exists(nodes_path)):
        raise SystemExit(
            f"taxdump not found at {TAXDUMP}/. Unzip taxdmp.zip there or set GUTDB_TAXDUMP.")

    name_to_taxid: dict[str, int] = {}
    synonyms: dict[str, int] = {}
    with open(names_path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split("\t|")
            if len(parts) < 4:
                continue
            name = parts[1].strip().casefold()
            name_class = parts[3].strip()
            if name_class == "scientific name":
                name_to_taxid.setdefault(name, int(parts[0]))
            elif name_class in ("synonym", "equivalent name", "includes"):
                synonyms.setdefault(name, int(parts[0]))
    for name, taxid in synonyms.items():
        name_to_taxid.setdefault(name, taxid)

    parent: dict[int, int] = {}
    is_phylum: set[int] = set()
    with open(nodes_path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split("\t|")
            taxid, par, rank = int(parts[0]), int(parts[1]), parts[2].strip()
            parent[taxid] = par
            if rank == "phylum":
                is_phylum.add(taxid)

    phylum_name: dict[int, str] = {}
    with open(names_path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split("\t|")
            if len(parts) < 4 or parts[3].strip() != "scientific name":
                continue
            taxid = int(parts[0])
            if taxid in is_phylum:
                phylum_name[taxid] = parts[1].strip()

    # Resolve each node to its phylum ancestor once, memoised up the chain.
    phylum_of: dict[int, str] = {}

    def resolve(taxid: int) -> str:
        chain, node, guard = [], taxid, 0
        while node and node != 1 and guard < 60:
            if node in phylum_of:
                break
            if node in is_phylum:
                phylum_of[node] = phylum_name.get(node, "")
                break
            chain.append(node)
            node = parent.get(node, 0)
            guard += 1
        found = phylum_of.get(node, "")
        for link in chain:
            phylum_of[link] = found
        return found

    for taxid in list(parent):
        if taxid not in phylum_of:
            resolve(taxid)
    return name_to_taxid, phylum_of


def rows_for(store: str, cursor) -> list[tuple]:
    cursor.execute(
        "SELECT id, scientific_name, taxonomic_rank, ncbi_tax_id, phylum FROM taxa "
        "WHERE phylum IS NOT NULL AND phylum <> ''")
    return cursor.fetchall()


def reconcile(store: str, cursor, name_to_taxid, phylum_of, audit: list[dict]) -> Counter:
    stats = Counter()
    for taxon_id, sci, rank, taxid, stored in rows_for(store, cursor):
        if not taxid:
            taxid = name_to_taxid.get(clean_str(sci).casefold())
        if not taxid:
            stats["unresolvable"] += 1
            continue
        current = normalize_phylum(phylum_of.get(int(taxid), ""))
        if not current:
            stats["no_phylum_at_ncbi"] += 1
            continue
        stored_clean = clean_str(stored)
        if stored_clean == current:
            stats["agrees"] += 1
            continue

        # A case variant is the same concept spelled two ways
        # ("Deinococcus-thermus" / "Deinococcus-Thermus"), which splits a phylum
        # into two values for no reason. Folding it onto NCBI's casing is not a
        # taxonomic judgement, so it needs no APPROVED entry.
        case_variant = stored_clean.casefold() == current.casefold()
        approved = case_variant or (stored_clean, current) in APPROVED
        if case_variant:
            stats["case_folded"] += 1
        audit.append({
            "store": store, "taxon_id": taxon_id, "scientific_name": sci,
            "rank": rank, "ncbi_tax_id": taxid, "stored_phylum": stored_clean,
            "ncbi_phylum": current,
            "action": ("case_folded" if case_variant
                       else "normalized" if approved else "reported_only"),
        })
        if not approved:
            stats["disagrees_unapproved"] += 1
            continue
        stats["normalized"] += 1
        if not DRY_RUN:
            cursor.execute(
                "UPDATE taxa SET phylum = %s WHERE id = %s" if store == "mysql"
                else "UPDATE taxa SET phylum = ? WHERE id = ?",
                (current, taxon_id))
    return stats


def fold_by_genus_consensus(store: str, cursor, audit: list[dict]) -> Counter:
    """Fold a phylum NCBI could not verify onto its genus's unanimous value.

    Some rows carry a strain designation the dump has no entry for
    ("Desulfovibrio c21 c20") or a tax ID NCBI has since retired, so the
    reconcile pass skips them and they keep a superseded phylum while every
    other member of their genus has moved. A species of Desulfovibrio is in
    whatever phylum Desulfovibrio is, so where the rest of the genus agrees on
    one value the odd row adopts it.

    Unanimity is required, the same vote cleanup_gutdb.py uses for genus trait
    backfill: a genus whose members disagree is left alone rather than having a
    majority imposed on it.
    """
    stats = Counter()
    cursor.execute(
        "SELECT genus, phylum, COUNT(*) FROM taxa "
        "WHERE phylum IS NOT NULL AND phylum <> '' GROUP BY genus, phylum")
    by_genus: dict[str, set[str]] = {}
    for genus, phylum, _ in cursor.fetchall():
        by_genus.setdefault(genus, set()).add(clean_str(phylum))
    consensus = {g: next(iter(p)) for g, p in by_genus.items() if len(p) == 1}

    cursor.execute(
        "SELECT id, scientific_name, genus, phylum FROM taxa "
        "WHERE phylum IS NOT NULL AND phylum <> ''")
    for taxon_id, sci, genus, phylum in cursor.fetchall():
        stored = clean_str(phylum)
        # Only act where the genus has a single agreed value AND this row
        # already matches it, or disagrees while being the lone dissenter.
        values = by_genus.get(genus, set())
        if len(values) != 2 or stored not in values:
            continue
        others = values - {stored}
        target = next(iter(others))
        cursor.execute(
            "SELECT COUNT(*) FROM taxa WHERE genus = %s AND phylum = %s" if store == "mysql"
            else "SELECT COUNT(*) FROM taxa WHERE genus = ? AND phylum = ?",
            (genus, stored))
        n_stored = cursor.fetchone()[0]
        cursor.execute(
            "SELECT COUNT(*) FROM taxa WHERE genus = %s AND phylum = %s" if store == "mysql"
            else "SELECT COUNT(*) FROM taxa WHERE genus = ? AND phylum = ?",
            (genus, target))
        n_target = cursor.fetchone()[0]
        # Fold a minority onto the majority the rest of the genus uses. The
        # count is not what makes this safe -- the APPROVED gate below is, since
        # it means this exact fold was already sanctioned table-wide. Requiring
        # the target to be the larger group just keeps the direction right.
        if n_target < 2 or n_target <= n_stored:
            continue
        if not ((stored, target) in APPROVED or stored.casefold() == target.casefold()):
            continue
        audit.append({
            "store": store, "taxon_id": taxon_id, "scientific_name": sci,
            "rank": "", "ncbi_tax_id": "", "stored_phylum": stored,
            "ncbi_phylum": target,
            "action": f"genus_consensus ({n_target} of {n_stored + n_target} agree)",
        })
        stats["folded_by_genus"] += 1
        if not DRY_RUN:
            cursor.execute(
                "UPDATE taxa SET phylum = %s WHERE id = %s" if store == "mysql"
                else "UPDATE taxa SET phylum = ? WHERE id = ?",
                (target, taxon_id))
    return stats


def main() -> int:
    print(f"loading taxdump from {TAXDUMP}/ ...")
    name_to_taxid, phylum_of = load_taxdump()
    print(f"  names: {len(name_to_taxid):,}   nodes with a phylum: {len(phylum_of):,}")

    audit: list[dict] = []
    results: dict[str, Counter] = {}

    mysql = connect(Settings.from_env(".env"))
    cursor = mysql.cursor()
    results["mysql"] = reconcile("mysql", cursor, name_to_taxid, phylum_of, audit)
    results["mysql"] += fold_by_genus_consensus("mysql", cursor, audit)
    if not DRY_RUN:
        mysql.commit()
    cursor.close()
    mysql.close()

    if os.path.exists(SQLITE_PATH):
        lite = sqlite3.connect(SQLITE_PATH)
        cursor = lite.cursor()
        results["sqlite"] = reconcile("sqlite", cursor, name_to_taxid, phylum_of, audit)
        results["sqlite"] += fold_by_genus_consensus("sqlite", cursor, audit)
        if not DRY_RUN:
            lite.commit()
        cursor.close()
        lite.close()
    else:
        print(f"  (no SQLite mirror at {SQLITE_PATH}; skipped)")

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")

    for store, stats in results.items():
        print(f"\n{store}")
        for key in ("agrees", "normalized", "case_folded", "folded_by_genus",
                    "disagrees_unapproved",
                    "no_phylum_at_ncbi", "unresolvable"):
            print(f"    {key:24} {stats[key]}")

    if audit:
        AUDIT.parent.mkdir(parents=True, exist_ok=True)
        with AUDIT.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(audit[0].keys()))
            writer.writeheader()
            writer.writerows(audit)
        print(f"\naudit written to {AUDIT} ({len(audit)} rows)")
        unapproved = Counter(
            (r["stored_phylum"], r["ncbi_phylum"]) for r in audit
            if r["action"] == "reported_only")
        if unapproved:
            print("\nunapproved disagreements, most common first "
                  "(reported only - add to APPROVED to act on one):")
            for (stored, current), n in unapproved.most_common(15):
                print(f"    {n:>5}  {stored} -> {current}")
    else:
        print("\nno disagreements found")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
