#!/usr/bin/env python3
"""Repair cross-kingdom homonym lineages in gut_microbiome_new_full.csv.

The curated table's Phylum/Class/Order/Family were populated by resolving each
genus name without constraining the search to a domain, so every genus that is
a homonym shared with an insect or a plant picked up the wrong lineage:

    Bacillus     -> Arthropoda / Insecta / Phasmatodea / Bacillidae  (stick insect)
    Serratia     -> Arthropoda / Insecta / Coleoptera / Lampyridae   (firefly)
    Yersinia     -> Arthropoda / Insecta / Mantodea / Amelidae       (mantis)
    Rhodococcus  -> Arthropoda / Insecta / Hemiptera / Coccidae      (scale insect)

This is the same failure mode TaxonomyClient.resolve_name's domain guard was
later written to prevent, but the guard cannot help a value already sitting in
the file -- and because this CSV is the mirror's base taxon source and MySQL's
`load-csv` input, the bad lineage is reproduced on every rebuild. Fixing it
here is what makes the downstream fixes durable.

The test is not "this phylum looks non-microbial": the curated table
deliberately holds real helminths and plants from the MiMeDB reference set, and
Nematoda or Streptophyta is correct for those. Instead each row's organism is
resolved in the local NCBI dump and its stored lineage compared with NCBI's.
A row is rewritten only where the dump resolves the organism unambiguously and
disagrees; anything unresolvable is reported and left alone.

Writes the repaired file in place and a before/after audit. Dry run by default;
pass --apply to write.
"""
from __future__ import annotations

import csv
import os
import shutil
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.transform import clean_str, normalize_phylum

DRY_RUN = "--apply" not in sys.argv
TAXDUMP = os.environ.get("GUTDB_TAXDUMP", "taxdump")
CSV_PATH = Path("gut_microbiome_new_full.csv")
AUDIT = Path("data/curated_csv_lineage_audit.csv")

# CSV column -> NCBI rank
COLUMN_RANK = {"Phylum": "phylum", "Class": "class", "Order": "order", "Family": "family"}
SUPERKINGDOM_RANKS = ("superkingdom", "domain", "kingdom")

# Within-domain reclassifications this script may apply. Everything not listed
# is reported only: a dated name is a taxonomic decision, not a repair, and the
# default is to surface it rather than normalize it away. Each entry mirrors one
# in scripts/normalize_reclassified_phyla.py so the CSV and the two databases
# agree about which moves have been sanctioned.
APPROVED_RECLASSIFICATIONS = {
    ("Proteobacteria", "Campylobacterota"),
}


def load_taxdump() -> tuple[dict[str, int], dict[int, dict[str, str]]]:
    """Return (name -> taxid, taxid -> {rank: name}) for the four ranks above."""
    names_path, nodes_path = f"{TAXDUMP}/names.dmp", f"{TAXDUMP}/nodes.dmp"
    if not (os.path.exists(names_path) and os.path.exists(nodes_path)):
        raise SystemExit(
            f"taxdump not found at {TAXDUMP}/. Unzip taxdmp.zip there or set GUTDB_TAXDUMP.")

    # Every candidate, not the first: "Bacillus" is both a bacterial genus and a
    # stick insect, and keeping whichever line names.dmp happened to list first
    # is precisely the homonym bug this script exists to repair.
    name_to_taxid: dict[str, list[int]] = {}
    synonyms: dict[str, list[int]] = {}
    with open(names_path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split("\t|")
            if len(parts) < 4:
                continue
            name = parts[1].strip().casefold()
            klass = parts[3].strip()
            if klass == "scientific name":
                name_to_taxid.setdefault(name, []).append(int(parts[0]))
            elif klass in ("synonym", "equivalent name", "includes"):
                synonyms.setdefault(name, []).append(int(parts[0]))
    for name, taxids in synonyms.items():
        name_to_taxid.setdefault(name, []).extend(taxids)

    parent: dict[int, int] = {}
    rank_of: dict[int, str] = {}
    wanted_ranks = set(COLUMN_RANK.values()) | set(SUPERKINGDOM_RANKS)
    with open(nodes_path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split("\t|")
            taxid, par, rank = int(parts[0]), int(parts[1]), parts[2].strip()
            parent[taxid] = par
            if rank in wanted_ranks:
                rank_of[taxid] = rank

    node_name: dict[int, str] = {}
    with open(names_path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split("\t|")
            if len(parts) < 4 or parts[3].strip() != "scientific name":
                continue
            taxid = int(parts[0])
            if taxid in rank_of:
                node_name[taxid] = parts[1].strip()

    class Lineages:
        """taxid -> {rank: name}, resolved on demand.

        Memoised root-ward: each node stores the complete set of ranks at or
        above it, built by walking down from the root-most unmemoised ancestor
        so a deeper node of the same rank overwrites a shallower one. Resolving
        all 3M nodes up front in dict order was wrong -- a node could be stored
        with a partial lineage and then be inherited by its descendants, which
        is how species rows ended up with an order and family but no phylum.
        """

        def __init__(self) -> None:
            self._memo: dict[int, dict[str, str]] = {}

        def get(self, taxid: int, default=None):
            if not taxid:
                return default if default is not None else {}
            if taxid in self._memo:
                return self._memo[taxid]
            chain: list[int] = []
            node = taxid
            while node and node != 1 and node not in self._memo:
                chain.append(node)
                node = parent.get(node, 0)
            accumulated = dict(self._memo.get(node, {}))
            for link in reversed(chain):
                if link in rank_of:
                    accumulated[rank_of[link]] = node_name.get(link, "")
                self._memo[link] = dict(accumulated)
            return self._memo.get(taxid, accumulated)

    return name_to_taxid, Lineages()


def superkingdom_of(lineage: dict[str, str]) -> str:
    """The top-level domain in a resolved lineage, whichever rank NCBI used."""
    for rank in ("superkingdom", "domain"):
        if lineage.get(rank):
            return lineage[rank]
    return lineage.get("kingdom", "")


def domain_signature(lineage: dict[str, str]) -> tuple[str, str]:
    """(superkingdom, kingdom), the pair that decides "different organism".

    Superkingdom alone is too coarse inside Eukaryota: an insect and a water
    lily are both Eukaryota, so a moth lineage sitting on Victoria amazonica
    looks like agreement. Kingdom separates Metazoa from Viridiplantae.
    """
    top = ""
    for rank in ("superkingdom", "domain"):
        if lineage.get(rank):
            top = lineage[rank]
            break
    return top.casefold(), lineage.get("kingdom", "").casefold()


PROKARYOTE_DOMAINS = {"bacteria", "archaea"}


def resolve_row(row: dict, name_to_taxid, lineage_of) -> tuple[int | None, str, str]:
    """Resolve a row's organism, refusing to guess between homonyms.

    Returns (taxid, matched_on, note). A binomial is specific enough to be
    unambiguous in practice and is tried first. A genus-only row is resolved
    against prokaryotes only: this is a gut-microbiome taxon table, so a bare
    genus here is a bacterium or archaeon, while the plants and helminths the
    MiMeDB reference set deliberately holds all carry a species epithet and so
    never reach the genus fallback. If the restriction still leaves more than
    one candidate, nothing is written -- a tie is reported, not broken.
    """
    genus = clean_str(row.get("Genus"))
    species = clean_str(row.get("Species"))

    if genus and species:
        candidates = name_to_taxid.get(f"{genus} {species}".casefold(), [])
        if len(candidates) == 1:
            return candidates[0], f"{genus} {species}".casefold(), "binomial"
        if len(candidates) > 1:
            prokaryotes = [
                t for t in candidates
                if superkingdom_of(lineage_of.get(t, {})).casefold() in PROKARYOTE_DOMAINS]
            if len(prokaryotes) == 1:
                return prokaryotes[0], f"{genus} {species}".casefold(), "binomial/prokaryote"
            return None, "", f"ambiguous binomial ({len(candidates)} candidates)"

    if genus:
        candidates = name_to_taxid.get(genus.casefold(), [])
        prokaryotes = [
            t for t in candidates
            if superkingdom_of(lineage_of.get(t, {})).casefold() in PROKARYOTE_DOMAINS]
        if len(prokaryotes) == 1:
            return prokaryotes[0], genus.casefold(), "genus/prokaryote"
        if len(prokaryotes) > 1:
            return None, "", f"ambiguous genus ({len(prokaryotes)} prokaryote candidates)"
        if candidates:
            return None, "", "genus resolves to no prokaryote"
    return None, "", "not found"


def main() -> int:
    if not CSV_PATH.exists():
        raise SystemExit(f"{CSV_PATH} not found; run from the repo root.")

    print(f"loading taxdump from {TAXDUMP}/ ...")
    name_to_taxid, lineage_of = load_taxdump()
    print(f"  distinct names: {len(name_to_taxid):,}")

    with CSV_PATH.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []
        rows = list(reader)
    print(f"  rows: {len(rows)}")

    audit: list[dict] = []
    stats = Counter()
    unresolved_notes = Counter()

    for row in rows:
        taxid, matched_on, note = resolve_row(row, name_to_taxid, lineage_of)
        if not taxid:
            stats["unresolvable"] += 1
            unresolved_notes[note] += 1
            continue

        lineage = lineage_of.get(taxid, {})
        if not lineage:
            stats["no_lineage_at_ncbi"] += 1
            continue

        changes = {}
        for column, rank in COLUMN_RANK.items():
            current = clean_str(row.get(column))
            correct = lineage.get(rank, "")
            if column == "Phylum":
                correct = normalize_phylum(correct)
                current = normalize_phylum(current)
            if correct and current and current != correct:
                changes[column] = (current, correct)

        if not changes:
            stats["agrees"] += 1
            continue

        # Only a CROSS-KINGDOM mismatch is rewritten. If the stored phylum and
        # the organism sit in different superkingdoms, the stored lineage belongs
        # to a different organism entirely -- a homonym, unambiguously wrong. If
        # they share a superkingdom the disagreement is a reclassification within
        # the right domain (Euryarchaeota -> Methanobacteriota), which is a
        # taxonomic judgement, not a repair, so it is reported rather than made.
        if "Phylum" in changes:
            stored_phylum = changes["Phylum"][0]
            stored_candidates = name_to_taxid.get(stored_phylum.casefold(), [])
            stored_sig = next(
                (domain_signature(lineage_of.get(t, {})) for t in stored_candidates
                 if any(domain_signature(lineage_of.get(t, {})))), ("", ""))
            correct_sig = domain_signature(lineage)
            correct_phylum = changes["Phylum"][1]

            # Three ways the stored value cannot be right, none of them a
            # taxonomic judgement:
            #   - it names an organism in a different kingdom (a homonym)
            #   - it is the same name in different capitalization
            #   - it is not a phylum NCBI recognizes at all (a misspelling)
            different_organism = (
                any(stored_sig) and any(correct_sig) and stored_sig != correct_sig)
            case_variant = stored_phylum.casefold() == correct_phylum.casefold()
            unknown_name = not stored_candidates
            approved = (stored_phylum, correct_phylum) in APPROVED_RECLASSIFICATIONS
            if not (different_organism or case_variant or unknown_name or approved):
                stats["reclassification_not_homonym"] += 1
                audit.append({
                    "id": row.get("id", ""),
                    "organism": f"{row.get('Genus','')} {row.get('Species','')}".strip(),
                    "matched_on": matched_on, "ncbi_tax_id": taxid,
                    "action": "reported_only_reclassification",
                    "column": "Phylum",
                    "before": f"{stored_phylum} ({stored_sig[0] or '?'})",
                    "after": f"{correct_phylum} ({correct_sig[0] or '?'})",
                })
                continue

        if "Phylum" not in changes:
            stats["sub_phylum_disagreement_only"] += 1
            audit.append({
                "id": row.get("id", ""), "organism": f"{row.get('Genus','')} {row.get('Species','')}".strip(),
                "matched_on": matched_on, "ncbi_tax_id": taxid, "action": "reported_only",
                "column": " ".join(changes), "before": " | ".join(v[0] for v in changes.values()),
                "after": " | ".join(v[1] for v in changes.values()),
            })
            continue

        stats["lineage_rewritten"] += 1
        audit.append({
            "id": row.get("id", ""), "organism": f"{row.get('Genus','')} {row.get('Species','')}".strip(),
            "matched_on": matched_on, "ncbi_tax_id": taxid, "action": "rewritten",
            "column": " ".join(changes),
            "before": " | ".join(f"{c}={v[0]}" for c, v in changes.items()),
            "after": " | ".join(f"{c}={v[1]}" for c, v in changes.items()),
        })
        if not DRY_RUN:
            for column, (_, correct) in changes.items():
                row[column] = correct

    # A row whose organism the dump cannot name keeps its wrong lineage, so an
    # insect order can survive on a bacterium ("Rhodococcus fascians" is in the
    # dump only as "Bacterium fascians"). Where the rest of the genus agrees on
    # a rank, the row adopts it; where the genus disagrees, the value is cleared
    # rather than left describing a different organism.
    ANIMAL_PHYLA = {"arthropoda", "chordata", "mollusca", "cnidaria", "annelida",
                    "echinodermata", "porifera"}
    by_genus: dict[str, list[dict]] = {}
    for row in rows:
        by_genus.setdefault(clean_str(row.get("Genus")), []).append(row)
    for row in rows:
        if clean_str(row.get("Phylum")).casefold() not in ANIMAL_PHYLA:
            continue
        siblings = [r for r in by_genus.get(clean_str(row.get("Genus")), []) if r is not row]
        if not siblings:
            continue
        agreed, cleared = {}, []
        for column in COLUMN_RANK:
            values = {clean_str(r.get(column)) for r in siblings if clean_str(r.get(column))}
            values -= {v for v in values if v.casefold() in ANIMAL_PHYLA}
            if len(values) == 1:
                agreed[column] = next(iter(values))
            else:
                cleared.append(column)
        if "Phylum" not in agreed:
            continue
        stats["genus_consensus"] += 1
        audit.append({
            "id": row.get("id", ""),
            "organism": f"{row.get('Genus','')} {row.get('Species','')}".strip(),
            "matched_on": "genus consensus", "ncbi_tax_id": "", "action": "genus_consensus",
            "column": " ".join(list(agreed) + [f"{c}(cleared)" for c in cleared]),
            "before": " | ".join(f"{c}={clean_str(row.get(c))}" for c in COLUMN_RANK),
            "after": " | ".join(
                [f"{c}={v}" for c, v in agreed.items()] + [f"{c}=" for c in cleared]),
        })
        if not DRY_RUN:
            for column, value in agreed.items():
                row[column] = value
            for column in cleared:
                row[column] = ""

    if not DRY_RUN:
        backup = CSV_PATH.with_suffix(".csv.before_lineage_fix")
        shutil.copy2(CSV_PATH, backup)
        with CSV_PATH.open("w", newline="", encoding="utf-8") as fh:
            # csv.writer defaults to CRLF; this file is LF, and rewriting 2,116
            # line endings would bury 125 changed cells in a whole-file diff.
            writer = csv.DictWriter(fh, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nrewrote {CSV_PATH} (previous version kept at {backup.name})")
    else:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")

    for key in ("agrees", "lineage_rewritten", "genus_consensus",
                "reclassification_not_homonym",
                "sub_phylum_disagreement_only", "no_lineage_at_ncbi", "unresolvable"):
        print(f"    {key:30} {stats[key]}")
    if unresolved_notes:
        print("\n  why rows went unresolved:")
        for note, n in unresolved_notes.most_common():
            print(f"    {n:>5}  {note}")

    if audit:
        AUDIT.parent.mkdir(parents=True, exist_ok=True)
        with AUDIT.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(audit[0].keys()))
            writer.writeheader()
            writer.writerows(audit)
        print(f"\naudit written to {AUDIT} ({len(audit)} rows)")
        rewritten = [a for a in audit if a["action"] == "rewritten"]
        if rewritten:
            print(f"\nphylum rewrites ({len(rewritten)}):")
            for entry in rewritten[:40]:
                print(f"    {entry['organism']:34} {entry['before'][:60]}")
                print(f"    {'':34} -> {entry['after'][:60]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
