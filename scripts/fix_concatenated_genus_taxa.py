#!/usr/bin/env python3
"""Repair taxa rows whose stored name concatenates two genera.

69 rows read like `Lactobacillus limosilactobacillus reuteri`,
`Clavispora candida auris`, `Phocaeicola bacteroides vulgatus` -- a genus from
one classification followed by a full binomial from another. No organism bears
such a name: NCBI lookups miss, genus-level grouping files the row under a genus
it does not belong to (`Clavispora candida auris` under Clavispora while
`Candidozyma auris` groups elsewhere), and where a correctly named row also
exists one organism occupies two ML features.

The cause is ours, not MiMeDB's. `_collapse_mimedb` builds the stored name as
`genus` + `extract_species_epithet(genus, species)`, and that function only
strips the leading token of the `species` field when it equals the `genus`
column. 157 of the 2,648 MiMeDB source rows hold a `species` field whose leading
token is a *different* genus, because one of the two fields was updated for a
reclassification and the other was not, so both genera survive into the name.
`merge_duplicate_taxa.py` cleaned the subset that shared a tax ID with a
correctly named row; rows whose tax ID is unique in the table were never
grouped, and still carry the mangled name.

WHAT DECIDES THE NAME

The tax ID, as everywhere else here: it is this database's canonical key, and
`merge_duplicate_taxa.py` canonicalizes on NCBI's current scientific name for
it. A rank='species' row holding a *strain* tax ID is repointed to the species
node, which is the `retaxid` action that script's OVERRIDES performs twice
(F. prausnitzii, L. salivarius) and which `fix_sphingobacterium_paucimobilis.py`
performed for the single orphan of this class. Names are resolved offline from
the local taxdump; neither the stale genus nor the current one wins by rule.

Two guards stop that rule where it cannot be trusted, and between them they
caught every row in this set needing a hand decision:

  * EPITHET DRIFT. `same_epithet` tolerates Latin gender agreement
    (`rectale`/`rectalis`) but not a different epithet. A renamed strain is
    indistinguishable from a renamed species by name alone -- the reason tax ID
    411483 "F. prausnitzii A2-165" must not relabel F. prausnitzii as
    F. duncaniae -- so a row whose epithet disagrees with NCBI's for its tax ID
    is skipped unless OVERRIDES settles it.
  * SAME SPECIES NODE. A merge requires the surviving row's tax ID to resolve to
    the same species node as the mangled row's, so a merge can never join two
    organisms.

DECISION CLASSES

  merge                 A correctly named row already holds the target key, so
                        the mangled row cannot be renamed into it (the unique
                        key is `(genus_key, species_key)`). Dependents are
                        repointed, traits are backfilled where the survivor is
                        empty, and the mangled row is deleted. 52 rows, and all
                        52 survivors carry the evidence this corruption split
                        off -- row 1569 `Phocaeicola vulgatus` holds 10,263
                        abundances while `Phocaeicola bacteroides vulgatus`
                        holds none.
  merge_into_taxid_row  No row holds the target key, but exactly one row already
                        carries the target *tax ID* under a stale name. Same tax
                        ID is the same organism, so that row is the survivor and
                        is renamed to NCBI's current name -- the rule
                        `merge_duplicate_taxa.py` applies to a group it could
                        not see here only because the mangled row carried a
                        strain tax ID instead of the species one.
  rename                Nothing holds the key or the tax ID: NCBI moved the
                        genus and this row is the only one for the organism.
                        Renamed in place, tax ID repointed to the species node,
                        lineage repointed with it.
  skip                  A guard fired. Reported, never written.

ABUNDANCES

Where a sample holds both the mangled row and the survivor the values are
SUMMED, as `fix_candidatus_taxa.py` does: overwriting silently drops that
sample's share of its rank total, and `validate` cannot see it because the
per-rank sums were renormalized after the duplication. On the data as measured
this path does not fire -- all 69 rows carry zero abundances, zero associations
and zero adjudications -- but this database changes hourly, so the repair does
not assume what the dry run saw.

MIMEDB PROVENANCE

Nothing in the schema records MiMeDB's microbe_id (`_upsert_source_record` only
promotes `source_database`), so deleting a merged-away row drops no provenance
that was ever stored. Traits move to the survivor instead of being dropped,
which is `merge_duplicate_taxa.py`'s rule, with one documented exception below.

THE SQLITE MIRROR

The mirror is built from the curated CSV and GMrepo and holds no MiMeDB rows at
all, so none of these keys exist there. It is still checked, keyed on
`(genus_key, species_key)` -- never on `taxa.id`, which is an independent id
space in that store -- so a future mirror that does carry them is repaired
rather than silently skipped.

Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import csv
import os
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.pipeline import open_dict_rows
from gutdb.traits import dissents_from_genus, load_genus_consensus

# Trait moves refused because the survivor's genus contradicts them. Reported
# at the end; the merge still stands, one cell just does not travel with it.
REFUSED_TRAITS: list[str] = []
from gutdb.transform import (
    clean_str,
    extract_species_epithet,
    first_value,
    normalize_genus,
    normalize_phylum,
    parse_scientific_name,
    resolve_source_binomial,
    same_epithet,
)

DRY_RUN = "--apply" not in sys.argv
TAXDUMP = os.environ.get("GUTDB_TAXDUMP", "taxdump")
SQLITE_PATH = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")
MIMEDB_CSV = Path(os.environ.get("GUTDB_MIMEDB_CSV", "mimedb_microbes_v2.csv"))
AUDIT = Path("data/concatenated_genus_audit.csv")

# Moved from the mangled row onto the survivor when the survivor is empty.
# Lineage and traits only -- never keys, names, rank or tax ID. Same list as
# merge_duplicate_taxa.py, for the same reason: MiMeDB measured these traits and
# the GMrepo row that holds the evidence usually has none.
TRAIT_COLUMNS = (
    "superkingdom", "phylum", "class_name", "order_name", "family",
    "gram_stain", "oxygen_requirement", "ph_preference", "sporulation",
    "shape", "cell_arrangement", "mobility", "flagella_presence",
    "number_of_membranes", "biotic_relationship", "habitat",
    "temperature_range", "optimal_temperature", "metabolism",
    "energy_source", "human_pathogen",
)
NA_COLUMNS = ("energy_mode", "primary_food_source")

DEPENDENTS = (
    ("sample_taxon_abundances", "taxon_id"),
    ("taxon_disease_associations", "taxon_id"),
    ("taxon_disease_adjudication", "taxon_id"),
)

# NCBI rank -> taxa column, for the lineage a renamed row carries.
LINEAGE_COLUMNS = (
    ("superkingdom", ("superkingdom", "domain", "kingdom")),
    ("phylum", ("phylum",)),
    ("class_name", ("class",)),
    ("order_name", ("order",)),
    ("family", ("family",)),
)

# Rows the two guards flag, resolved by hand. Keyed on what makes the row
# corrupt -- its (genus_key, species_key) -- rather than on taxa.id, so the
# entry cannot silently land on a different row.
OVERRIDES: dict[tuple[str, str], dict] = {
    # MiMeDB record MMDBm0002385's `name` ("Propionibacterium freudenreichii
    # subsp. shermanii CIRM-BIA1") and its ncbi_tax_id (754252) both say
    # P. freudenreichii; only `species` says Porphyromonas gingivalis, a wholly
    # different organism. Its traits settle it independently: gram POSITIVE,
    # facultative, bacilli, which is P. freudenreichii -- P. gingivalis is a
    # gram-negative anaerobe. So the `species` cell is the corrupt one and the
    # epithet guard is firing on corruption, not on a reclassification. The row
    # is not P. gingivalis and must not be merged into it.
    ("propionibacterium", "porphyromonas gingivalis"): {
        "taxid": 1744,
        "note": "MiMeDB species cell holds an unrelated organism; name, tax ID "
                "and traits all say P. freudenreichii",
    },
    # tax 411466 is strain ATCC 17982, which NCBI has moved to Schaalia
    # dentiphila (species 3050224). The species Schaalia odontolytica still
    # stands as 1660 on a different type strain (ATCC 17929), and MiMeDB's
    # `name` and `species` both say odontolytica. This is the F. prausnitzii
    # shape exactly: a rank='species' row holding a strain tax ID whose strain
    # was reidentified. The row names the species, so it is given the species
    # tax ID rather than following the strain to a different organism.
    ("actinomyces", "schaalia odontolytica"): {
        "taxid": 1660,
        "no_trait_backfill": True,
        "note": "strain ATCC 17982 was reidentified as S. dentiphila; the row "
                "names the species S. odontolytica, so it takes taxid 1660. "
                "Traits are ATCC 17982's, i.e. another species', so they are "
                "not backfilled onto the survivor",
    },
    # Not drift: names.dmp records Clostridium methoxybenzovorans as a synonym
    # of 69825 Lacrimispora indolis, and strain SR3 (1286820) sits under that
    # species. NCBI itself made the epithet change.
    ("lachnoclostridium", "clostridium methoxybenzovorans"): {
        "taxid": 69825,
        "note": "C. methoxybenzovorans is a synonym of Lacrimispora indolis in "
                "names.dmp; strain SR3 sits under it",
    },
    # tax 1232666 (strain Z8) is not in the dump -- a deleted node, so the tax
    # ID cannot resolve. MiMeDB's `name` ("Staphylococcus sciuri subsp. sciuri
    # Z8") and `species` ("Mammaliicoccus sciuri") agree on the epithet, and
    # NCBI moved S. sciuri to Mammaliicoccus in 2020.
    ("staphylococcus", "mammaliicoccus sciuri"): {
        "taxid": 1296,
        "note": "taxid 1232666 is absent from the dump; both MiMeDB fields name "
                "the organism NCBI calls Mammaliicoccus sciuri (1296)",
    },
    # The merge target, row 1370, carries 1328388 -- also a deleted node -- so
    # the same-species-node guard cannot confirm the merge from the survivor's
    # side. The mangled row's own 1328337 resolves to species 573, and the two
    # names are identical, so the organism is not in doubt.
    ("hafnia", "klebsiella pneumoniae"): {
        "taxid": 573,
        "allow_unresolvable_survivor": True,
        "note": "survivor row carries deleted node 1328388; mangled row's "
                "1328337 resolves to species 573 and the names agree",
    },
}


class Taxdump:
    """Ranks, parents and names from nodes.dmp/names.dmp, read once."""

    def __init__(self, directory: str) -> None:
        names_path, nodes_path = f"{directory}/names.dmp", f"{directory}/nodes.dmp"
        if not (os.path.exists(names_path) and os.path.exists(nodes_path)):
            raise SystemExit(
                f"taxdump not found at {directory}/. Unzip taxdmp.zip there or "
                "set GUTDB_TAXDUMP.")
        self.rank_of: dict[int, str] = {}
        self.parent: dict[int, int] = {}
        with open(nodes_path, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = line.split("\t|")
                taxid = int(parts[0])
                self.parent[taxid] = int(parts[1])
                self.rank_of[taxid] = sys.intern(parts[2].strip())
        self.scientific: dict[int, str] = {}
        self.any_name: dict[str, list[int]] = {}
        with open(names_path, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = line.split("\t|")
                if len(parts) < 4:
                    continue
                taxid, name, kind = int(parts[0]), parts[1].strip(), parts[3].strip()
                if kind == "scientific name":
                    self.scientific[taxid] = name
                if kind in ("scientific name", "synonym", "equivalent name",
                            "includes", "genbank synonym"):
                    self.any_name.setdefault(name.casefold(), []).append(taxid)
        self.genus_names = {
            name for name, taxids in self.any_name.items()
            if any(self.rank_of.get(t) == "genus" for t in taxids)
            and any(self.scientific.get(t, "").casefold() == name for t in taxids)
        }

    def ancestor(self, taxid: int | None, rank: str) -> int | None:
        node, steps = taxid or 0, 0
        while node and node != 1 and steps < 60:
            if self.rank_of.get(node) == rank:
                return node
            node = self.parent.get(node, 0)
            steps += 1
        return None

    def species_node(self, taxid: int | None) -> int | None:
        """The species node at or above a tax ID, which is what a rank='species'
        row must carry; a strain or subspecies ID walks up to it."""
        if not taxid or taxid not in self.rank_of:
            return None
        return self.ancestor(taxid, "species")

    def lineage(self, taxid: int | None) -> dict[str, str]:
        out: dict[str, str] = {}
        for column, ranks in LINEAGE_COLUMNS:
            for rank in ranks:
                node = self.ancestor(taxid, rank)
                if node:
                    out[column] = self.scientific.get(node, "")
                    break
        if out.get("phylum"):
            out["phylum"] = normalize_phylum(out["phylum"])
        return {k: v for k, v in out.items() if v}

    def is_genus_name(self, token: str) -> bool:
        return clean_str(token).strip("[]").casefold() in self.genus_names

    def has_name(self, name: str) -> bool:
        return clean_str(name).casefold() in self.any_name


def loader_would_rekey(path: Path) -> dict[tuple[str, str], dict]:
    """Old key -> what the fixed loader derives instead, over every source row.

    This is the authoritative scope: a row whose key the fixed loader no longer
    produces, left unrepaired, gets a fresh duplicate beside it on the next
    ingest, because the loader recomputes the new key and inserts rather than
    updates.
    """
    rekeyed: dict[tuple[str, str], dict] = {}
    rows = 0
    for row in open_dict_rows(path):
        rows += 1
        genus = normalize_genus(first_value(row, "genus"))
        field = first_value(row, "species", "name")
        old = (genus.casefold(), extract_species_epithet(genus, field).casefold())
        new_genus, new_species = resolve_source_binomial(genus, field)
        new = (new_genus.casefold(), new_species.casefold())
        if old == new:
            continue
        entry = rekeyed.setdefault(old, {"new": new, "taxids": set(), "microbe_ids": []})
        taxid = first_value(row, "ncbi_tax_id")
        if taxid.isdigit():
            entry["taxids"].add(int(taxid))
        entry["microbe_ids"].append(first_value(row, "microbe_id"))
    print(f"  source rows read: {rows}")
    print(f"  keys the fixed loader would no longer produce: {len(rekeyed)}")
    return rekeyed


def scan_for_two_genus_names(rows: list[dict], dump: Taxdump) -> set[int]:
    """Rows whose species field begins with a genus other than their own.

    The independent check: a row is flagged when the first token of `species` is
    a genus scientific name in the dump and differs from the `genus` column.
    `Candidatus` is skipped (it is a status prefix, not a genus), and so is a
    legitimate infraspecific name -- `Kocuria marina subsp. indica` would
    otherwise flag because *marina* is coincidentally a genus, and it is
    recognised by `genus + first token` ("Kocuria marina") being a name NCBI
    knows, which makes the first token this row's own epithet.
    """
    found: set[int] = set()
    for row in rows:
        species, genus = row["species"], row["genus"]
        if " " not in species:
            continue
        first = species.split()[0]
        if species.casefold().startswith("candidatus ") or \
                genus.casefold().startswith("candidatus"):
            continue
        if not dump.is_genus_name(first) or first.casefold() == genus.casefold():
            continue
        if dump.has_name(f"{genus} {first}"):
            continue
        found.add(row["id"])
    return found


def load_taxa(cursor, placeholder: str) -> list[dict]:
    cursor.execute(
        "SELECT id, scientific_name, genus, species, genus_key, species_key, "
        "ncbi_tax_id, taxonomic_rank, source_database, superkingdom, phylum, "
        "class_name, order_name, family FROM taxa")
    columns = [d[0] for d in cursor.description]
    return [dict(zip(columns, row)) for row in cursor.fetchall()]


def dependents(cursor, placeholder: str, tables: set[str], taxon_id: int) -> dict[str, int]:
    counts = {}
    for table, column in DEPENDENTS:
        if table not in tables:
            counts[table] = 0
            continue
        cursor.execute(
            f"SELECT COUNT(*) FROM `{table}` WHERE `{column}`={placeholder}", (taxon_id,))
        counts[table] = cursor.fetchone()[0]
    return counts


def resolve_target(row: dict, dump: Taxdump) -> tuple[dict | None, str]:
    """NCBI's name for this row's organism, or None with the reason it is skipped."""
    key = (row["genus_key"], row["species_key"])
    override = OVERRIDES.get(key)
    taxid = override["taxid"] if override else row["ncbi_tax_id"]
    note = override["note"] if override else ""

    if not taxid:
        return None, "row carries no ncbi_tax_id, so nothing decides the name"
    node = taxid if dump.rank_of.get(taxid) == "species" else dump.species_node(taxid)
    if not node:
        return None, f"tax ID {taxid} has no species node in the dump"
    name = dump.scientific.get(node, "")
    if not name:
        return None, f"species node {node} has no scientific name in the dump"
    genus, species = parse_scientific_name(name)
    if not genus or not species:
        return None, f"NCBI name {name!r} does not parse as a binomial"

    stored_epithet = row["species"].split()[-1] if row["species"] else ""
    if not override and not same_epithet(stored_epithet, species):
        return None, (f"epithet drift: row says {stored_epithet!r}, NCBI says "
                      f"{species!r} for tax ID {taxid} -- needs an OVERRIDES entry")
    return {
        "name": name, "genus": genus, "species": species,
        "genus_key": genus.casefold(), "species_key": species.casefold(),
        "taxid": node, "note": note,
        "no_trait_backfill": bool(override and override.get("no_trait_backfill")),
        "allow_unresolvable_survivor": bool(
            override and override.get("allow_unresolvable_survivor")),
    }, ""


def lineage_for(target: dict, row: dict, dump: Taxdump, known: dict[str, set[str]]) -> dict:
    """NCBI's lineage for the target, in the spellings this table already uses.

    A renamed row keeps a lineage it no longer belongs to otherwise. Where the
    dump's name for a level is already in the table it is used as the table
    spells it; where it is absent but the row's current value names the same
    NCBI node, the current value is kept, so renaming one row does not introduce
    an 87th spelling of a class the rest of the table agrees on.
    """
    proposed = dump.lineage(target["taxid"])
    out: dict[str, str] = {}
    for column, value in proposed.items():
        current = clean_str(row.get(column))
        existing = known.get(column, set())
        match = next((e for e in existing if e.casefold() == value.casefold()), None)
        if match:
            value = match
        elif current and current.casefold() in existing:
            current_nodes = set(dump.any_name.get(current.casefold(), []))
            value_nodes = set(dump.any_name.get(value.casefold(), []))
            if current_nodes & value_nodes:
                value = current
        if value != current:
            out[column] = value
    return out


def novel_lineage(lineage: dict, known: dict[str, set[str]]) -> list[str]:
    """Lineage values no other row in the table uses yet.

    Not an error: NCBI's placement for a row whose name comes from NCBI. But a
    value of one row where the siblings agree on an older spelling is the shape
    of a lineage splitting in two, so it is surfaced rather than written
    quietly. The table already carries such minorities (3 rows Serinales beside
    71 Saccharomycetales), so this reports and proceeds.
    """
    return [f"{column}={value}" for column, value in lineage.items()
            if value.casefold() not in {k.casefold() for k in known.get(column, set())}]


def plan(rows: list[dict], scope: set[int], dump: Taxdump, cursor, placeholder: str,
         tables: set[str]) -> tuple[list[dict], list[dict]]:
    by_key = {(r["genus_key"], r["species_key"]): r for r in rows}
    by_id = {r["id"]: r for r in rows}
    known = {column: {clean_str(r[column]) for r in rows if clean_str(r[column])}
             for column, _ in LINEAGE_COLUMNS}

    targets: dict[int, dict] = {}
    skips: list[dict] = []
    for taxon_id in sorted(scope):
        row = by_id[taxon_id]
        target, reason = resolve_target(row, dump)
        if not target:
            skips.append({"row": row, "reason": reason})
            continue
        targets[taxon_id] = target

    groups: dict[tuple[str, str], list[int]] = {}
    for taxon_id, target in targets.items():
        groups.setdefault((target["genus_key"], target["species_key"]), []).append(taxon_id)

    actions: list[dict] = []
    for key, members in sorted(groups.items()):
        target = targets[members[0]]
        survivor_row = by_key.get(key)
        rename_survivor = False

        if survivor_row is None:
            # Nothing holds the name. Does anything already hold the tax ID?
            holders = [r for r in rows
                       if r["ncbi_tax_id"] == target["taxid"] and r["id"] not in members]
            if len(holders) > 1:
                for taxon_id in members:
                    skips.append({"row": by_id[taxon_id], "reason":
                                  f"tax ID {target['taxid']} is held by "
                                  f"{len(holders)} rows; ambiguous survivor"})
                continue
            if holders:
                survivor_row = holders[0]
                if not same_epithet(survivor_row["species"].split()[-1]
                                    if survivor_row["species"] else "",
                                    target["species"]):
                    for taxon_id in members:
                        skips.append({"row": by_id[taxon_id], "reason":
                                      f"row {survivor_row['id']} holds tax ID "
                                      f"{target['taxid']} but its epithet differs"})
                    continue
                rename_survivor = True
            else:
                # Rename one member in place; prefer the row that already holds
                # the species tax ID, so the row needing no retaxid survives.
                members.sort(key=lambda i: (by_id[i]["ncbi_tax_id"] != target["taxid"], i))
                survivor_row = by_id[members[0]]
                members = members[1:]
                lineage = lineage_for(target, survivor_row, dump, known)
                actions.append({
                    "action": "rename", "row": survivor_row, "target": target,
                    "lineage": lineage, "novel": novel_lineage(lineage, known),
                    "deps": dependents(cursor, placeholder, tables, survivor_row["id"]),
                })

        if rename_survivor:
            lineage = lineage_for(target, survivor_row, dump, known)
            actions.append({
                "action": "rename_survivor", "row": survivor_row, "target": target,
                "lineage": lineage, "novel": novel_lineage(lineage, known),
                "deps": dependents(cursor, placeholder, tables, survivor_row["id"]),
            })

        for taxon_id in members:
            if taxon_id == survivor_row["id"]:
                continue
            row = by_id[taxon_id]
            survivor_node = dump.species_node(survivor_row["ncbi_tax_id"])
            if survivor_node != target["taxid"] and not target["allow_unresolvable_survivor"]:
                skips.append({"row": row, "reason":
                              f"survivor row {survivor_row['id']} "
                              f"(tax {survivor_row['ncbi_tax_id']}) resolves to species "
                              f"node {survivor_node}, not {target['taxid']}"})
                continue
            actions.append({
                # the row's own target, not the group's: the key is shared but
                # the reasoning that got each row here is not.
                "action": "merge", "row": row, "target": targets[taxon_id],
                "survivor": survivor_row,
                "deps": dependents(cursor, placeholder, tables, taxon_id),
            })
    return actions, skips


def backfill_traits(cursor, placeholder: str, survivor_id: int, loser_id: int) -> list[str]:
    """Move the loser's traits onto empty survivor fields, minus implausible ones.

    A value the survivor's genus contradicts is refused. This is the path that
    carried MiMeDB's gram-positive call for Bacteroides vulgatus onto
    Phocaeicola vulgatus and its 10,263 abundance rows: the rule is right, but
    it moves a trait off a row nobody queries onto the row everybody does, so a
    wrong cell stops being harmless exactly when it starts to matter. An empty
    field is honest; a wrong one is not.
    """
    consensus = load_genus_consensus(cursor, placeholder)
    cursor.execute(f"SELECT genus_key FROM taxa WHERE id={placeholder}", (survivor_id,))
    survivor_genus = (cursor.fetchone() or ("",))[0]

    cursor.execute(
        f"SELECT {', '.join(TRAIT_COLUMNS + NA_COLUMNS)} FROM taxa "
        f"WHERE id IN ({placeholder}, {placeholder}) ORDER BY id={placeholder} DESC",
        (survivor_id, loser_id, survivor_id))
    survivor, loser = cursor.fetchall()
    moved = []
    for index, column in enumerate(TRAIT_COLUMNS + NA_COLUMNS):
        empty = survivor[index] is None or (
            column in NA_COLUMNS and str(survivor[index]) in ("", "N/A"))
        if empty and loser[index] is not None and str(loser[index]) not in ("", "N/A"):
            agreed = dissents_from_genus(
                consensus, survivor_genus, column, loser[index])
            if agreed is not None:
                REFUSED_TRAITS.append(
                    f"taxon {survivor_id}: refused {column}="
                    f"{str(loser[index]).strip()!r} -- genus {survivor_genus!r} "
                    f"agrees on {agreed!r}")
                continue
            cursor.execute(
                f"UPDATE taxa SET {column}={placeholder} WHERE id={placeholder}",
                (loser[index], survivor_id))
            moved.append(column)
    return moved


def repoint(cursor, placeholder: str, tables: set[str], loser_id: int,
            survivor_id: int) -> dict[str, int]:
    """Move dependents onto the survivor, summing abundances where both are held.

    Abundances are keyed (sample, taxon): where a sample holds both rows the two
    values are the same organism measured once under two names, so they are
    summed. Overwriting would drop that sample's share of its rank total, and
    the per-rank sums were renormalized after the duplication, so `validate`
    would not see the loss. Done row by row because MySQL cannot reference the
    target table inside ON DUPLICATE KEY UPDATE when it is also the source.
    """
    moved = {"abundances_summed": 0, "abundances_moved": 0,
             "associations_moved": 0, "adjudications_moved": 0}
    if "sample_taxon_abundances" in tables:
        cursor.execute(
            f"SELECT sample_id, relative_abundance FROM sample_taxon_abundances "
            f"WHERE taxon_id={placeholder}", (loser_id,))
        for sample_id, value in cursor.fetchall():
            cursor.execute(
                f"SELECT relative_abundance FROM sample_taxon_abundances "
                f"WHERE sample_id={placeholder} AND taxon_id={placeholder}",
                (sample_id, survivor_id))
            held = cursor.fetchone()
            if held:
                cursor.execute(
                    f"UPDATE sample_taxon_abundances SET relative_abundance={placeholder} "
                    f"WHERE sample_id={placeholder} AND taxon_id={placeholder}",
                    (float(held[0]) + float(value), sample_id, survivor_id))
                moved["abundances_summed"] += 1
            else:
                cursor.execute(
                    f"INSERT INTO sample_taxon_abundances (sample_id, taxon_id, "
                    f"relative_abundance) VALUES ({placeholder},{placeholder},{placeholder})",
                    (sample_id, survivor_id, value))
                moved["abundances_moved"] += 1
        cursor.execute(
            f"DELETE FROM sample_taxon_abundances WHERE taxon_id={placeholder}", (loser_id,))

    ignore = "IGNORE " if placeholder == "%s" else "OR IGNORE "
    for table, counter in (("taxon_disease_associations", "associations_moved"),
                           ("taxon_disease_adjudication", "adjudications_moved")):
        if table not in tables:
            continue
        cursor.execute(
            f"UPDATE {ignore}{table} SET taxon_id={placeholder} WHERE taxon_id={placeholder}",
            (survivor_id, loser_id))
        moved[counter] += cursor.rowcount
        cursor.execute(f"DELETE FROM {table} WHERE taxon_id={placeholder}", (loser_id,))
    return moved


def apply_action(cursor, placeholder: str, tables: set[str], action: dict) -> dict:
    row, target = action["row"], action["target"]
    result = {"traits": [], "moved": {}}
    if action["action"] in ("rename", "rename_survivor"):
        updates = {
            "scientific_name": target["name"], "genus": target["genus"],
            "species": target["species"], "genus_key": target["genus_key"],
            "species_key": target["species_key"], "ncbi_tax_id": target["taxid"],
        }
        updates.update(action["lineage"])
        assignments = ", ".join(f"{c}={placeholder}" for c in updates)
        cursor.execute(f"UPDATE taxa SET {assignments} WHERE id={placeholder}",
                       tuple(updates.values()) + (row["id"],))
    else:
        survivor = action["survivor"]
        if not target["no_trait_backfill"]:
            result["traits"] = backfill_traits(cursor, placeholder, survivor["id"], row["id"])
        result["moved"] = repoint(cursor, placeholder, tables, row["id"], survivor["id"])
        cursor.execute(f"DELETE FROM taxa WHERE id={placeholder}", (row["id"],))
    return result


def run_store(store: str, cursor, tables: set[str], rekeyed: dict,
              dump: Taxdump, audit: list[dict]) -> dict[str, int]:
    placeholder = "%s" if store == "mysql" else "?"
    rows = load_taxa(cursor, placeholder)
    print(f"  taxa rows: {len(rows)}")

    by_rekeyed = {r["id"] for r in rows
                  if (r["genus_key"], r["species_key"]) in rekeyed}
    by_scan = scan_for_two_genus_names(rows, dump)
    print(f"  rows the fixed loader would re-key: {len(by_rekeyed)}")
    print(f"  rows the two-genus scan finds:      {len(by_scan)}")
    only_loader = by_rekeyed - by_scan
    only_scan = by_scan - by_rekeyed
    by_id = {r["id"]: r for r in rows}
    for taxon_id in sorted(only_loader):
        print(f"    loader only: {taxon_id} {by_id[taxon_id]['scientific_name']!r}")
    for taxon_id in sorted(only_scan):
        print(f"    scan only:   {taxon_id} {by_id[taxon_id]['scientific_name']!r}")
    scope = by_rekeyed | by_scan
    print(f"  rows in scope: {len(scope)}")
    if not scope:
        return {}

    actions, skips = plan(rows, scope, dump, cursor, placeholder, tables)
    stats = {"renamed": 0, "renamed_survivor": 0, "merged": 0, "skipped": len(skips),
             "traits_backfilled": 0, "abundances_summed": 0, "abundances_moved": 0,
             "associations_moved": 0, "adjudications_moved": 0}

    for action in sorted(actions, key=lambda a: (a["action"] != "rename_survivor",
                                                 a["row"]["id"])):
        row, target = action["row"], action["target"]
        deps = action["deps"]
        found = ("loader+scan" if row["id"] in by_rekeyed and row["id"] in by_scan
                 else "loader" if row["id"] in by_rekeyed else "scan")
        line = f"    {action['action']:<16} {row['id']:>8} {row['scientific_name'][:46]:<46}"
        if action["action"] == "merge":
            print(f"{line} -> {action['survivor']['id']} {target['name']}"
                  f"   deps={sum(deps.values())}")
        else:
            print(f"{line} -> {target['name']}  tax {row['ncbi_tax_id']} -> "
                  f"{target['taxid']}  deps={sum(deps.values())}  "
                  f"lineage={action['lineage'] or '-'}"
                  + (f"  NEW TO TABLE: {', '.join(action['novel'])}"
                     if action.get("novel") else ""))
        result = {"traits": [], "moved": {}}
        if not DRY_RUN:
            result = apply_action(cursor, placeholder, tables, action)
        stats["renamed"] += action["action"] == "rename"
        stats["renamed_survivor"] += action["action"] == "rename_survivor"
        stats["merged"] += action["action"] == "merge"
        stats["traits_backfilled"] += len(result["traits"])
        for key, value in result["moved"].items():
            stats[key] += value
        audit.append({
            "store": store, "action": action["action"], "taxon_id": row["id"],
            "scientific_name": row["scientific_name"],
            "ncbi_tax_id": row["ncbi_tax_id"],
            "target_name": target["name"], "target_taxid": target["taxid"],
            "survivor_id": action.get("survivor", {}).get("id", ""),
            "survivor_name": action.get("survivor", {}).get("scientific_name", ""),
            "abundances": deps["sample_taxon_abundances"],
            "associations": deps["taxon_disease_associations"],
            "adjudications": deps["taxon_disease_adjudication"],
            "lineage_changes": "; ".join(f"{k}={v}" for k, v in
                                         action.get("lineage", {}).items()),
            "traits_backfilled": "; ".join(result["traits"]),
            "lineage_new_to_table": "; ".join(action.get("novel", [])),
            "found_by": found, "note": target["note"],
        })

    for skip in skips:
        row = skip["row"]
        print(f"    {'SKIPPED':<16} {row['id']:>8} {row['scientific_name'][:46]:<46}"
              f" -- {skip['reason']}")
        audit.append({
            "store": store, "action": "skipped", "taxon_id": row["id"],
            "scientific_name": row["scientific_name"],
            "ncbi_tax_id": row["ncbi_tax_id"], "target_name": "", "target_taxid": "",
            "survivor_id": "", "survivor_name": "", "abundances": "",
            "associations": "", "adjudications": "", "lineage_changes": "",
            "traits_backfilled": "", "lineage_new_to_table": "",
            "found_by": "", "note": skip["reason"],
        })
    return stats


def main() -> int:
    if not MIMEDB_CSV.exists():
        raise SystemExit(f"{MIMEDB_CSV} not found; set GUTDB_MIMEDB_CSV or link it here.")
    print(f"reading {MIMEDB_CSV} ...")
    rekeyed = loader_would_rekey(MIMEDB_CSV)
    print(f"loading taxdump from {TAXDUMP}/ ...")
    dump = Taxdump(TAXDUMP)
    print(f"  nodes {len(dump.rank_of):,}   genus scientific names {len(dump.genus_names):,}")

    audit: list[dict] = []

    print("\nmysql")
    conn = connect(Settings.from_env(".env"))
    cursor = conn.cursor()
    cursor.execute("SHOW TABLES")
    tables = {r[0] for r in cursor.fetchall()}
    mysql_stats = run_store("mysql", cursor, tables, rekeyed, dump, audit)
    if not DRY_RUN:
        conn.commit()
    cursor.close()
    conn.close()

    print("\nsqlite")
    sqlite_stats = {}
    if Path(SQLITE_PATH).exists():
        lite = sqlite3.connect(SQLITE_PATH)
        cursor = lite.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = {r[0] for r in cursor.fetchall()}
        sqlite_stats = run_store("sqlite", cursor, tables, rekeyed, dump, audit)
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
    else:
        print("\nno concatenated-genus rows found - nothing to repair")

    for store, stats in (("mysql", mysql_stats), ("sqlite", sqlite_stats)):
        if stats:
            print(f"\n{store}")
            for key, value in stats.items():
                print(f"  {key:26} {value}")

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
