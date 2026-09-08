"""Clean up and back-fill gut_microbiome.sqlite.

Run after build_gutdb_sqlite.py. Three kinds of work, in this order:

1. Deduplication — merge disease rows that are the same MeSH concept under
   different spellings, after back-filling their MeSH ids from the source CSVs.
2. Back-fill — fill NULL/empty cells from sources the build did not read:
   gutdb/final_taxa.tsv, mimedb_microbes_v2.csv,
   gut_microbiome_925_oxygen_pH_spore.csv, the NCBI taxonomy dump, the
   literature marker CSVs' notes/citation columns, and the sample metadata
   CSV's MeSH ids.
3. Vocabulary normalization — collapse case and synonym variants in the
   controlled-vocabulary trait columns ('Anaerobe'/'anaerobe'/'Anaerobic' are one
   value), and move cell-arrangement text out of `sporulation`.

Nothing is invented: every filled cell traces to a source file or to an
unambiguous relationship inside the database. Cells with no source stay NULL,
and `coverage_before_after.csv` reports what remains empty and why.

Idempotent — safe to re-run on an already-cleaned database.
"""
from __future__ import annotations

import csv
import glob
import os
import re
import sqlite3
import sys
from collections import Counter, defaultdict

import pandas as pd

# Repo root: see build_gutdb_sqlite.py — derived from this file's location,
# overridable with GUTDB_REPO.
REPO = os.environ.get(
    "GUTDB_REPO",
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
)
sys.path.insert(0, REPO)
from gutdb import transform as T  # noqa: E402  (needs REPO on the path)

DB = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")
# NCBI taxonomy dump (names.dmp + nodes.dmp), unzipped; see the README.
TAXDUMP = os.environ.get("GUTDB_TAXDUMP", "taxdump")

# --- controlled vocabularies -------------------------------------------------
# Case and synonym variants collapse onto one lowercase canonical form. The
# right-hand sides are the vocabulary; anything not listed is left verbatim so
# no information is silently dropped.
OXYGEN = {
    "aerobe": "aerobe", "aerobic": "aerobe", "obligate aerobe": "aerobe",
    "anaerobe": "anaerobe", "anaerobic": "anaerobe", "obligate anaerobe": "anaerobe",
    "facultative": "facultative anaerobe", "facultative anaerobe": "facultative anaerobe",
    "facultative aerobe": "facultative anaerobe",
    "aerobe/facultative anaerobe": "facultative anaerobe",
    "aerobe; facultative anaerobe": "facultative anaerobe",
    "microaerophile": "microaerophile", "microaerophilic": "microaerophile",
    "aerobe; anaerobe": "variable",
    "obligate aerobic": "aerobe", "facultatively anaerobe": "facultative anaerobe",
    "facultative anaerobe, air+co2": "facultative anaerobe",
}
GRAM = {
    "positive": "positive", "negative": "negative", "variable": "variable",
    "uncharacterized": "uncharacterized",
    # A cell wall that is structurally positive but does not hold the stain is
    # reported either way in the literature; "variable" is the honest label.
    "structurally positive but stains negative": "variable",
    "structurally positive but may stain negative": "variable",
    "negative due to the absence of a cell wall": "negative",
}
TEMPERATURE = {
    "mesophilic": "mesophilic", "mesophil": "mesophilic", "thermophilic": "thermophilic",
    "hyperthermophilic": "hyperthermophilic", "psychrophilic": "psychrophilic",
    "psychrotolerant": "psychrotolerant",
}
YES_NO = {"yes": "yes", "no": "no", "yes?": "yes", "no?": "no", "1": "yes", "0": "no",
          "true": "yes", "false": "no", "not known": ""}
PH = {"alkaliphile": "alkaliphile", "acidophile": "acidophile", "neutrophile": "neutrophile"}
# `superkingdom` collects three different vocabularies: MiMeDB's domains, the
# curated table's `kingdom` (Eubacteria, Fungi), and NCBI's post-2024 kingdom
# names (Bacillati, Pseudomonadati). All fold onto the three domains + Viruses.
SUPERKINGDOM = {
    "bacteria": "Bacteria", "eubacteria": "Bacteria", "bacillati": "Bacteria",
    "pseudomonadati": "Bacteria", "fusobacteriati": "Bacteria", "thermotogati": "Bacteria",
    "archaea": "Archaea", "methanobacteriati": "Archaea", "thermoproteati": "Archaea",
    "eukaryota": "Eukaryota", "fungi": "Eukaryota", "metazoa": "Eukaryota",
    "viridiplantae": "Eukaryota", "viruses": "Viruses",
}

ARRANGEMENT_TERMS = {"singles", "pairs", "chains", "clusters", "tetrads", "filaments",
                     "v-shaped forms"}


def canon(value, mapping: dict[str, str] | None = None) -> str:
    """Clean a source value; map it through a vocabulary when one is given.

    Under a vocabulary, a value the vocabulary does not list is kept but
    lowercased, so 'Aerotolerant' and 'aerotolerant' stop being two values while
    a category nobody anticipated is still never discarded.
    """
    text = T.clean_str(value).replace("\xa0", " ")
    text = re.sub(r"\s+", " ", text).strip()
    if not text:
        return ""
    if mapping is None:
        return text
    return mapping.get(text.casefold(), text if mapping is SUPERKINGDOM else text.casefold())


def split_spore(value) -> tuple[str, str]:
    """Return (sporulation, cell_arrangement) from one mixed `Spore`-style value."""
    text = canon(value)
    if not text:
        return "", ""
    folded = text.casefold()
    if folded in YES_NO:
        return YES_NO[folded], ""
    if "sporulat" in folded or "spotulating" in folded:
        return ("no" if folded.startswith("non") else "yes"), ""
    if "budding" in folded:
        return ("no" if folded.startswith(("no;", "no ")) else "yes"), ""
    tokens = {p.strip() for p in folded.replace("/", "-").split("-")}
    if tokens and tokens <= ARRANGEMENT_TERMS:
        return "", text
    return "", ""


# --- coverage reporting ------------------------------------------------------
def coverage(con: sqlite3.Connection) -> pd.DataFrame:
    tables = [r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'")]
    rows = []
    for table in tables:
        n = con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        cols = [c[1] for c in con.execute(f"PRAGMA table_info({table})")]
        if not n:
            continue
        expr = ", ".join(
            f"SUM(CASE WHEN {c} IS NULL OR TRIM(CAST({c} AS TEXT))='' THEN 1 ELSE 0 END)"
            for c in cols)
        for col, missing in zip(cols, con.execute(f"SELECT {expr} FROM {table}").fetchone()):
            rows.append({"table": table, "column": col, "n_rows": n, "n_missing": missing})
    return pd.DataFrame(rows)


# --- step 1: diseases -------------------------------------------------------
def source_mesh_map() -> dict[str, str]:
    """disease/phenotype name (casefolded) -> MeSH id, from every source CSV."""
    mapping: dict[str, str] = {}
    files = [os.path.join(REPO, "data", f) for f in (
        "gmrepo_all_disease_associations.csv", "gmrepo_ibs_additional_projects.csv",
        "gmrepo_PRJNA705217_ibs_markers.csv", "gmrepo_project_samples.csv")]
    files += sorted(glob.glob(os.path.join(REPO, "data", "literature_*_markers.csv")))
    pairs = [("disease_name", "mesh_id"), ("phenotype_a", "phenotype_a_mesh_id"),
             ("phenotype_b", "phenotype_b_mesh_id")]
    for path in files:
        if not os.path.exists(path):
            continue
        with open(path, newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                for name_col, mesh_col in pairs:
                    name, mesh = canon(row.get(name_col)), canon(row.get(mesh_col))
                    if name and mesh:
                        mapping.setdefault(name.casefold(), mesh)
    return mapping


def fill_disease_mesh(con: sqlite3.Connection) -> tuple[int, list[tuple[int, int]]]:
    """Fill missing disease MeSH ids from the source CSVs.

    A name whose MeSH id is already held by another row is not a fillable gap but
    a duplicate concept: 'Crohn's disease' resolves to D003424, which 'Crohn
    Disease' already owns. Those are returned as (dup_id, keeper_id) pairs for the
    merge step rather than written, since mesh_id is unique.
    """
    mapping = source_mesh_map()
    taken = {mesh: did for did, mesh in con.execute(
        "SELECT id, mesh_id FROM diseases WHERE mesh_id IS NOT NULL AND TRIM(mesh_id) <> ''")}
    filled, duplicates = 0, []
    for did, name in con.execute(
            "SELECT id, name FROM diseases WHERE mesh_id IS NULL OR TRIM(mesh_id)=''").fetchall():
        mesh = mapping.get(canon(name).casefold())
        if not mesh:
            continue
        if mesh in taken and taken[mesh] != did:
            duplicates.append((did, taken[mesh]))
        else:
            con.execute("UPDATE diseases SET mesh_id = ? WHERE id = ?", (mesh, did))
            taken[mesh] = did
            filled += 1
    con.commit()
    return filled, duplicates


def absorb_disease(con: sqlite3.Connection, dup: int, keeper: int) -> tuple[str, str, int, int]:
    """Repoint every reference from `dup` to `keeper`, then delete `dup`."""
    dup_name = con.execute("SELECT name FROM diseases WHERE id = ?", (dup,)).fetchone()[0]
    keeper_name = con.execute("SELECT name FROM diseases WHERE id = ?", (keeper,)).fetchone()[0]
    # Association rows that would collide with a keeper row on the evidence
    # unique key are duplicates of it, so they go rather than move.
    con.execute("""
        DELETE FROM taxon_disease_associations WHERE disease_id = ? AND id IN (
            SELECT d.id FROM taxon_disease_associations d
            JOIN taxon_disease_associations k
              ON k.disease_id = ? AND k.taxon_id = d.taxon_id
             AND k.comparison_id = d.comparison_id AND k.effect_type = d.effect_type
            WHERE d.disease_id = ?)
    """, (dup, keeper, dup))
    moved = con.execute(
        "UPDATE taxon_disease_associations SET disease_id = ? WHERE disease_id = ?",
        (keeper, dup)).rowcount
    moved_samples = con.execute(
        "UPDATE samples SET disease_id = ? WHERE disease_id = ?", (keeper, dup)).rowcount
    for col in ("phenotype_a_id", "phenotype_b_id",
                "positive_score_enriched_in_id", "negative_score_enriched_in_id"):
        con.execute(f"UPDATE phenotype_comparisons SET {col} = ? WHERE {col} = ?", (keeper, dup))
    con.execute("DELETE FROM diseases WHERE id = ?", (dup,))
    con.commit()
    return dup_name, keeper_name, max(moved, 0), max(moved_samples, 0)


def merge_duplicate_diseases(con: sqlite3.Connection) -> list[tuple[str, str, int, int]]:
    """Merge disease rows sharing a MeSH id; keep the best-attested name.

    'Crohn Disease' and \"Crohn's disease\" arrive from different sources as
    separate rows with the same MeSH concept, which splits their associations and
    breaks every per-disease rollup. The keeper is the row with the most
    associations (the MeSH-preferred spelling in practice); references in
    associations, samples and comparisons are repointed, and association rows
    that would collide on the evidence unique key after repointing are dropped as
    the duplicates they are.
    """
    merges = []
    groups = con.execute("""
        SELECT mesh_id, GROUP_CONCAT(id) FROM diseases
        WHERE mesh_id IS NOT NULL AND TRIM(mesh_id) <> ''
        GROUP BY mesh_id HAVING COUNT(*) > 1
    """).fetchall()
    for mesh, ids in groups:
        candidates = [int(i) for i in ids.split(",")]
        counts = {
            i: con.execute("SELECT COUNT(*) FROM taxon_disease_associations WHERE disease_id = ?",
                           (i,)).fetchone()[0] for i in candidates}
        keeper = max(candidates, key=lambda i: (counts[i], -i))
        for dup in candidates:
            if dup != keeper:
                merges.append(absorb_disease(con, dup, keeper))
    return merges


# --- step 2: taxon trait and lineage back-fill ------------------------------
TRAIT_COLUMNS = [
    "superkingdom", "phylum", "class_name", "order_name", "family", "ncbi_tax_id",
    "gram_stain", "oxygen_requirement", "ph_preference", "sporulation", "shape",
    "cell_arrangement", "mobility", "flagella_presence", "number_of_membranes",
    "biotic_relationship", "habitat", "temperature_range", "optimal_temperature",
    "metabolism", "energy_source", "human_pathogen",
]


def load_local_trait_sources() -> tuple[dict, dict]:
    """Return (species_records, genus_records) keyed by taxon key.

    species: (genus_key, species_key) -> {column: value}
    genus:   genus_key -> {column: value}, only where every source species in
             that genus agrees, so a genus never inherits a trait its members
             disagree about.
    """
    species: dict[tuple[str, str], dict[str, str]] = defaultdict(dict)
    genus_votes: dict[str, dict[str, Counter]] = defaultdict(lambda: defaultdict(Counter))
    genus_direct: dict[str, dict[str, str]] = defaultdict(dict)

    def put(target: dict, column: str, value: str) -> None:
        if value != "" and value is not None and column not in target:
            target[column] = value

    # (a) the project's own curated taxon table
    final_taxa = pd.read_csv(os.path.join(REPO, "gutdb", "final_taxa.tsv"), sep="\t")
    for row in final_taxa.itertuples():
        spor, arrangement = split_spore(getattr(row, "sporulation", ""))
        values = {
            "phylum": canon(T.normalize_phylum(row.phylum)),
            "class_name": canon(row.class_name), "order_name": canon(row.order_name),
            "family": canon(row.family), "superkingdom": canon(getattr(row, "kingdom", ""), SUPERKINGDOM),
            "ncbi_tax_id": canon(T.nullable_int(row.ncbi_tax_id)),
            "gram_stain": canon(row.gram_stain, GRAM),
            "oxygen_requirement": canon(row.oxygen_requirement, OXYGEN),
            "ph_preference": canon(row.ph_preference, PH),
            "sporulation": spor, "cell_arrangement": arrangement,
        }
        gk, sk = canon(row.genus_key).casefold(), canon(row.species_key).casefold()
        if row.taxonomic_rank == "species" and gk and sk:
            for col, val in values.items():
                put(species[(gk, sk)], col, val)
                if val:
                    genus_votes[gk][col][val] += 1
        elif gk:
            for col, val in values.items():
                put(genus_direct[gk], col, val)

    # (b) MiMeDB microbe attributes
    mimedb = pd.read_csv(os.path.join(REPO, "mimedb_microbes_v2.csv"))
    for row in mimedb.itertuples():
        gk, sk = T.taxon_key(row.genus, row.species)
        spor, arrangement = split_spore(row.sporulation)
        membranes = T.nullable_int(row.number_of_membranes)
        values = {
            "superkingdom": canon(row.superkingdom, SUPERKINGDOM),
            "phylum": canon(T.normalize_phylum(row.phylum)),
            "class_name": canon(row.klass), "order_name": canon(getattr(row, "order")),
            "family": canon(row.family),
            "ncbi_tax_id": canon(T.nullable_int(row.ncbi_tax_id)),
            "gram_stain": canon(row.gram, GRAM),
            "oxygen_requirement": canon(row.oxygen_requirement, OXYGEN),
            "sporulation": spor, "cell_arrangement": arrangement or canon(row.cell_arrangement),
            "shape": canon(row.shape), "mobility": canon(row.mobility, YES_NO),
            "flagella_presence": canon(row.flagella_presence, YES_NO),
            "number_of_membranes": canon(membranes),
            "biotic_relationship": canon(row.biotic_relationship),
            "habitat": canon(row.habitat),
            "temperature_range": canon(row.temperature_range, TEMPERATURE),
            "optimal_temperature": canon(row.optimal_temperature),
            "metabolism": canon(row.metabolism), "energy_source": canon(row.energy_source),
            "human_pathogen": canon(T.normalize_pathogen_flag(T.nullable_int(row.human_pathogen))),
        }
        if gk and sk:
            for col, val in values.items():
                put(species[(gk, sk)], col, val)
                if val:
                    genus_votes[gk][col][val] += 1
        elif gk:
            for col, val in values.items():
                put(genus_direct[gk], col, val)

    # (c) the 925-taxon oxygen/pH/spore table
    extra = pd.read_csv(os.path.join(REPO, "gut_microbiome_925_oxygen_pH_spore.csv"))
    for row in extra.itertuples():
        gk, sk = T.taxon_key(row.Genus, row.Species)
        spor, arrangement = split_spore(row.Spore)
        values = {
            "oxygen_requirement": canon(row.Oxygen, OXYGEN),
            "ph_preference": canon(row.pH, PH),
            "sporulation": spor, "cell_arrangement": arrangement,
            "phylum": canon(T.normalize_phylum(row.Phylum)),
            "class_name": canon(row.Class), "order_name": canon(row.Order),
            "family": canon(row.Family),
        }
        if gk and sk:
            for col, val in values.items():
                put(species[(gk, sk)], col, val)
                if val:
                    genus_votes[gk][col][val] += 1
        elif gk:
            for col, val in values.items():
                put(genus_direct[gk], col, val)

    # genus-level values: direct genus rows win; otherwise a unanimous vote
    genus: dict[str, dict[str, str]] = {}
    for gk, columns in genus_votes.items():
        record = {}
        for col, votes in columns.items():
            if len(votes) == 1:
                record[col] = next(iter(votes))
        genus[gk] = record
    for gk, record in genus_direct.items():
        genus.setdefault(gk, {}).update(
            {c: v for c, v in record.items() if c not in genus.get(gk, {})})
    # A genus has no species epithet or species-specific tax id of its own.
    for record in genus.values():
        record.pop("ncbi_tax_id", None)
        record.pop("human_pathogen", None)
    return species, genus


def fill_taxa_from_local(con: sqlite3.Connection) -> Counter:
    species, genus = load_local_trait_sources()
    filled = Counter()
    rows = con.execute(
        "SELECT id, genus_key, species_key, taxonomic_rank, "
        + ", ".join(TRAIT_COLUMNS) + " FROM taxa").fetchall()
    for row in rows:
        taxon_id, gk, sk, rank = row[0], row[1] or "", row[2] or "", row[3]
        current = dict(zip(TRAIT_COLUMNS, row[4:]))
        record = species.get((gk, sk), {}) if rank == "species" else {}
        fallback = genus.get(gk, {})
        updates = {}
        for col in TRAIT_COLUMNS:
            if current.get(col) not in (None, ""):
                continue
            value = record.get(col) or fallback.get(col)
            if value:
                updates[col] = value
        if updates:
            con.execute(f"UPDATE taxa SET {', '.join(c + ' = ?' for c in updates)} WHERE id = ?",
                        (*updates.values(), taxon_id))
            filled.update(updates.keys())
    con.commit()
    return filled


# --- step 3: NCBI taxonomy dump --------------------------------------------
RANK_COLUMN = {"superkingdom": "superkingdom", "domain": "superkingdom", "phylum": "phylum",
               "class": "class_name", "order": "order_name", "family": "family"}


def fill_from_ncbi_taxonomy(con: sqlite3.Connection) -> Counter:
    """Fill tax ids and missing lineage ranks from the NCBI taxonomy dump.

    Authoritative for lineage, which the association exports do not carry at all.
    Names are matched on the NCBI scientific name first, then on synonyms and
    equivalent names, so taxa renamed since the export still resolve. Phylum
    names pass through gutdb.transform.normalize_phylum so post-2021 NCBI names
    (Bacillota, Bacteroidota) fold onto the classic names already in the table.
    """
    names_path, nodes_path = f"{TAXDUMP}/names.dmp", f"{TAXDUMP}/nodes.dmp"
    if not (os.path.exists(names_path) and os.path.exists(nodes_path)):
        return Counter()

    taxa = con.execute(
        "SELECT id, scientific_name, taxonomic_rank, ncbi_tax_id, superkingdom, phylum,"
        " class_name, order_name, family FROM taxa").fetchall()
    wanted = {canon(r[1]).casefold(): [] for r in taxa}
    for row in taxa:
        wanted[canon(row[1]).casefold()].append(row)

    # pass 1 — resolve our names to tax ids
    resolved: dict[str, int] = {}
    synonyms: dict[str, int] = {}
    with open(names_path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split("\t|")
            if len(parts) < 4:
                continue
            name = parts[1].strip().casefold()
            if name not in wanted:
                continue
            name_class = parts[3].strip()
            if name_class == "scientific name":
                resolved.setdefault(name, int(parts[0]))
            elif name_class in ("synonym", "equivalent name", "includes"):
                synonyms.setdefault(name, int(parts[0]))
    for name, taxid in synonyms.items():
        resolved.setdefault(name, taxid)

    # pass 2 — parent/rank for every node (array-indexed: the dump is 2.6M nodes)
    parent: dict[int, int] = {}
    rank: dict[int, str] = {}
    with open(nodes_path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split("\t|")
            taxid, par, rnk = int(parts[0]), int(parts[1]), parts[2].strip()
            parent[taxid] = par
            if rnk in RANK_COLUMN:
                rank[taxid] = rnk

    # ancestors we need names for
    needed: set[int] = set()
    chains: dict[int, list[int]] = {}
    for taxid in set(resolved.values()):
        chain, node, guard = [], taxid, 0
        while node and node != 1 and guard < 60:
            if node in rank:
                chain.append(node)
                needed.add(node)
            node = parent.get(node, 0)
            guard += 1
        chains[taxid] = chain

    # pass 3 — scientific names for those ancestors
    ancestor_name: dict[int, str] = {}
    with open(names_path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split("\t|")
            if len(parts) < 4 or parts[3].strip() != "scientific name":
                continue
            taxid = int(parts[0])
            if taxid in needed:
                ancestor_name[taxid] = parts[1].strip()

    filled = Counter()
    for name, taxid in resolved.items():
        lineage = {}
        for node in chains.get(taxid, []):
            column = RANK_COLUMN[rank[node]]
            value = ancestor_name.get(node, "")
            if column == "phylum":
                value = T.normalize_phylum(value)
            if column == "superkingdom":
                # walking outward, so the last (root-most) domain node wins
                lineage[column] = canon(value, SUPERKINGDOM)
            else:
                lineage.setdefault(column, value)
        for row in wanted.get(name, []):
            taxon_id, _, taxon_rank, current_taxid = row[0], row[1], row[2], row[3]
            current = dict(zip(("superkingdom", "phylum", "class_name", "order_name", "family"),
                               row[4:]))
            updates = {c: v for c, v in lineage.items() if v and current.get(c) in (None, "")}
            # Only a species row gets a species-level tax id; a genus row keeps
            # the genus id it resolved to.
            if current_taxid in (None, "") and (
                    taxon_rank == "species" or rank.get(taxid) == "genus"):
                updates["ncbi_tax_id"] = taxid
            if updates:
                con.execute(
                    f"UPDATE taxa SET {', '.join(c + ' = ?' for c in updates)} WHERE id = ?",
                    (*updates.values(), taxon_id))
                filled.update(updates.keys())
    con.commit()
    return filled


def propagate_genus_lineage(con: sqlite3.Connection) -> int:
    """Give genus rows the modal lineage of the species indexed under them."""
    filled = 0
    for level in ("superkingdom", "phylum", "class_name", "order_name", "family"):
        rows = con.execute(
            f"SELECT genus_key, {level}, COUNT(*) FROM taxa WHERE taxonomic_rank = 'species'"
            f" AND {level} IS NOT NULL AND TRIM({level}) <> ''"
            " GROUP BY genus_key, 2 ORDER BY genus_key, 3 DESC").fetchall()
        modal: dict[str, str] = {}
        for genus_key, value, _ in rows:
            modal.setdefault(genus_key, value)
        cur = con.executemany(
            f"UPDATE taxa SET {level} = ? WHERE genus_key = ? AND taxonomic_rank = 'genus'"
            f" AND ({level} IS NULL OR TRIM({level}) = '')",
            [(v, k) for k, v in modal.items()])
        filled += max(cur.rowcount, 0)
    con.commit()
    return filled


# --- step 4: vocabulary normalization on rows already in the table ----------
def normalize_existing_values(con: sqlite3.Connection) -> Counter:
    changed = Counter()
    specs = [("oxygen_requirement", OXYGEN), ("gram_stain", GRAM), ("ph_preference", PH),
             ("superkingdom", SUPERKINGDOM),
             ("temperature_range", TEMPERATURE), ("mobility", YES_NO),
             ("flagella_presence", YES_NO), ("human_pathogen", None),
             ("shape", None), ("metabolism", None), ("habitat", None),
             ("energy_source", None), ("biotic_relationship", None),
             ("optimal_temperature", None), ("cell_arrangement", None)]
    for column, mapping in specs:
        for taxon_id, value in con.execute(
                f"SELECT id, {column} FROM taxa WHERE {column} IS NOT NULL"
                f" AND TRIM({column}) <> ''"):
            new = canon(value, mapping)
            # Compare cleaned-to-cleaned: an INTEGER cell reads back as 1, not
            # '1', and would otherwise look changed on every run.
            if new != canon(value):
                con.execute(f"UPDATE taxa SET {column} = ? WHERE id = ?",
                            (new or None, taxon_id))
                changed[column] += 1
    # sporulation carries three different kinds of value; re-split it
    for taxon_id, value, arrangement in con.execute(
            "SELECT id, sporulation, cell_arrangement FROM taxa"
            " WHERE sporulation IS NOT NULL AND TRIM(sporulation) <> ''"):
        spor, arr = split_spore(value)
        if spor != value or (arr and not arrangement):
            con.execute("UPDATE taxa SET sporulation = ?, cell_arrangement = COALESCE(?, "
                        "cell_arrangement) WHERE id = ?",
                        (spor or None, arr or None, taxon_id))
            changed["sporulation"] += 1
    con.commit()
    return changed


# --- step 5: samples and comparisons ---------------------------------------
def fill_samples(con: sqlite3.Connection) -> dict[str, int]:
    """Resolve missing sample disease labels via the source MeSH id; set body_site.

    `body_site` is absent from every source file but is not unknown: GMrepo
    curates human *gut* metagenomes, and the literature panel is a gut-microbiome
    set, so the column is a constant of the source's scope rather than
    per-sample metadata. It is filled as 'gut' and flagged as such here.
    """
    mesh_to_disease = {m: i for i, m in con.execute(
        "SELECT id, mesh_id FROM diseases WHERE mesh_id IS NOT NULL AND TRIM(mesh_id) <> ''")}
    by_run: dict[str, str] = {}
    with open(os.path.join(REPO, "data", "gmrepo_project_samples.csv"),
              newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            run = canon(row.get("run_id"))
            mesh = canon(row.get("mesh_id"))
            if run and mesh:
                by_run[run] = mesh
    out = {"disease_id": 0, "body_site": 0}
    for sample_id, run in con.execute(
            "SELECT id, run_accession FROM samples WHERE disease_id IS NULL"):
        disease_id = mesh_to_disease.get(by_run.get(canon(run), ""))
        if disease_id:
            con.execute("UPDATE samples SET disease_id = ? WHERE id = ?", (disease_id, sample_id))
            out["disease_id"] += 1
    out["body_site"] = max(con.execute(
        "UPDATE samples SET body_site = 'gut' WHERE body_site IS NULL"
        " OR TRIM(body_site) = ''").rowcount, 0)
    con.commit()
    return out


def fill_comparison_context(con: sqlite3.Connection) -> dict[str, int]:
    """Fill comparison `notes` and `method` from the literature marker CSVs.

    Those CSVs carry a per-study `notes` field (cohort caveats, how the p-value
    was obtained) and a `citation`, both of which the build discarded because the
    association table has nowhere to put them. `method` is taken from the source's
    `effect_type` only when that column names an actual method (LEfSe); rows whose
    effect_type is 'reported' get 'reported in publication', which is what the
    source says and no more.
    """
    per_project: dict[str, dict[str, str]] = {}
    for path in sorted(glob.glob(os.path.join(REPO, "data", "literature_*_markers.csv"))):
        with open(path, newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                project = canon(row.get("project_id"))
                if not project or project in per_project:
                    continue
                note = " | ".join(x for x in (canon(row.get("notes")),
                                              canon(row.get("citation")),
                                              canon(row.get("association_url"))) if x)
                per_project[project] = {"notes": note,
                                        "effect_type": canon(row.get("effect_type"))}
    out = {"notes": 0, "method": 0}
    for comparison_id, project, method, notes in con.execute("""
            SELECT c.id, s.project_accession, c.method, c.notes
            FROM phenotype_comparisons c JOIN studies s ON s.id = c.study_id"""):
        record = per_project.get(project)
        if not record:
            continue
        if record["notes"] and not notes:
            con.execute("UPDATE phenotype_comparisons SET notes = ? WHERE id = ?",
                        (record["notes"], comparison_id))
            out["notes"] += 1
        if not method:
            effect_type = record["effect_type"]
            value = ("reported in publication" if effect_type.casefold() == "reported"
                     else effect_type or None)
            if value:
                con.execute("UPDATE phenotype_comparisons SET method = ? WHERE id = ?",
                            (value, comparison_id))
                out["method"] += 1
    con.commit()
    return out


def main() -> None:
    con = sqlite3.connect(DB)
    con.execute("PRAGMA foreign_keys = ON")
    before = coverage(con)

    n_mesh, dup_pairs = fill_disease_mesh(con)
    print("disease mesh_id filled:", n_mesh)
    merges = [absorb_disease(con, dup, keeper) for dup, keeper in dup_pairs]
    merges += merge_duplicate_diseases(con)
    for dup, keeper, moved, moved_samples in merges:
        print(f"  merged {dup!r} -> {keeper!r} "
              f"({moved} associations, {moved_samples} samples repointed)")

    print("taxa cells filled from local sources:", dict(fill_taxa_from_local(con)))
    print("taxa cells filled from NCBI taxonomy:", dict(fill_from_ncbi_taxonomy(con)))
    print("genus lineage cells propagated:", propagate_genus_lineage(con))
    print("values normalized:", dict(normalize_existing_values(con)))
    print("samples filled:", fill_samples(con))
    print("comparison context filled:", fill_comparison_context(con))

    after = coverage(con)
    report = before.merge(after, on=["table", "column", "n_rows"],
                          suffixes=("_before", "_after"))
    report["n_filled"] = report.n_missing_before - report.n_missing_after
    report["pct_missing_after"] = (100 * report.n_missing_after / report.n_rows).round(1)
    report = report.sort_values(["n_filled", "n_missing_after"], ascending=[False, False])
    report.to_csv("coverage_before_after.csv", index=False)
    print("\nwrote coverage_before_after.csv")
    print(report[report.n_filled > 0].to_string(index=False))

    print("\nforeign_key_check violations:", len(con.execute("PRAGMA foreign_key_check").fetchall()))
    con.execute("ANALYZE")
    con.commit()
    con.close()


if __name__ == "__main__":
    main()
