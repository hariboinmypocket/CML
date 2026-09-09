"""Build a SQLite mirror of the gut-microbiome-disease database from the local CSVs.

The project's own pipeline targets MySQL 8 (gutdb/db.py, gutdb/schema.sql); this
script loads the same sources into the same normalized schema in a single-file
SQLite database so the data can be queried without a server. Normalization is not
re-implemented: gutdb.transform is imported and reused, so taxon keys, phylum
folding, effect-size quantization and direction derivation match the MySQL loads.

Idempotent: rerunning rebuilds the file from scratch.
"""
from __future__ import annotations

import csv
import glob
import os
import sqlite3
import sys
from datetime import datetime, timezone

# Repo root: the directory this script's parent lives in (scripts/ -> repo root),
# overridable with GUTDB_REPO when the checkout is somewhere else.
PROJECT = os.environ.get(
    "GUTDB_REPO",
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
)
sys.path.insert(0, PROJECT)
csv.field_size_limit(10_000_000)

from gutdb import transform as T  # noqa: E402

DB_PATH = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")
SCHEMA = os.environ.get("GUTDB_SCHEMA", os.path.join(PROJECT, "gutdb", "schema_sqlite.sql"))

TAXA_BASE = os.path.join(PROJECT, "gut_microbiome_new_full.csv")
ASSOC_FILES = [
    os.path.join(PROJECT, "data/gmrepo_all_disease_associations.csv"),
    os.path.join(PROJECT, "data/gmrepo_ibs_additional_projects.csv"),
    os.path.join(PROJECT, "data/gmrepo_PRJNA705217_ibs_markers.csv"),
] + sorted(glob.glob(os.path.join(PROJECT, "data/literature_*_markers.csv")))
SAMPLES_FILE = os.path.join(PROJECT, "data/gmrepo_project_samples.csv")

# Control arm of a literature marker comparison, dropped from those CSVs as
# constant columns; see Loader.comparison.
CONTROL_PHENOTYPE = "Health"
CONTROL_MESH_ID = "D006262"
ABUNDANCE_FILE = os.path.join(PROJECT, "data/gmrepo_species_abundances.csv")


ENERGY_MODE_SYNONYMS = {"respirer": "respirator"}
# The `Spore` column in gut_microbiome_new_full.csv carries three different
# things: a 0/1 sporulation flag, free-text sporulation phrases, and
# cell-arrangement values (Singles/Pairs/Chains/...) that belong in
# cell_arrangement. Route each to the column it describes rather than writing
# arrangement text into sporulation.
ARRANGEMENT_TERMS = {
    "singles", "pairs", "chains", "clusters", "tetrads", "filaments",
    "v-shaped forms",
}


def classify_spore(value: str) -> tuple[str, str]:
    """Return (sporulation, cell_arrangement) for one raw `Spore` value."""
    text = T.clean_str(value)
    if not text:
        return "", ""
    folded = text.casefold()
    if folded == "1":
        return "Yes", ""
    if folded == "0":
        return "No", ""
    if "sporulat" in folded or "spotulating" in folded or "budding" in folded:
        return text, ""
    tokens = {part.strip() for part in folded.replace("/", "-").split("-")}
    if tokens and tokens <= ARRANGEMENT_TERMS:
        return "", text
    return "", ""


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def read_rows(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            yield row


class Loader:
    def __init__(self, con: sqlite3.Connection):
        self.con = con
        self.disease_ids: dict[str, int] = {}
        self.study_ids: dict[str, int] = {}
        self.taxon_ids: dict[tuple[str, str], int] = {}
        self.comparison_ids: dict[str, int] = {}
        self.sample_ids: dict[str, int] = {}
        self.comparison_ranks: dict[int, set[str]] = {}

    # ---------- audit ----------
    def start_run(self, source_name: str, uri: str) -> int:
        cur = self.con.execute(
            "INSERT INTO ingestion_runs (source_name, source_uri, started_at, status)"
            " VALUES (?, ?, ?, 'running')",
            (source_name, uri, now()),
        )
        return int(cur.lastrowid)

    def finish_run(self, run_id: int, read: int, ins: int, upd: int, skip: int) -> None:
        self.con.execute(
            "UPDATE ingestion_runs SET completed_at = ?, status = 'completed',"
            " rows_read = ?, rows_inserted = ?, rows_updated = ?, rows_skipped = ?"
            " WHERE id = ?",
            (now(), read, ins, upd, skip, run_id),
        )
        self.con.commit()

    # ---------- dimension upserts ----------
    def disease(self, name, mesh_id="") -> int | None:
        name = T.clean_str(name)
        if not name:
            return None
        key = name.casefold()
        mesh = T.clean_str(mesh_id) or None
        if key in self.disease_ids:
            did = self.disease_ids[key]
            if mesh:
                row = self.con.execute(
                    "SELECT mesh_id FROM diseases WHERE id = ?", (did,)
                ).fetchone()
                if row[0] is None:
                    taken = self.con.execute(
                        "SELECT 1 FROM diseases WHERE mesh_id = ?", (mesh,)
                    ).fetchone()
                    if not taken:
                        self.con.execute(
                            "UPDATE diseases SET mesh_id = ? WHERE id = ?", (mesh, did)
                        )
            return did
        if mesh and self.con.execute(
            "SELECT 1 FROM diseases WHERE mesh_id = ?", (mesh,)
        ).fetchone():
            mesh = None  # MeSH already claimed by another spelling; keep name unique
        cur = self.con.execute(
            "INSERT INTO diseases (mesh_id, name, name_key) VALUES (?, ?, ?)",
            (mesh, name, key),
        )
        self.disease_ids[key] = int(cur.lastrowid)
        return self.disease_ids[key]

    def study(self, accession, title="", description="", data_type="",
              data_quality="", source_database="GMrepo") -> int | None:
        acc = T.clean_str(accession)
        if not acc:
            return None
        quality = T.clean_str(data_quality).casefold()
        quality = quality if quality in {"curated", "qualified"} else "unknown"
        if acc in self.study_ids:
            sid = self.study_ids[acc]
            self.con.execute(
                "UPDATE studies SET title = COALESCE(NULLIF(?, ''), title),"
                " description = COALESCE(NULLIF(?, ''), description),"
                " data_type = COALESCE(NULLIF(?, ''), data_type),"
                " data_quality = CASE WHEN data_quality = 'unknown' THEN ? ELSE data_quality END"
                " WHERE id = ?",
                (T.clean_str(title), T.clean_str(description), T.clean_str(data_type),
                 quality, sid),
            )
            return sid
        cur = self.con.execute(
            "INSERT INTO studies (project_accession, title, description, data_type,"
            " source_database, data_quality) VALUES (?, ?, ?, ?, ?, ?)",
            (acc, T.clean_str(title) or None, T.clean_str(description) or None,
             T.clean_str(data_type) or None, T.clean_str(source_database) or "GMrepo",
             quality),
        )
        self.study_ids[acc] = int(cur.lastrowid)
        return self.study_ids[acc]

    def taxon_from_name(self, scientific_name, rank="", ncbi_tax_id="",
                        source_database="GMrepo") -> int | None:
        genus, species = T.parse_scientific_name(scientific_name)
        if not genus:
            return None
        return self.taxon(genus, species, rank=rank, ncbi_tax_id=ncbi_tax_id,
                          source_database=source_database)

    def taxon(self, genus, species, rank="", ncbi_tax_id="", source_database="GMrepo",
              extra: dict | None = None) -> int | None:
        genus = T.normalize_genus(genus)
        if not genus:
            return None
        species = T.extract_species_epithet(genus, species)
        key = (genus.casefold(), species.casefold())
        rank = T.clean_str(rank).casefold()
        if rank not in {"species", "genus", "other"}:
            rank = "species" if species else "genus"
        tax_id = T.nullable_int(ncbi_tax_id)
        if key in self.taxon_ids:
            tid = self.taxon_ids[key]
            if tax_id is not None:
                self.con.execute(
                    "UPDATE taxa SET ncbi_tax_id = COALESCE(ncbi_tax_id, ?) WHERE id = ?",
                    (tax_id, tid),
                )
            if extra:
                self._fill_taxon(tid, extra)
            return tid
        cols = {
            "genus": genus,
            "species": species,
            "genus_key": key[0],
            "species_key": key[1],
            "scientific_name": f"{genus} {species}".strip(),
            "taxonomic_rank": rank,
            "ncbi_tax_id": tax_id,
            "source_database": T.clean_str(source_database) or "GMrepo",
        }
        if extra:
            cols.update({k: v for k, v in extra.items() if v not in ("", None)})
        names = ", ".join(cols)
        cur = self.con.execute(
            f"INSERT INTO taxa ({names}) VALUES ({', '.join('?' * len(cols))})",
            tuple(cols.values()),
        )
        self.taxon_ids[key] = int(cur.lastrowid)
        return self.taxon_ids[key]

    def _fill_taxon(self, taxon_id: int, extra: dict) -> None:
        usable = {k: v for k, v in extra.items() if v not in ("", None)}
        if not usable:
            return
        sets = ", ".join(f"{k} = COALESCE({k}, ?)" for k in usable)
        self.con.execute(f"UPDATE taxa SET {sets} WHERE id = ?",
                         (*usable.values(), taxon_id))

    def comparison(self, study_id, row) -> int:
        # The literature marker CSVs no longer carry the phenotype columns: every
        # one of their comparisons was disease-versus-healthy, so phenotype_a was
        # the constant "Health"/D006262 and phenotype_b restated disease_name/
        # mesh_id. The defaults below reconstruct them. The full GMrepo export
        # still carries real phenotype pairs, and supplies them, so no default
        # fires on that path.
        pheno_a = T.first_value(row, "phenotype_a") or CONTROL_PHENOTYPE
        pheno_b = T.first_value(row, "phenotype_b") or T.first_value(row, "disease_name")
        a_mesh = T.first_value(row, "phenotype_a_mesh_id") or (
            CONTROL_MESH_ID if pheno_a == CONTROL_PHENOTYPE else None)
        a_id = self.disease(pheno_a, a_mesh)
        b_id = self.disease(pheno_b, T.first_value(row, "phenotype_b_mesh_id")
                            or T.first_value(row, "mesh_id"))
        method = T.first_value(row, "method")
        effect_type = T.first_value(row, "effect_type") or "LDA"
        # A comparison is identified by study + the two phenotype groups + effect
        # type. `method` is deliberately NOT part of the key: the standalone
        # PRJNA705217 marker export omits the column while the full GMrepo export
        # records "LEfSe" for the same rows, and keying on it would split one
        # comparison in two and double-count all 34 associations.
        key = "|".join([
            T.clean_str(row.get("project_id")),
            (pheno_a or "").casefold(),
            (pheno_b or "").casefold(),
            effect_type.casefold(),
        ])
        if key in self.comparison_ids:
            comparison_id = self.comparison_ids[key]
            if method:
                self.con.execute(
                    "UPDATE phenotype_comparisons SET method = COALESCE(method, ?)"
                    " WHERE id = ?",
                    (method, comparison_id),
                )
            return comparison_id
        pos = self.disease(T.first_value(row, "positive_enriched_in"))
        neg = self.disease(T.first_value(row, "negative_enriched_in"))
        cur = self.con.execute(
            "INSERT INTO phenotype_comparisons (study_id, phenotype_a_id, phenotype_b_id,"
            " comparison_key, method, taxonomic_level, positive_score_enriched_in_id,"
            " negative_score_enriched_in_id) VALUES (?, ?, ?, ?, ?, 'species', ?, ?)",
            (study_id, a_id, b_id, key, method or None, pos, neg),
        )
        self.comparison_ids[key] = int(cur.lastrowid)
        return self.comparison_ids[key]

    # ---------- fact loaders ----------
    def load_base_taxa(self) -> None:
        run = self.start_run("gut_microbiome_new_full", TAXA_BASE)
        read = ins = 0
        for row in read_rows(TAXA_BASE):
            read += 1
            sporulation, arrangement = classify_spore(row.get("Spore"))
            energy_mode = T.clean_str(row.get("energy_mode")).casefold()
            energy_mode = ENERGY_MODE_SYNONYMS.get(energy_mode, energy_mode)
            if energy_mode not in {"fermenter", "respirator", "mixed"}:
                energy_mode = "N/A"
            extra = {
                "phylum": T.normalize_phylum(row.get("Phylum")),
                "class_name": T.clean_str(row.get("Class")),
                "order_name": T.clean_str(row.get("Order")),
                "family": T.clean_str(row.get("Family")),
                "oxygen_requirement": T.clean_str(row.get("Oxygen")),
                "ph_preference": T.clean_str(row.get("pH")),
                "sporulation": sporulation,
                "cell_arrangement": arrangement,
                "energy_mode": energy_mode,
                "primary_food_source": T.clean_str(row.get("primary_food_source")) or "N/A",
            }
            before = len(self.taxon_ids)
            self.taxon(row.get("Genus"), row.get("Species"),
                       source_database="curated_local", extra=extra)
            ins += len(self.taxon_ids) > before
        self.finish_run(run, read, ins, read - ins, 0)

    def load_associations(self, path: str) -> tuple[int, int, int]:
        run = self.start_run(f"associations:{os.path.basename(path)}", path)
        read = ins = skip = 0
        for row in read_rows(path):
            read += 1
            study_id = self.study(
                row.get("project_id"), row.get("project_title"),
                row.get("project_description"), row.get("data_type"),
                row.get("data_quality"),
                T.first_value(row, "source_database") or "GMrepo",
            )
            disease_name = T.first_value(row, "disease_name")
            disease_id = self.disease(disease_name, T.first_value(row, "mesh_id"))
            taxon_id = self.taxon_from_name(
                T.first_value(row, "scientific_name", "marker_taxon"),
                T.first_value(row, "taxonomic_rank"),
                T.first_value(row, "ncbi_tax_id"),
                T.first_value(row, "source_database") or "GMrepo",
            )
            if not (study_id and disease_id and taxon_id):
                skip += 1
                continue
            comparison_id = self.comparison(study_id, row)
            rank = self.con.execute(
                "SELECT taxonomic_rank FROM taxa WHERE id = ?", (taxon_id,)
            ).fetchone()[0]
            self.comparison_ranks.setdefault(comparison_id, set()).add(rank)
            lda = T.nullable_float(T.first_value(row, "lda_score", "effect_size"))
            direction = T.association_direction(
                lda, disease_name,
                T.first_value(row, "positive_enriched_in"),
                T.first_value(row, "negative_enriched_in"),
                T.first_value(row, "direction"),
            )
            cur = self.con.execute(
                "INSERT OR IGNORE INTO taxon_disease_associations (taxon_id, disease_id,"
                " comparison_id, direction, effect_type, effect_size, p_value, q_value,"
                " source_database) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (taxon_id, disease_id, comparison_id, direction,
                 T.first_value(row, "effect_type") or "LDA",
                 T.quantize_effect_size(lda),
                 T.nullable_float(T.first_value(row, "p_value", "pvalue")),
                 T.nullable_float(T.first_value(row, "q_value", "fdr")),
                 T.first_value(row, "source_database") or "GMrepo"),
            )
            if cur.rowcount:
                ins += 1
            else:
                skip += 1
        self.finish_run(run, read, ins, 0, skip)
        return read, ins, skip

    def finalize_comparison_levels(self) -> None:
        for comparison_id, ranks in self.comparison_ranks.items():
            level = "mixed" if len(ranks) > 1 else next(iter(ranks))
            if level not in {"species", "genus", "mixed"}:
                level = "mixed"
            self.con.execute(
                "UPDATE phenotype_comparisons SET taxonomic_level = ? WHERE id = ?",
                (level, comparison_id),
            )
        self.con.commit()

    def backfill_genus_lineage(self) -> int:
        """Give genus-rank taxa the lineage carried by their own species.

        Genus rows are created from the association exports, which ship no
        lineage columns, so `phylum` came out NULL for every genus — and the
        specificity/evidence views group by phylum. The curated species rows
        already carry a lineage, so each genus inherits the modal lineage of the
        species indexed under it. Genera with no curated species stay NULL rather
        than being guessed at.
        """
        filled = 0
        for level in ("phylum", "class_name", "order_name", "family"):
            rows = self.con.execute(
                f"SELECT genus_key, {level}, COUNT(*) AS n FROM taxa"
                f" WHERE taxonomic_rank = 'species' AND {level} IS NOT NULL"
                " GROUP BY genus_key, 2 ORDER BY genus_key, n DESC"
            ).fetchall()
            modal: dict[str, str] = {}
            for genus_key, value, _ in rows:
                modal.setdefault(genus_key, value)
            cur = self.con.executemany(
                f"UPDATE taxa SET {level} = ? WHERE genus_key = ?"
                f" AND taxonomic_rank = 'genus' AND {level} IS NULL",
                [(value, genus_key) for genus_key, value in modal.items()],
            )
            filled += cur.rowcount if cur.rowcount and cur.rowcount > 0 else 0
        self.con.commit()
        return filled

    def load_samples(self) -> tuple[int, int, int]:
        run = self.start_run("gmrepo_project_samples", SAMPLES_FILE)
        read = ins = skip = 0
        for row in read_rows(SAMPLES_FILE):
            read += 1
            study_id = self.study(row.get("project_id"), data_type=row.get("data_type"),
                                  data_quality=row.get("data_quality"))
            run_acc = T.clean_str(row.get("run_id"))
            if not (study_id and run_acc):
                skip += 1
                continue
            disease_id = self.disease(row.get("disease_name"), row.get("mesh_id"))
            cur = self.con.execute(
                "INSERT OR IGNORE INTO samples (study_id, disease_id, gmrepo_sample_id,"
                " run_accession, sex, age_years, bmi, country, qc_status)"
                " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (study_id, disease_id, T.clean_str(row.get("sample_id")) or None,
                 run_acc, T.clean_str(row.get("sex")) or None,
                 T.nullable_float(row.get("age_years")), T.nullable_float(row.get("bmi")),
                 T.clean_str(row.get("country")) or None,
                 T.clean_str(row.get("qc_status")) or None),
            )
            if cur.rowcount:
                ins += 1
                self.sample_ids[run_acc] = int(cur.lastrowid)
            else:
                skip += 1
        self.finish_run(run, read, ins, 0, skip)
        return read, ins, skip

    def load_abundances(self) -> tuple[int, int, int]:
        run = self.start_run("gmrepo_species_abundances", ABUNDANCE_FILE)
        read = ins = skip = 0
        batch: list[tuple] = []
        for row in read_rows(ABUNDANCE_FILE):
            read += 1
            sample_id = self.sample_ids.get(T.clean_str(row.get("run_id")))
            abundance = T.nullable_float(row.get("relative_abundance"))
            if sample_id is None or abundance is None:
                skip += 1
                continue
            taxon_id = self.taxon_from_name(
                row.get("scientific_name"), row.get("taxonomic_rank"),
                row.get("ncbi_tax_id"),
                T.clean_str(row.get("source_database")) or "GMrepo",
            )
            if taxon_id is None:
                skip += 1
                continue
            batch.append((sample_id, taxon_id, abundance))
            if len(batch) >= 20000:
                ins += self._flush_abundances(batch)
                batch.clear()
        ins += self._flush_abundances(batch)
        self.finish_run(run, read, ins, 0, read - ins)
        return read, ins, read - ins

    def _flush_abundances(self, batch) -> int:
        if not batch:
            return 0
        before = self.con.total_changes
        self.con.executemany(
            "INSERT OR IGNORE INTO sample_taxon_abundances (sample_id, taxon_id,"
            " relative_abundance) VALUES (?, ?, ?)",
            batch,
        )
        self.con.commit()
        return self.con.total_changes - before


def main() -> None:
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
    con = sqlite3.connect(DB_PATH)
    con.executescript(open(SCHEMA, encoding="utf-8").read())
    con.execute("PRAGMA foreign_keys = ON")
    loader = Loader(con)

    loader.load_base_taxa()
    print(f"base taxa loaded: {len(loader.taxon_ids)}")

    for path in ASSOC_FILES:
        read, ins, skip = loader.load_associations(path)
        print(f"{os.path.basename(path)}: read={read} inserted={ins} skipped={skip}")
    loader.finalize_comparison_levels()
    print("genus lineage rows backfilled:", loader.backfill_genus_lineage())

    read, ins, skip = loader.load_samples()
    print(f"samples: read={read} inserted={ins} skipped={skip}")

    read, ins, skip = loader.load_abundances()
    print(f"abundances: read={read} inserted={ins} skipped={skip}")

    con.commit()
    con.execute("ANALYZE")
    con.commit()
    for table in ("taxa", "diseases", "studies", "phenotype_comparisons",
                  "taxon_disease_associations", "samples", "sample_taxon_abundances"):
        n = con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        print(f"{table}: {n}")
    con.close()


if __name__ == "__main__":
    main()
