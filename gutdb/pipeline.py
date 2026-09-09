from __future__ import annotations

from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from dataclasses import dataclass
from datetime import datetime
import json
from pathlib import Path
import time
from typing import Any, Iterable, Iterator, Mapping

from mysql.connector import MySQLConnection

from .gmrepo import GMrepoClient
from .cmd_curation import (
    CMD_RAW,
    CMDCurationClient,
    parse_age_years,
    parse_bmi as parse_cmd_bmi,
    parse_sex,
    run_accessions,
)
from .mgnify import (
    MGnifyClient,
    best_analysis_per_sample,
    genus_abundances,
    sample_metadata_dict,
)
from .ncbi import (
    AGE_TAGS,
    NON_PROKARYOTIC_PHYLA,
    TaxonomyClient,
    BMI_TAGS,
    SEX_TAGS,
    SRAClient,
    batched,
    parse_bmi,
    parse_host_age_years,
    parse_host_sex,
)
from .transform import (
    association_direction,
    clean_str,
    normalize_pathogen_flag,
    normalize_phylum,
    normalize_yes_no,
    collapse_energy_modes,
    collapse_food_sources,
    derive_energy_mode,
    derive_primary_food_source,
    extract_species_epithet,
    first_value,
    normalize_genus,
    nullable_float,
    nullable_int,
    parse_scientific_name,
    quantize_effect_size,
)


@dataclass
class LoadStats:
    read: int = 0
    inserted: int = 0
    updated: int = 0
    skipped: int = 0

    def __iadd__(self, other: "LoadStats") -> "LoadStats":
        self.read += other.read
        self.inserted += other.inserted
        self.updated += other.updated
        self.skipped += other.skipped
        return self


def open_dict_rows(path: str | Path) -> Iterator[dict[str, str]]:
    source = Path(path)
    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(8192)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;|")
        except csv.Error:
            dialect = csv.excel_tab if source.suffix.lower() in {".tsv", ".tab"} else csv.excel
        reader = csv.DictReader(handle, dialect=dialect)
        if not reader.fieldnames:
            raise ValueError(f"{source} has no header row")
        yield from reader


def most_common(values: Iterable[Any]) -> str:
    cleaned = [clean_str(v) for v in values]
    cleaned = [v for v in cleaned if v]
    if not cleaned:
        return ""
    counts = Counter(cleaned)
    return min(counts, key=lambda value: (-counts[value], value.casefold()))


def _start_run(connection: MySQLConnection, source_name: str, source_uri: str) -> int:
    cursor = connection.cursor()
    cursor.execute(
        "INSERT INTO ingestion_runs (source_name, source_uri) VALUES (%s, %s)",
        (source_name, source_uri),
    )
    run_id = cursor.lastrowid
    cursor.close()
    connection.commit()
    return run_id


def _finish_run(
    connection: MySQLConnection, run_id: int, stats: LoadStats, error: Exception | None = None
) -> None:
    if error:
        connection.rollback()
    cursor = connection.cursor()
    cursor.execute(
        """
        UPDATE ingestion_runs
        SET completed_at = %s, status = %s, rows_read = %s, rows_inserted = %s,
            rows_updated = %s, rows_skipped = %s, error_message = %s
        WHERE id = %s
        """,
        (
            datetime.now(),
            "failed" if error else "completed",
            stats.read,
            stats.inserted,
            stats.updated,
            stats.skipped,
            str(error)[:65000] if error else None,
            run_id,
        ),
    )
    cursor.close()
    connection.commit()


def _taxon_exists(connection: MySQLConnection, genus_key: str, species_key: str) -> bool:
    cursor = connection.cursor()
    cursor.execute(
        "SELECT 1 FROM taxa WHERE genus_key = %s AND species_key = %s",
        (genus_key, species_key),
    )
    exists = cursor.fetchone() is not None
    cursor.close()
    return exists


NUMERIC_FILL_COLUMNS = frozenset({
    "ncbi_tax_id", "number_of_membranes", "optimal_temperature", "human_pathogen",
})


def normalize_source_database(value: Any) -> str:
    source = clean_str(value).casefold()
    if source == "gmrepo":
        return "GMrepo"
    if source == "pubmed":
        return "PubMed"
    if source == "mgnify":
        return "MGnify"
    return "MiMeDB"


def upsert_taxon(
    connection: MySQLConnection,
    taxon: Mapping[str, Any],
    overwrite: bool = False,
) -> tuple[int, bool]:
    genus = normalize_genus(taxon.get("genus"))
    species = extract_species_epithet(genus, taxon.get("species"))
    if not genus:
        raise ValueError("A taxon must have a genus")
    genus_key, species_key = genus.casefold(), species.casefold()
    existed = _taxon_exists(connection, genus_key, species_key)
    scientific_name = f"{genus} {species}".strip()
    rank = clean_str(taxon.get("taxonomic_rank")) or ("species" if species else "genus")
    if rank not in {"species", "genus", "other"}:
        rank = "other"
    energy = clean_str(taxon.get("energy_mode")) or "N/A"
    if energy not in {"fermenter", "respirator", "mixed", "N/A"}:
        energy = "N/A"
    values = (
        clean_str(taxon.get("superkingdom")) or None,
        clean_str(taxon.get("phylum")) or None,
        clean_str(taxon.get("class_name")) or None,
        clean_str(taxon.get("order_name")) or None,
        clean_str(taxon.get("family")) or None,
        genus,
        species,
        genus_key,
        species_key,
        scientific_name,
        rank,
        nullable_int(taxon.get("ncbi_tax_id")),
        clean_str(taxon.get("gram_stain")) or None,
        clean_str(taxon.get("oxygen_requirement")) or None,
        clean_str(taxon.get("ph_preference")) or None,
        clean_str(taxon.get("sporulation")) or None,
        clean_str(taxon.get("shape")) or None,
        clean_str(taxon.get("cell_arrangement")) or None,
        normalize_yes_no(taxon.get("mobility")) or None,
        normalize_yes_no(taxon.get("flagella_presence")) or None,
        nullable_int(taxon.get("number_of_membranes")),
        clean_str(taxon.get("biotic_relationship")) or None,
        clean_str(taxon.get("habitat")) or None,
        clean_str(taxon.get("temperature_range")) or None,
        nullable_float(taxon.get("optimal_temperature")),
        clean_str(taxon.get("metabolism"))[:255] or None,
        clean_str(taxon.get("energy_source")) or None,
        normalize_pathogen_flag(taxon.get("human_pathogen")),
        energy,
        clean_str(taxon.get("primary_food_source")) or "N/A",
        normalize_source_database(taxon.get("source_database")),
    )
    columns = (
        "superkingdom", "phylum", "class_name", "order_name", "family",
        "genus", "species",
        "genus_key", "species_key", "scientific_name", "taxonomic_rank", "ncbi_tax_id",
        "gram_stain", "oxygen_requirement", "ph_preference", "sporulation",
        "shape", "cell_arrangement", "mobility", "flagella_presence",
        "number_of_membranes", "biotic_relationship", "habitat", "temperature_range",
        "optimal_temperature", "metabolism", "energy_source", "human_pathogen",
        "energy_mode", "primary_food_source", "source_database",
    )
    if overwrite:
        updates = ", ".join(f"{c} = VALUES({c})" for c in columns if c not in {"genus_key", "species_key"})
    else:
        immutable = {"genus_key", "species_key", "genus", "species", "scientific_name", "taxonomic_rank"}
        update_parts = []
        for column in columns:
            if column in immutable:
                continue
            if column in NUMERIC_FILL_COLUMNS:
                # Numeric columns get a NULL-only guard. The generic branch below
                # compares against the string 'N/A', which MySQL coerces to 0 on a
                # numeric column and would therefore clobber a legitimate stored 0.
                update_parts.append(
                    f"{column} = CASE WHEN {column} IS NULL THEN VALUES({column}) ELSE {column} END"
                )
            elif column == "source_database":
                update_parts.append(
                    f"{column} = CASE WHEN VALUES({column}) = 'GMrepo' THEN 'GMrepo' ELSE {column} END"
                )
            else:
                update_parts.append(
                    f"{column} = CASE WHEN {column} IS NULL OR TRIM(CAST({column} AS CHAR)) = '' "
                    f"OR {column} = 'N/A' THEN VALUES({column}) ELSE {column} END"
                )
        updates = ", ".join(update_parts)
    cursor = connection.cursor()
    cursor.execute(
        f"INSERT INTO taxa ({', '.join(columns)}) VALUES ({', '.join(['%s'] * len(columns))}) "
        f"ON DUPLICATE KEY UPDATE {updates}",
        values,
    )
    cursor.execute(
        "SELECT id FROM taxa WHERE genus_key = %s AND species_key = %s",
        (genus_key, species_key),
    )
    taxon_id = cursor.fetchone()[0]
    cursor.close()
    return taxon_id, existed


def _upsert_source_record(
    connection: MySQLConnection,
    taxon_id: int,
    source_database: str,
    source_record_id: str,
    raw_record: Mapping[str, Any],
    source_url: str | None = None,
) -> None:
    """Keep taxon-level source provenance after merging source records into taxa."""
    cursor = connection.cursor()
    cursor.execute(
        """
        UPDATE taxa
        SET source_database = CASE
            WHEN %s = 'GMrepo' THEN 'GMrepo'
            ELSE source_database
        END
        WHERE id = %s
        """,
        (normalize_source_database(source_database), taxon_id),
    )
    cursor.close()


def load_baseline_taxa(
    connection: MySQLConnection, path: str | Path, overwrite: bool = False
) -> LoadStats:
    stats = LoadStats()
    run_id = _start_run(connection, "baseline_csv", str(Path(path).resolve()))
    try:
        for row_number, row in enumerate(open_dict_rows(path), start=2):
            stats.read += 1
            genus = first_value(row, "Genus", "genus")
            species = first_value(row, "Species", "species")
            if not genus:
                stats.skipped += 1
                continue
            taxon = {
                "phylum": first_value(row, "Phylum", "phylum"),
                "class_name": first_value(row, "Class", "class", "klass"),
                "order_name": first_value(row, "Order", "order"),
                "family": first_value(row, "Family", "family"),
                "genus": genus,
                "species": species,
                "oxygen_requirement": first_value(row, "Oxygen", "oxygen_requirement"),
                "ph_preference": first_value(row, "pH", "ph_preference"),
                "sporulation": first_value(row, "Spore", "sporulation"),
                "energy_mode": first_value(row, "energy_mode") or "N/A",
                "primary_food_source": first_value(row, "primary_food_source") or "N/A",
            }
            taxon_id, existed = upsert_taxon(connection, taxon, overwrite=overwrite)
            source_id = first_value(row, "id") or f"row-{row_number}"
            _upsert_source_record(connection, taxon_id, "baseline_csv", source_id, row)
            stats.updated += int(existed)
            stats.inserted += int(not existed)
        _finish_run(connection, run_id, stats)
        return stats
    except Exception as error:
        _finish_run(connection, run_id, stats, error)
        raise


def _collapse_mimedb(path: str | Path) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in open_dict_rows(path):
        genus = normalize_genus(first_value(row, "genus"))
        species = extract_species_epithet(genus, first_value(row, "species", "name"))
        if genus:
            groups[(genus.casefold(), species.casefold())].append(dict(row))

    collapsed: list[dict[str, Any]] = []
    for rows in groups.values():
        genus = normalize_genus(first_value(rows[0], "genus"))
        species = extract_species_epithet(genus, first_value(rows[0], "species", "name"))
        modes = [
            derive_energy_mode(
                first_value(r, "oxygen_requirement"),
                first_value(r, "metabolism"),
                first_value(r, "energy_source"),
                first_value(r, "activity"),
            )
            for r in rows
        ]
        foods = [
            derive_primary_food_source(
                first_value(r, "metabolism"),
                first_value(r, "activity"),
                first_value(r, "energy_source"),
            )
            for r in rows
        ]
        pathogen_flags = [
            normalize_pathogen_flag(first_value(r, "human_pathogen")) for r in rows
        ]
        collapsed.append(
            {
                "superkingdom": most_common(first_value(r, "superkingdom") for r in rows),
                "phylum": most_common(first_value(r, "phylum") for r in rows),
                "class_name": most_common(first_value(r, "klass", "class") for r in rows),
                "order_name": most_common(first_value(r, "order") for r in rows),
                "family": most_common(first_value(r, "family") for r in rows),
                "genus": genus,
                "species": species,
                "ncbi_tax_id": most_common(first_value(r, "ncbi_tax_id") for r in rows),
                "gram_stain": most_common(first_value(r, "gram") for r in rows),
                "oxygen_requirement": most_common(first_value(r, "oxygen_requirement") for r in rows),
                "sporulation": most_common(first_value(r, "sporulation") for r in rows),
                "shape": most_common(first_value(r, "shape") for r in rows),
                "cell_arrangement": most_common(first_value(r, "cell_arrangement") for r in rows),
                "mobility": most_common(normalize_yes_no(first_value(r, "mobility")) for r in rows),
                "flagella_presence": most_common(
                    normalize_yes_no(first_value(r, "flagella_presence")) for r in rows
                ),
                "number_of_membranes": most_common(
                    first_value(r, "number_of_membranes") for r in rows
                ),
                "biotic_relationship": most_common(
                    first_value(r, "biotic_relationship") for r in rows
                ),
                "habitat": most_common(first_value(r, "habitat") for r in rows),
                "temperature_range": most_common(first_value(r, "temperature_range") for r in rows),
                "optimal_temperature": most_common(
                    first_value(r, "optimal_temperature") for r in rows
                ),
                "metabolism": most_common(first_value(r, "metabolism") for r in rows),
                "energy_source": most_common(first_value(r, "energy_source") for r in rows),
                # any strain documented as a pathogen makes the species positive;
                # unannotated strains stay silent rather than voting "not a pathogen"
                "human_pathogen": 1 if any(f == 1 for f in pathogen_flags) else None,
                "energy_mode": collapse_energy_modes(modes),
                "primary_food_source": collapse_food_sources(foods),
                "source_rows": rows,
            }
        )
    return collapsed


def load_mimedb(
    connection: MySQLConnection, path: str | Path, overwrite: bool = False
) -> LoadStats:
    stats = LoadStats()
    run_id = _start_run(connection, "MiMeDB", str(Path(path).resolve()))
    try:
        for taxon in _collapse_mimedb(path):
            stats.read += len(taxon["source_rows"])
            taxon_id, existed = upsert_taxon(connection, taxon, overwrite=overwrite)
            for raw in taxon.pop("source_rows"):
                source_id = first_value(raw, "microbe_id", "id", "name")
                _upsert_source_record(
                    connection,
                    taxon_id,
                    "MiMeDB",
                    source_id,
                    raw,
                    f"https://mimedb.org/microbes/{source_id}" if source_id else None,
                )
            stats.updated += int(existed)
            stats.inserted += int(not existed)
        _finish_run(connection, run_id, stats)
        return stats
    except Exception as error:
        _finish_run(connection, run_id, stats, error)
        raise


def upsert_disease(
    connection: MySQLConnection,
    name: str,
    mesh_id: str = "",
) -> int:
    name = clean_str(name)
    if not name:
        raise ValueError("Disease/phenotype name is required")
    name_key = name.casefold()
    cursor = connection.cursor()
    cursor.execute(
        """
        INSERT INTO diseases (mesh_id, name, name_key)
        VALUES (%s, %s, %s)
        ON DUPLICATE KEY UPDATE
            mesh_id = COALESCE(VALUES(mesh_id), mesh_id)
        """,
        (clean_str(mesh_id) or None, name, name_key),
    )
    cursor.execute(
        "SELECT id FROM diseases WHERE name_key = %s OR (mesh_id IS NOT NULL AND mesh_id = %s) LIMIT 1",
        (name_key, clean_str(mesh_id) or None),
    )
    disease_id = cursor.fetchone()[0]
    cursor.close()
    return disease_id


def upsert_study(connection: MySQLConnection, row: Mapping[str, Any]) -> int:
    project = first_value(row, "project_id", "project_accession", "project", "study_id")
    if not project:
        raise ValueError("project_id/project_accession is required")
    cursor = connection.cursor()
    cursor.execute(
        """
        INSERT INTO studies
            (project_accession, title, description, data_type, source_database,
             data_quality)
        VALUES (%s, %s, %s, %s, %s, %s)
        ON DUPLICATE KEY UPDATE
            title = COALESCE(NULLIF(VALUES(title), ''), title),
            description = COALESCE(NULLIF(VALUES(description), ''), description),
            data_type = COALESCE(NULLIF(VALUES(data_type), ''), data_type),
            data_quality = VALUES(data_quality)
        """,
        (
            project,
            first_value(row, "project_title", "title") or None,
            first_value(row, "project_description", "description") or None,
            first_value(row, "data_type", "experiment_type") or None,
            first_value(row, "source_database") or "GMrepo",
            first_value(row, "data_quality") if first_value(row, "data_quality") in {"curated", "qualified"} else "unknown",
        ),
    )
    cursor.execute("SELECT id FROM studies WHERE project_accession = %s", (project,))
    study_id = cursor.fetchone()[0]
    cursor.close()
    return study_id


def _upsert_comparison(
    connection: MySQLConnection,
    study_id: int,
    project_id: str,
    phenotype_a_id: int,
    phenotype_a: str,
    phenotype_b_id: int,
    phenotype_b: str,
    positive_id: int | None,
    negative_id: int | None,
    method: str,
    rank: str,
    notes: str = "",
) -> int:
    comparison_key = "|".join(
        [project_id.casefold(), phenotype_a.casefold(), phenotype_b.casefold(), method.casefold(), rank]
    )
    cursor = connection.cursor()
    cursor.execute(
        """
        INSERT INTO phenotype_comparisons
            (study_id, phenotype_a_id, phenotype_b_id, comparison_key, method,
             taxonomic_level, positive_score_enriched_in_id, negative_score_enriched_in_id, notes)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
        ON DUPLICATE KEY UPDATE
            positive_score_enriched_in_id = VALUES(positive_score_enriched_in_id),
            negative_score_enriched_in_id = VALUES(negative_score_enriched_in_id),
            notes = COALESCE(notes, VALUES(notes))
        """,
        (
            study_id, phenotype_a_id, phenotype_b_id, comparison_key, method, rank,
            positive_id, negative_id, notes or None,
        ),
    )
    cursor.execute("SELECT id FROM phenotype_comparisons WHERE comparison_key = %s", (comparison_key,))
    comparison_id = cursor.fetchone()[0]
    cursor.close()
    return comparison_id


def load_associations(
    connection: MySQLConnection, path: str | Path, overwrite: bool = True
) -> LoadStats:
    stats = LoadStats()
    run_id = _start_run(connection, "GMrepo_associations", str(Path(path).resolve()))
    try:
        for row in open_dict_rows(path):
            stats.read += 1
            scientific_name = first_value(row, "scientific_name", "marker_taxon", "taxon", "microbe")
            row_rank = first_value(row, "taxonomic_rank", "rank").casefold() or "species"
            genus = first_value(row, "genus")
            species = first_value(row, "species", "species_epithet")
            if scientific_name and not genus:
                if row_rank == "genus":
                    genus, species = normalize_genus(scientific_name), ""
                else:
                    genus, species = parse_scientific_name(scientific_name)
            if not genus:
                stats.skipped += 1
                continue

            disease_name = first_value(row, "disease_name", "disease", "phenotype")
            phenotype_a = first_value(row, "phenotype_a", "phenotype1") or "Health"
            phenotype_b = first_value(row, "phenotype_b", "phenotype2") or disease_name
            if not disease_name:
                disease_name = phenotype_b if phenotype_b.casefold() != "health" else phenotype_a
            if not disease_name or not phenotype_a or not phenotype_b:
                stats.skipped += 1
                continue

            disease_id = upsert_disease(
                connection,
                disease_name,
                first_value(row, "mesh_id", "disease_mesh_id"),
            )
            # The literature marker CSVs no longer carry the phenotype columns
            # (phenotype_a was the constant "Health"/D006262 and phenotype_b
            # restated disease_name/mesh_id); the defaults on lines 556-557 and
            # here reconstruct them. The GMrepo export still supplies its own.
            a_mesh = first_value(row, "phenotype_a_mesh_id") or (
                "D006262" if phenotype_a == "Health" else None)
            a_id = upsert_disease(connection, phenotype_a, a_mesh)
            b_id = upsert_disease(connection, phenotype_b, first_value(row, "phenotype_b_mesh_id", "mesh_id"))
            positive_name = first_value(row, "positive_enriched_in")
            negative_name = first_value(row, "negative_enriched_in")
            positive_id = upsert_disease(connection, positive_name) if positive_name else None
            negative_id = upsert_disease(connection, negative_name) if negative_name else None
            study_id = upsert_study(connection, row)
            project_id = first_value(row, "project_id", "project_accession", "project", "study_id")
            rank = row_rank
            if rank not in {"species", "genus", "mixed"}:
                rank = "species"
            method = first_value(row, "method", "effect_type") or "LEfSe"
            comparison_id = _upsert_comparison(
                connection, study_id, project_id, a_id, phenotype_a, b_id, phenotype_b,
                positive_id, negative_id, method, rank,
                first_value(row, "notes", "comparison_notes"),
            )
            taxon_id, taxon_existed = upsert_taxon(
                connection,
                {
                    "genus": genus,
                    "species": species,
                    "taxonomic_rank": "genus" if rank == "genus" else "species",
                    "ncbi_tax_id": first_value(row, "ncbi_tax_id", "taxon_id"),
                    "source_database": first_value(row, "source_database") or "GMrepo",
                },
                overwrite=False,
            )
            lda_raw = nullable_float(first_value(row, "lda_score", "effect_size", "score"))
            # Direction is derived from the FULL-PRECISION value: quantizing first
            # could send a tiny positive score to 0.00, which association_direction
            # reads as "marker" rather than "enriched".
            direction = association_direction(
                lda_raw, disease_name, positive_name, negative_name,
                first_value(row, "direction"),
            )
            lda = quantize_effect_size(lda_raw)
            cursor = connection.cursor()
            cursor.execute(
                """
                SELECT 1 FROM taxon_disease_associations
                WHERE taxon_id = %s AND disease_id = %s AND comparison_id = %s
                    AND effect_type = %s
                """,
                (taxon_id, disease_id, comparison_id, first_value(row, "effect_type") or "LDA"),
            )
            existed = cursor.fetchone() is not None
            cursor.execute(
                """
                INSERT INTO taxon_disease_associations
                    (taxon_id, disease_id, comparison_id, direction, effect_type, effect_size,
                     source_database)
                VALUES (%s, %s, %s, %s, %s, %s, %s)
                ON DUPLICATE KEY UPDATE
                    direction = VALUES(direction),
                    effect_size = IF(%s, VALUES(effect_size), COALESCE(effect_size, VALUES(effect_size))),
                    source_database = VALUES(source_database)
                """,
                (
                    taxon_id, disease_id, comparison_id, direction,
                    first_value(row, "effect_type") or "LDA", lda,
                    first_value(row, "source_database") or "GMrepo",
                    overwrite,
                ),
            )
            cursor.close()
            stats.updated += int(existed)
            stats.inserted += int(not existed)
            if not taxon_existed:
                _upsert_source_record(
                    connection, taxon_id, first_value(row, "source_database") or "GMrepo",
                    f"{project_id}:{scientific_name or f'{genus} {species}'}", row,
                    first_value(row, "association_url", "source_url") or None,
                )
        _finish_run(connection, run_id, stats)
        return stats
    except Exception as error:
        _finish_run(connection, run_id, stats, error)
        raise


GMREPO_ASSOCIATION_FIELDS = [
    "project_id", "project_title", "project_description", "data_type", "data_quality",
    "scientific_name", "taxonomic_rank", "ncbi_tax_id", "disease_name", "mesh_id",
    "phenotype_a", "phenotype_a_mesh_id", "phenotype_b", "phenotype_b_mesh_id",
    "negative_enriched_in", "positive_enriched_in", "method", "effect_type", "lda_score",
    "direction", "association_url", "source_database",
]


def sync_gmrepo_comparison(
    connection: MySQLConnection,
    client: GMrepoClient,
    mesh_id1: str,
    mesh_id2: str,
    output_path: str | Path,
    exclude_projects: Iterable[str] = (),
    overwrite: bool = True,
) -> tuple[LoadStats, list[str]]:
    """Fetch a curated cross-project marker comparison, save it, and load it."""
    payload = client.phenotype_comparison(mesh_id1, mesh_id2)
    stats = payload.get("stats", {})
    phenotype1 = clean_str(stats.get("phenotype1_term"))
    phenotype2 = clean_str(stats.get("phenotype2_term"))
    if not phenotype1 or not phenotype2:
        raise ValueError("GMrepo comparison response did not include phenotype names")
    excluded = {clean_str(project).casefold() for project in exclude_projects}
    rows = []
    project_ids = set()
    for marker in payload.get("alldata", []):
        project_id = clean_str(marker.get("project_id"))
        if not project_id or project_id.casefold() in excluded:
            continue
        project_ids.add(project_id)
        rows.append(
            {
                "project_id": project_id,
                "project_title": f"GMrepo curated {phenotype1} vs. {phenotype2}",
                "project_description": "Curated project-level marker taxa from GMrepo",
                "data_type": clean_str(marker.get("experiment_type")),
                "data_quality": "curated",
                "scientific_name": clean_str(marker.get("scientific_name")),
                "taxonomic_rank": clean_str(marker.get("taxon_rank_level")),
                "ncbi_tax_id": clean_str(marker.get("ncbi_taxon_id")),
                "disease_name": phenotype2,
                "mesh_id": mesh_id2,
                "phenotype_a": phenotype1,
                "phenotype_a_mesh_id": mesh_id1,
                "phenotype_b": phenotype2,
                "phenotype_b_mesh_id": mesh_id2,
                "negative_enriched_in": phenotype1,
                "positive_enriched_in": phenotype2,
                "method": "LEfSe",
                "effect_type": "LDA",
                "lda_score": marker.get("LDA"),
                "direction": "",
                "association_url": (
                    f"{client.base_url}/data/project/{project_id}/{mesh_id1}/{mesh_id2}"
                ),
                "source_database": "GMrepo",
            }
        )
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GMREPO_ASSOCIATION_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    return load_associations(connection, destination, overwrite), sorted(project_ids)


def _comparison_rows(
    client: GMrepoClient,
    mesh_id1: str,
    mesh_id2: str,
    include_health: bool = False,
) -> tuple[list[dict[str, Any]], set[str]]:
    payload = client.phenotype_comparison(mesh_id1, mesh_id2)
    stats = payload.get("stats", {})
    phenotype1 = clean_str(stats.get("phenotype1_term"))
    phenotype2 = clean_str(stats.get("phenotype2_term"))
    if not phenotype1 or not phenotype2:
        raise ValueError(f"Missing phenotype names for {mesh_id1} vs {mesh_id2}")
    targets = [(phenotype1, mesh_id1, -1), (phenotype2, mesh_id2, 1)]
    if not include_health:
        targets = [target for target in targets if target[0].casefold() != "health"]
    rows: list[dict[str, Any]] = []
    projects: set[str] = set()
    for marker in payload.get("alldata", []):
        project_id = clean_str(marker.get("project_id"))
        scientific_name = clean_str(marker.get("scientific_name"))
        lda = nullable_float(marker.get("LDA"))
        if not project_id or not scientific_name or lda is None:
            continue
        projects.add(project_id)
        for disease_name, disease_mesh_id, enriched_sign in targets:
            direction = "enriched" if (lda > 0) == (enriched_sign > 0) else "depleted"
            rows.append(
                {
                    "project_id": project_id,
                    "project_title": f"GMrepo curated {phenotype1} vs. {phenotype2}",
                    "project_description": "Curated project-level marker taxa from GMrepo",
                    "data_type": clean_str(marker.get("experiment_type")),
                    "data_quality": "curated",
                    "scientific_name": scientific_name,
                    "taxonomic_rank": clean_str(marker.get("taxon_rank_level")),
                    "ncbi_tax_id": clean_str(marker.get("ncbi_taxon_id")),
                    "disease_name": disease_name,
                    "mesh_id": disease_mesh_id,
                    "phenotype_a": phenotype1,
                    "phenotype_a_mesh_id": mesh_id1,
                    "phenotype_b": phenotype2,
                    "phenotype_b_mesh_id": mesh_id2,
                    "negative_enriched_in": phenotype1,
                    "positive_enriched_in": phenotype2,
                    "method": "LEfSe",
                    "effect_type": "LDA",
                    "lda_score": lda,
                    "direction": direction,
                    "association_url": (
                        f"{client.base_url}/data/project/{project_id}/{mesh_id1}/{mesh_id2}"
                    ),
                    "source_database": "GMrepo",
                }
            )
    return rows, projects


def sync_all_gmrepo_comparisons(
    connection: MySQLConnection,
    client: GMrepoClient,
    output_path: str | Path,
    workers: int = 6,
    include_health: bool = False,
    overwrite: bool = True,
) -> tuple[LoadStats, dict[str, int]]:
    """Traverse every current curated phenotype comparison and load project evidence."""
    catalog = client.all_phenotype_comparisons()
    comparisons = [
        (clean_str(row.get("phenotype1")), clean_str(row.get("phenotype2")))
        for row in catalog.get("data", [])
        if clean_str(row.get("phenotype1")) and clean_str(row.get("phenotype2"))
    ]
    all_rows: list[dict[str, Any]] = []
    all_projects: set[str] = set()
    failures: list[str] = []

    def fetch(pair: tuple[str, str]) -> tuple[list[dict[str, Any]], set[str]]:
        local_client = GMrepoClient(client.base_url, timeout=client.timeout)
        return _comparison_rows(local_client, pair[0], pair[1], include_health)

    with ThreadPoolExecutor(max_workers=max(1, workers)) as executor:
        futures = {executor.submit(fetch, pair): pair for pair in comparisons}
        for future in as_completed(futures):
            pair = futures[future]
            try:
                rows, projects = future.result()
                all_rows.extend(rows)
                all_projects.update(projects)
            except Exception as error:
                failures.append(f"{pair[0]} vs {pair[1]}: {error}")
    if failures:
        raise RuntimeError(
            f"Failed to fetch {len(failures)} of {len(comparisons)} comparisons: "
            + "; ".join(failures[:10])
        )

    all_rows.sort(
        key=lambda row: (
            row["disease_name"].casefold(), row["project_id"],
            row["scientific_name"].casefold(), row["direction"],
        )
    )
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GMREPO_ASSOCIATION_FIELDS)
        writer.writeheader()
        writer.writerows(all_rows)
    loaded = load_associations(connection, destination, overwrite)
    pruned = prune_gmrepo_associations(connection, destination)
    metadata = {
        "comparisons": len(comparisons),
        "projects": len(all_projects),
        "exported_rows": len(all_rows),
        "pruned_rows": pruned,
    }
    return loaded, metadata


GMREPO_SAMPLE_FIELDS = [
    "project_id", "data_type", "data_quality", "run_id", "sample_id", "disease_name",
    "mesh_id", "sex", "age_years", "country", "qc_status", "instrument_model",
    "nr_reads_sequenced", "longitude", "latitude", "phenotypes", "nr_phenotypes",
]


def _parse_gmrepo_phenotypes(value: Any) -> list[dict[str, Any]]:
    if isinstance(value, list):
        return [row for row in value if isinstance(row, dict)]
    text = clean_str(value)
    if not text:
        return []
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        return []
    return [row for row in parsed if isinstance(row, dict)] if isinstance(parsed, list) else []


def _choose_sample_phenotype(value: Any) -> tuple[str, str]:
    phenotypes = _parse_gmrepo_phenotypes(value)
    if not phenotypes:
        return "", ""
    chosen = next(
        (
            row for row in phenotypes
            if clean_str(row.get("term")).casefold() != "health"
        ),
        phenotypes[0],
    )
    return clean_str(chosen.get("term")), clean_str(chosen.get("disease"))


def _study_projects(connection: MySQLConnection) -> list[str]:
    cursor = connection.cursor()
    cursor.execute("SELECT project_accession FROM studies ORDER BY project_accession")
    projects = [row[0] for row in cursor.fetchall()]
    cursor.close()
    return projects


def sync_gmrepo_samples(
    connection: MySQLConnection,
    client: GMrepoClient,
    output_path: str | Path,
    limit: int = 1000,
) -> tuple[LoadStats, dict[str, int]]:
    """Fetch GMrepo run/sample metadata for every imported study project."""
    projects = _study_projects(connection)
    rows: list[dict[str, Any]] = []
    for project_id in projects:
        skip = 0
        while True:
            page = client.project_runs(project_id, limit=limit, skip=skip)
            if not page:
                break
            for raw in page:
                disease_name, mesh_id = _choose_sample_phenotype(raw.get("phenotypes"))
                rows.append(
                    {
                        "project_id": clean_str(raw.get("project_id")) or project_id,
                        "data_type": clean_str(raw.get("experiment_type")),
                        "data_quality": "curated",
                        "run_id": clean_str(raw.get("run_id")),
                        "sample_id": clean_str(raw.get("sample_id")),
                        "disease_name": disease_name,
                        "mesh_id": mesh_id,
                        "sex": clean_str(raw.get("sex")),
                        "age_years": raw.get("host_age"),
                        "country": clean_str(raw.get("country")),
                        "qc_status": raw.get("QCStatus"),
                        "instrument_model": clean_str(raw.get("instrument_model")),
                        "nr_reads_sequenced": raw.get("nr_reads_sequenced"),
                        "longitude": raw.get("longitude"),
                        "latitude": raw.get("latitude"),
                        "phenotypes": raw.get("phenotypes"),
                        "nr_phenotypes": raw.get("nr_phenotypes"),
                    }
                )
            if len(page) < limit:
                break
            skip += limit

    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GMREPO_SAMPLE_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    loaded = load_samples(connection, destination)
    metadata = {
        "projects": len(projects),
        "exported_rows": len(rows),
    }
    return loaded, metadata


def prune_gmrepo_associations(connection: MySQLConnection, current_export: str | Path) -> int:
    """Remove GMrepo evidence rows absent from the current authoritative export."""
    current_keys = set()
    for row in open_dict_rows(current_export):
        rank = first_value(row, "taxonomic_rank", "rank").casefold() or "species"
        scientific_name = first_value(row, "scientific_name", "marker_taxon", "taxon")
        if rank == "genus":
            genus, species = normalize_genus(scientific_name), ""
        else:
            genus, species = parse_scientific_name(scientific_name)
        current_keys.add(
            (
                first_value(row, "project_id", "project_accession").casefold(),
                first_value(row, "phenotype_a", "phenotype1").casefold(),
                first_value(row, "phenotype_b", "phenotype2").casefold(),
                rank,
                (first_value(row, "method") or "LEfSe").casefold(),
                first_value(row, "disease_name", "disease").casefold(),
                normalize_genus(genus).casefold(),
                extract_species_epithet(genus, species).casefold(),
                (first_value(row, "effect_type") or "LDA").casefold(),
            )
        )
    cursor = connection.cursor(dictionary=True)
    cursor.execute(
        """
        SELECT a.id, s.project_accession, da.name AS phenotype_a, db.name AS phenotype_b,
               pc.taxonomic_level, pc.method, d.name AS disease_name, t.genus_key, t.species_key,
               a.effect_type
        FROM gut_microbiome.taxon_disease_associations a
        JOIN gut_microbiome.taxa t ON t.id = a.taxon_id
        JOIN gut_microbiome.diseases d ON d.id = a.disease_id
        JOIN gut_microbiome.phenotype_comparisons pc ON pc.id = a.comparison_id
        JOIN gut_microbiome.diseases da ON da.id = pc.phenotype_a_id
        JOIN gut_microbiome.diseases db ON db.id = pc.phenotype_b_id
        JOIN gut_microbiome.studies s ON s.id = pc.study_id
        WHERE a.source_database = 'GMrepo'
        """
    )
    stale_ids = []
    for row in cursor.fetchall():
        key = (
            clean_str(row["project_accession"]).casefold(),
            clean_str(row["phenotype_a"]).casefold(),
            clean_str(row["phenotype_b"]).casefold(),
            clean_str(row["taxonomic_level"]).casefold(),
            clean_str(row["method"]).casefold(),
            clean_str(row["disease_name"]).casefold(),
            clean_str(row["genus_key"]).casefold(),
            clean_str(row["species_key"]).casefold(),
            clean_str(row["effect_type"]).casefold(),
        )
        if key not in current_keys:
            stale_ids.append((row["id"],))
    if stale_ids:
        cursor.executemany(
            "DELETE FROM gut_microbiome.taxon_disease_associations WHERE id = %s",
            stale_ids,
        )
        cursor.execute(
            """
            DELETE pc FROM gut_microbiome.phenotype_comparisons pc
            LEFT JOIN gut_microbiome.taxon_disease_associations a ON a.comparison_id = pc.id
            WHERE a.id IS NULL
            """
        )
    cursor.close()
    connection.commit()
    return len(stale_ids)


def sync_legacy_microbiome_tables(connection: MySQLConnection) -> dict[str, int]:
    """Mirror normalized taxa and aggregated disease evidence into the legacy schema."""
    cursor = connection.cursor()
    cursor.execute(
        """
        SELECT COUNT(*) FROM information_schema.TABLES
        WHERE TABLE_SCHEMA = 'microbiome_dataset' AND TABLE_NAME = 'gut_microbiome_new'
        """
    )
    if cursor.fetchone()[0] != 1:
        cursor.close()
        raise RuntimeError("microbiome_dataset.gut_microbiome_new is missing")

    # The child table uses an empty species epithet for genus-level taxa. Normalize
    # older NULL parent keys so the composite foreign key remains satisfiable.
    cursor.execute(
        "UPDATE microbiome_dataset.gut_microbiome_new SET Species = '' WHERE Species IS NULL"
    )

    cursor.execute("SELECT COUNT(*) FROM microbiome_dataset.gut_microbiome_new")
    taxa_before = int(cursor.fetchone()[0])
    cursor.execute(
        """
        INSERT INTO microbiome_dataset.gut_microbiome_new
            (Phylum, Class, `Order`, Family, Genus, Species, Oxygen, pH, Spore,
             energy_mode, primary_food_source)
        SELECT
            t.phylum, t.class_name, t.order_name, t.family, t.genus,
            CASE WHEN t.species = '' THEN ''
                 ELSE CONCAT(UPPER(LEFT(t.species, 1)), SUBSTRING(t.species, 2)) END,
            t.oxygen_requirement, t.ph_preference, t.sporulation,
            t.energy_mode, t.primary_food_source
        FROM gut_microbiome.taxa t
        LEFT JOIN microbiome_dataset.gut_microbiome_new legacy
          ON LOWER(TRIM(legacy.Genus)) = t.genus_key
         AND LOWER(TRIM(COALESCE(legacy.Species, ''))) = t.species_key
        WHERE legacy.id IS NULL
          AND CHAR_LENGTH(t.genus) <= 100
          AND CHAR_LENGTH(t.species) <= 150
        ON DUPLICATE KEY UPDATE
            Phylum = COALESCE(NULLIF(gut_microbiome_new.Phylum, ''), VALUES(Phylum)),
            Class = COALESCE(NULLIF(gut_microbiome_new.Class, ''), VALUES(Class)),
            `Order` = COALESCE(NULLIF(gut_microbiome_new.`Order`, ''), VALUES(`Order`)),
            Family = COALESCE(NULLIF(gut_microbiome_new.Family, ''), VALUES(Family)),
            Oxygen = COALESCE(NULLIF(gut_microbiome_new.Oxygen, ''), VALUES(Oxygen)),
            pH = COALESCE(NULLIF(gut_microbiome_new.pH, ''), VALUES(pH)),
            Spore = COALESCE(NULLIF(gut_microbiome_new.Spore, ''), VALUES(Spore)),
            energy_mode = COALESCE(NULLIF(gut_microbiome_new.energy_mode, ''), VALUES(energy_mode)),
            primary_food_source = COALESCE(
                NULLIF(gut_microbiome_new.primary_food_source, ''), VALUES(primary_food_source)
            )
        """
    )
    cursor.execute("SELECT COUNT(*) FROM microbiome_dataset.gut_microbiome_new")
    taxa_after = int(cursor.fetchone()[0])

    cursor.execute(
        """
        DELETE md FROM microbiome_dataset.microbe_disease md
        LEFT JOIN (
            SELECT DISTINCT t.genus_key, t.species_key, d.mesh_id
            FROM gut_microbiome.taxon_disease_associations a
            JOIN gut_microbiome.taxa t ON t.id = a.taxon_id
            JOIN gut_microbiome.diseases d ON d.id = a.disease_id
            WHERE d.mesh_id IS NOT NULL AND d.name <> 'Health'
        ) active
          ON LOWER(TRIM(md.genus)) = active.genus_key
         AND LOWER(TRIM(COALESCE(md.species, ''))) = active.species_key
         AND md.disease_mesh_id = active.mesh_id
        WHERE active.mesh_id IS NULL
        """
    )

    cursor.execute(
        """
        INSERT INTO microbiome_dataset.microbe_disease
            (genus, species, disease_mesh_id, disease_name, ncbi_taxon_id, effect,
             nr_studies_enriched, nr_studies_depleted, mean_rel_abundance, date_retrieved)
        SELECT
            legacy.Genus,
            COALESCE(legacy.Species, ''),
            d.mesh_id,
            d.name,
            MAX(t.ncbi_tax_id),
            CASE
                WHEN COUNT(DISTINCT CASE WHEN a.direction = 'enriched' THEN s.id END) >
                     COUNT(DISTINCT CASE WHEN a.direction = 'depleted' THEN s.id END)
                    THEN 'enriched'
                WHEN COUNT(DISTINCT CASE WHEN a.direction = 'depleted' THEN s.id END) >
                     COUNT(DISTINCT CASE WHEN a.direction = 'enriched' THEN s.id END)
                    THEN 'depleted'
                ELSE 'unclear'
            END,
            COUNT(DISTINCT CASE WHEN a.direction = 'enriched' THEN s.id END),
            COUNT(DISTINCT CASE WHEN a.direction = 'depleted' THEN s.id END),
            NULL,
            CURRENT_DATE
        FROM gut_microbiome.taxon_disease_associations a
        JOIN gut_microbiome.taxa t ON t.id = a.taxon_id
        JOIN gut_microbiome.diseases d ON d.id = a.disease_id
        JOIN gut_microbiome.phenotype_comparisons pc ON pc.id = a.comparison_id
        JOIN gut_microbiome.studies s ON s.id = pc.study_id
        JOIN microbiome_dataset.gut_microbiome_new legacy
          ON LOWER(TRIM(legacy.Genus)) = t.genus_key
         AND LOWER(TRIM(COALESCE(legacy.Species, ''))) = t.species_key
        WHERE d.mesh_id IS NOT NULL AND d.name <> 'Health'
        GROUP BY legacy.Genus, COALESCE(legacy.Species, ''), d.mesh_id, d.name
        ON DUPLICATE KEY UPDATE
            disease_name = VALUES(disease_name),
            ncbi_taxon_id = COALESCE(VALUES(ncbi_taxon_id), ncbi_taxon_id),
            effect = VALUES(effect),
            nr_studies_enriched = VALUES(nr_studies_enriched),
            nr_studies_depleted = VALUES(nr_studies_depleted),
            date_retrieved = VALUES(date_retrieved)
        """
    )
    cursor.execute("SELECT COUNT(*) FROM microbiome_dataset.microbe_disease")
    disease_rows = int(cursor.fetchone()[0])
    cursor.execute(
        """
        CREATE OR REPLACE VIEW microbiome_dataset.gut_microbiome_new_full AS
        SELECT Phylum, Class, `Order`, Family, Genus, Species, Oxygen, pH, Spore,
               energy_mode, primary_food_source, id
        FROM microbiome_dataset.gut_microbiome_new
        """
    )
    cursor.close()
    connection.commit()
    return {
        "legacy_taxa_before": taxa_before,
        "legacy_taxa_after": taxa_after,
        "legacy_taxa_inserted": taxa_after - taxa_before,
        "legacy_associations": disease_rows,
    }


def sync_gmrepo_phenotypes(connection: MySQLConnection, client: GMrepoClient) -> LoadStats:
    stats = LoadStats()
    run_id = _start_run(connection, "GMrepo_API_phenotypes", f"{client.base_url}/api/get_all_phenotypes")
    try:
        for row in client.all_phenotypes():
            stats.read += 1
            name = first_value(row, "disease_name", "phenotype", "name", "mesh_name")
            mesh_id = first_value(row, "mesh_id", "meshid", "mesh")
            if not name:
                stats.skipped += 1
                continue
            cursor = connection.cursor()
            cursor.execute(
                "SELECT 1 FROM diseases WHERE name_key = %s OR mesh_id = %s LIMIT 1",
                (name.casefold(), mesh_id or None),
            )
            existed = cursor.fetchone() is not None
            cursor.close()
            upsert_disease(connection, name, mesh_id)
            stats.updated += int(existed)
            stats.inserted += int(not existed)
        _finish_run(connection, run_id, stats)
        return stats
    except Exception as error:
        _finish_run(connection, run_id, stats, error)
        raise


def load_samples(connection: MySQLConnection, path: str | Path) -> LoadStats:
    """Load normalized clinical/run metadata exported from GMrepo or assembled locally."""
    stats = LoadStats()
    run_id = _start_run(connection, "GMrepo_samples", str(Path(path).resolve()))
    try:
        for row in open_dict_rows(path):
            stats.read += 1
            project = first_value(row, "project_id", "project_accession", "project")
            run_accession = first_value(row, "run_id", "run_accession", "run")
            if not project or not run_accession:
                stats.skipped += 1
                continue
            study_id = upsert_study(connection, row)
            disease_name = first_value(row, "disease_name", "disease", "phenotype")
            disease_id = (
                upsert_disease(connection, disease_name, first_value(row, "mesh_id"))
                if disease_name
                else None
            )
            cursor = connection.cursor()
            cursor.execute("SELECT 1 FROM samples WHERE run_accession = %s", (run_accession,))
            existed = cursor.fetchone() is not None
            cursor.execute(
                """
                INSERT INTO samples
                    (study_id, disease_id, gmrepo_sample_id, run_accession,
                     sex, age_years, bmi, country, qc_status, body_site)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                ON DUPLICATE KEY UPDATE
                    study_id = VALUES(study_id), disease_id = COALESCE(VALUES(disease_id), disease_id),
                    gmrepo_sample_id = COALESCE(NULLIF(VALUES(gmrepo_sample_id), ''), gmrepo_sample_id),
                    sex = COALESCE(NULLIF(VALUES(sex), ''), sex),
                    age_years = COALESCE(VALUES(age_years), age_years),
                    bmi = COALESCE(VALUES(bmi), bmi),
                    country = COALESCE(NULLIF(VALUES(country), ''), country),
                    qc_status = COALESCE(NULLIF(VALUES(qc_status), ''), qc_status),
                    body_site = COALESCE(NULLIF(VALUES(body_site), ''), body_site)
                """,
                (
                    study_id, disease_id,
                    first_value(row, "sample_id", "gmrepo_sample_id") or None,
                    run_accession,
                    first_value(row, "sex", "gender") or None,
                    nullable_float(first_value(row, "age_years", "age")),
                    nullable_float(first_value(row, "bmi")),
                    first_value(row, "country", "population") or None,
                    first_value(row, "qc_status", "QCStatus", "quality") or None,
                    first_value(row, "body_site", "body site", "sample_type") or None,
                ),
            )
            cursor.close()
            stats.updated += int(existed)
            stats.inserted += int(not existed)
        _finish_run(connection, run_id, stats)
        return stats
    except Exception as error:
        _finish_run(connection, run_id, stats, error)
        raise


def load_abundances(connection: MySQLConnection, path: str | Path) -> LoadStats:
    """Load long-format rows: run_id, scientific_name, relative_abundance."""
    stats = LoadStats()
    run_id = _start_run(connection, "GMrepo_abundances", str(Path(path).resolve()))
    try:
        for row in open_dict_rows(path):
            stats.read += 1
            run_accession = first_value(row, "run_id", "run_accession", "run")
            scientific_name = first_value(row, "scientific_name", "taxon", "microbe")
            abundance = nullable_float(first_value(row, "relative_abundance", "abundance"))
            if not run_accession or not scientific_name or abundance is None:
                stats.skipped += 1
                continue
            genus, species = parse_scientific_name(scientific_name)
            taxon_id, _ = upsert_taxon(
                connection,
                {
                    "genus": genus,
                    "species": species,
                    "taxonomic_rank": first_value(row, "taxonomic_rank", "rank") or "species",
                    "ncbi_tax_id": first_value(row, "ncbi_tax_id", "taxon_id"),
                    "source_database": first_value(row, "source_database") or "GMrepo",
                },
            )
            cursor = connection.cursor()
            cursor.execute("SELECT id FROM samples WHERE run_accession = %s", (run_accession,))
            sample = cursor.fetchone()
            if not sample:
                cursor.close()
                stats.skipped += 1
                continue
            cursor.execute(
                "SELECT 1 FROM sample_taxon_abundances WHERE sample_id = %s AND taxon_id = %s",
                (sample[0], taxon_id),
            )
            existed = cursor.fetchone() is not None
            cursor.execute(
                """
                INSERT INTO sample_taxon_abundances
                    (sample_id, taxon_id, relative_abundance, detection_threshold)
                VALUES (%s, %s, %s, %s)
                ON DUPLICATE KEY UPDATE relative_abundance = VALUES(relative_abundance),
                    detection_threshold = VALUES(detection_threshold)
                """,
                (
                    sample[0], taxon_id, abundance,
                    nullable_float(first_value(row, "detection_threshold")),
                ),
            )
            cursor.close()
            stats.updated += int(existed)
            stats.inserted += int(not existed)
        _finish_run(connection, run_id, stats)
        return stats
    except Exception as error:
        _finish_run(connection, run_id, stats, error)
        raise


def _write_demographics_audit(path: str | Path, rows: list[dict[str, Any]]) -> None:
    fieldnames = [
        "run_accession", "sex", "age_years", "age_note", "bmi", "has_biosample_attrs",
    ]
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def enrich_sample_demographics(
    connection: MySQLConnection,
    client: SRAClient,
    limit: int | None = None,
    batch_size: int = 100,
    sleep_seconds: float = 0.4,
    audit_csv: str | Path | None = None,
) -> tuple[LoadStats, dict[str, int]]:
    """Backfill sex/age_years/bmi from NCBI BioSample attributes embedded in SRA records.

    Only fills currently NULL fields via COALESCE; never overwrites a curated
    GMrepo value. host_age is unit-aware (days/weeks/months/years) and rejected
    outright if it converts to an implausible age, rather than trusting the
    wrong SAMPLE_ATTRIBUTE tag silently.
    """
    stats = LoadStats()
    source_uri = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra"
    run_id = _start_run(connection, "NCBI_BioSample_demographics", source_uri)
    counters = {"sex_filled": 0, "age_filled": 0, "bmi_filled": 0, "no_biosample_attrs": 0}
    audit_rows: list[dict[str, Any]] = []
    try:
        cursor = connection.cursor()
        query = (
            "SELECT run_accession FROM samples "
            "WHERE sex IS NULL OR age_years IS NULL OR bmi IS NULL ORDER BY id"
        )
        if limit:
            query += " LIMIT %s"
            cursor.execute(query, (limit,))
        else:
            cursor.execute(query)
        run_accessions = [row[0] for row in cursor.fetchall()]
        cursor.close()

        for batch in batched(run_accessions, batch_size):
            stats.read += len(batch)
            attrs_by_run = client.fetch_run_attributes(batch)
            time.sleep(sleep_seconds)

            current_cursor = connection.cursor()
            placeholders = ", ".join(["%s"] * len(batch))
            current_cursor.execute(
                f"SELECT run_accession, sex, age_years, bmi FROM samples "
                f"WHERE run_accession IN ({placeholders})",
                tuple(batch),
            )
            currently_null = {
                run_accession: {"sex": sex is None, "age_years": age is None, "bmi": bmi is None}
                for run_accession, sex, age, bmi in current_cursor.fetchall()
            }
            current_cursor.close()

            for run_accession in batch:
                attrs = attrs_by_run.get(run_accession)
                if not attrs:
                    counters["no_biosample_attrs"] += 1
                    stats.skipped += 1
                    continue
                null_state = currently_null.get(run_accession, {})

                sex_value = None
                if null_state.get("sex"):
                    for tag in SEX_TAGS:
                        if tag in attrs:
                            sex_value = parse_host_sex(attrs[tag])
                            if sex_value:
                                break

                age_value, age_note = None, ""
                if null_state.get("age_years"):
                    for tag in AGE_TAGS:
                        if tag in attrs:
                            age_value, age_note = parse_host_age_years(attrs[tag])
                            if age_value is not None:
                                break

                bmi_value = None
                if null_state.get("bmi"):
                    for tag in BMI_TAGS:
                        if tag in attrs:
                            bmi_value = parse_bmi(attrs[tag])
                            if bmi_value is not None:
                                break

                if sex_value is None and age_value is None and bmi_value is None:
                    stats.skipped += 1
                    if audit_csv:
                        audit_rows.append({
                            "run_accession": run_accession, "sex": "", "age_years": "",
                            "age_note": age_note, "bmi": "", "has_biosample_attrs": True,
                        })
                    continue

                update_cursor = connection.cursor()
                update_cursor.execute(
                    """
                    UPDATE samples
                    SET sex = COALESCE(sex, %s),
                        age_years = COALESCE(age_years, %s),
                        bmi = COALESCE(bmi, %s)
                    WHERE run_accession = %s
                    """,
                    (sex_value, age_value, bmi_value, run_accession),
                )
                update_cursor.close()
                stats.updated += 1
                counters["sex_filled"] += int(sex_value is not None)
                counters["age_filled"] += int(age_value is not None)
                counters["bmi_filled"] += int(bmi_value is not None)
                if audit_csv:
                    audit_rows.append({
                        "run_accession": run_accession,
                        "sex": sex_value or "",
                        "age_years": age_value if age_value is not None else "",
                        "age_note": age_note,
                        "bmi": bmi_value if bmi_value is not None else "",
                        "has_biosample_attrs": True,
                    })

        if audit_csv:
            _write_demographics_audit(audit_csv, audit_rows)
        _finish_run(connection, run_id, stats)
        return stats, counters
    except Exception as error:
        _finish_run(connection, run_id, stats, error)
        raise


def sync_mgnify_study(
    connection: MySQLConnection,
    client: MGnifyClient,
    study_accession: str,
    project_title: str,
    disease_stage_map: dict[str, tuple[str, str]],
    default_body_site: str,
    samples_output: str | Path,
    abundances_output: str | Path,
    sleep_seconds: float = 0.3,
) -> tuple[LoadStats, LoadStats, dict[str, int]]:
    """Pull a completed MGnify study's sample metadata and precomputed genus/species
    taxonomy, write them as CSVs matching this project's existing templates, and load
    them via the existing load_samples/load_abundances loaders.

    disease_stage_map maps the raw 'gastrointestinal tract disorder' metadata value
    (casefolded) to (disease_name, mesh_id). A stage absent from the map is loaded
    with disease_name left blank rather than guessed.
    """
    samples = client.study_samples(study_accession)
    analyses = client.study_analyses(study_accession)
    analysis_by_sample = best_analysis_per_sample(analyses)

    sample_rows: list[dict[str, Any]] = []
    abundance_rows: list[dict[str, Any]] = []
    counters = {"samples_seen": 0, "samples_no_analysis": 0, "taxa_rows": 0}

    for sample in samples:
        counters["samples_seen"] += 1
        run_accession = sample.get("id")
        if not run_accession:
            continue
        metadata = sample_metadata_dict(sample)
        stage = metadata.get("gastrointestinal tract disorder", "")
        disease_name, mesh_id = disease_stage_map.get(stage.casefold(), ("", ""))
        age_years, _ = parse_host_age_years(metadata.get("host age", ""))
        sex = parse_host_sex(metadata.get("host sex", ""))
        country = metadata.get("geographic location (country and/or sea,region)", "")

        sample_rows.append({
            "project_id": study_accession,
            "project_title": project_title,
            "data_type": "16S rRNA",
            "data_quality": "qualified",
            "source_database": "MGnify",
            "run_id": run_accession,
            "sample_id": clean_str(sample.get("attributes", {}).get("biosample")),
            "disease_name": disease_name,
            "mesh_id": mesh_id,
            "sex": sex or "",
            "age_years": age_years if age_years is not None else "",
            "country": country,
            "qc_status": "",
            "body_site": metadata.get("environment (feature)") or default_body_site,
        })

        analysis = analysis_by_sample.get(run_accession)
        if not analysis:
            counters["samples_no_analysis"] += 1
            continue
        taxonomy = client.analysis_taxonomy_ssu(analysis["id"])
        time.sleep(sleep_seconds)
        taxa = genus_abundances(taxonomy)
        total = sum(count for _, count in taxa)
        for name, count in taxa:
            if total <= 0:
                continue
            abundance_rows.append({
                "run_id": run_accession,
                "scientific_name": name,
                "taxonomic_rank": "genus",
                "ncbi_tax_id": "",
                "relative_abundance": round(count / total, 6),
                "source_database": "MGnify",
            })
            counters["taxa_rows"] += 1

    samples_path = Path(samples_output)
    samples_path.parent.mkdir(parents=True, exist_ok=True)
    sample_fieldnames = [
        "project_id", "project_title", "data_type", "data_quality", "source_database",
        "run_id", "sample_id", "disease_name", "mesh_id", "sex", "age_years", "country",
        "qc_status", "body_site",
    ]
    with samples_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=sample_fieldnames)
        writer.writeheader()
        writer.writerows(sample_rows)

    abundances_path = Path(abundances_output)
    abundances_path.parent.mkdir(parents=True, exist_ok=True)
    abundance_fieldnames = [
        "run_id", "scientific_name", "taxonomic_rank", "ncbi_tax_id",
        "relative_abundance", "source_database",
    ]
    with abundances_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=abundance_fieldnames)
        writer.writeheader()
        writer.writerows(abundance_rows)

    sample_stats = load_samples(connection, samples_path)
    abundance_stats = load_abundances(connection, abundances_path)
    return sample_stats, abundance_stats, counters


def sync_gmrepo_abundances(
    connection: MySQLConnection,
    client: GMrepoClient,
    output_path: str | Path,
    rank: str = "genus",
    limit: int | None = None,
    workers: int = 4,
    only_complete_demographics: bool = False,
    data_type: str = "",
) -> tuple[LoadStats, dict[str, int]]:
    """Fetch per-run taxonomic profiles from GMrepo and load them as sample abundances.

    Only runs that have no abundance rows yet are fetched, so this is resumable
    and safe to re-run after an interruption.

    A single taxonomic rank is loaded per run. GMrepo profiles carry several ranks
    at once, and loading more than one would make a sample's abundances sum to the
    number of ranks rather than to 1. Abundances are renormalized per run so every
    sample sums to 1.0 regardless of whether GMrepo reported percentages or
    fractions, matching the convention already used by the MGnify loader.
    """
    wanted_rank = clean_str(rank).casefold() or "species"
    cursor = connection.cursor()
    # Resumability is per RANK, not per sample: an mNGS run that already holds a
    # genus layer still needs its species layer. Selecting on "has no abundances
    # at all" would skip every sample that has ever been loaded.
    query = """
        SELECT s.run_accession
        FROM samples s
        WHERE NOT EXISTS (
            SELECT 1 FROM sample_taxon_abundances a
            JOIN taxa t ON t.id = a.taxon_id
            WHERE a.sample_id = s.id AND t.taxonomic_rank = %s
        )
    """
    params: list[Any] = [wanted_rank]
    if only_complete_demographics:
        query += " AND s.sex IS NOT NULL AND s.age_years IS NOT NULL"
    if data_type:
        # GMrepo reports no species list at all for 16S runs, so asking for
        # species there is ~18,000 guaranteed-empty requests.
        query += " AND s.study_id IN (SELECT id FROM studies WHERE data_type = %s)"
        params.append(data_type)
    query += " ORDER BY s.id"
    if limit:
        query += f" LIMIT {int(limit)}"
    cursor.execute(query, tuple(params))
    run_accessions = [row[0] for row in cursor.fetchall()]
    cursor.close()

    counters = {
        "runs_requested": len(run_accessions),
        "runs_with_profile": 0,
        "runs_empty": 0,
        "runs_failed": 0,
        "runs_aggregated_to_genus": 0,
        "rows_exported": 0,
    }
    export_rows: list[dict[str, Any]] = []
    unclassified_by_run: dict[str, float] = {}
    failures: list[str] = []

    def fetch(run_accession: str) -> tuple[str, list[dict[str, Any]] | None, str]:
        local = GMrepoClient(client.base_url, timeout=client.timeout)
        try:
            payload = local.full_taxonomic_profile(run_accession)
            return run_accession, GMrepoClient.parse_taxonomic_profile(payload), ""
        except Exception as error:  # network error or unrecognized shape
            return run_accession, None, str(error)[:200]

    with ThreadPoolExecutor(max_workers=max(1, workers)) as executor:
        futures = {executor.submit(fetch, acc): acc for acc in run_accessions}
        for future in as_completed(futures):
            run_accession, rows, error = future.result()
            if rows is None:
                counters["runs_failed"] += 1
                if len(failures) < 10:
                    failures.append(f"{run_accession}: {error}")
                continue
            # GMrepo reports an "Unknown" bin for reads it could not assign. It is
            # not a taxon, so it must not become a row, but it also must not sit in
            # the normalization denominator or the classified taxa would sum to less
            # than 1 by exactly the unassigned share.
            at_rank = [
                r for r in rows
                if clean_str(r.get("taxonomic_rank")).casefold() == wanted_rank
            ]
            unknown_pct = sum(
                nullable_float(r.get("relative_abundance")) or 0.0
                for r in at_rank
                if str(r.get("scientific_name", "")).strip().casefold()
                in {"unknown", "unclassified", "unassigned", "other"}
            )
            selected = [
                r for r in at_rank
                if str(r.get("scientific_name", "")).strip().casefold()
                not in {"unknown", "unclassified", "unassigned", "other"}
            ]
            if selected:
                gross = sum(
                    nullable_float(r.get("relative_abundance")) or 0.0 for r in at_rank
                )
                if gross > 0:
                    unclassified_by_run[run_accession] = round(unknown_pct / gross, 6)
            if not selected and wanted_rank == "genus":
                # GMrepo may report species only. Genus abundance is then the sum of
                # its species, which is well-defined and keeps 16S and shotgun runs on
                # the same footing; species calls from 16S are not reliable enough to
                # use directly, but their genus totals are.
                aggregated: dict[str, float] = {}
                for r in rows:
                    genus, _ = parse_scientific_name(r.get("scientific_name"))
                    value = nullable_float(r.get("relative_abundance"))
                    if genus and value and value > 0:
                        aggregated[genus] = aggregated.get(genus, 0.0) + value
                selected = [
                    {"scientific_name": g, "relative_abundance": v,
                     "taxonomic_rank": "genus", "ncbi_tax_id": None}
                    for g, v in aggregated.items()
                ]
                if selected:
                    counters["runs_aggregated_to_genus"] += 1
            if not selected:
                # requested rank unavailable for this run; skip rather than mixing ranks
                counters["runs_empty"] += 1
                continue
            total = 0.0
            for r in selected:
                value = nullable_float(r.get("relative_abundance"))
                if value and value > 0:
                    total += value
            if total <= 0:
                counters["runs_empty"] += 1
                continue
            counters["runs_with_profile"] += 1
            # GMrepo can list the same taxon more than once in a run, and distinct
            # source names can normalize to one taxon. Both are summed here: the
            # loader upserts on (sample, taxon) and would otherwise let the last
            # row overwrite the earlier one, silently dropping that abundance.
            merged: dict[tuple[str, str], dict[str, Any]] = {}
            for r in selected:
                value = nullable_float(r.get("relative_abundance"))
                if not value or value <= 0:
                    continue
                genus, species = parse_scientific_name(r["scientific_name"])
                if wanted_rank == "genus":
                    # GMrepo's genus array contains a few binomials, e.g.
                    # "[Clostridium] methylpentosum". Stored verbatim they create a
                    # species-rank taxon inside what is meant to be a pure genus
                    # layer, so that sample's genus abundances no longer sum to 1.
                    # The requested rank must be what actually gets written.
                    key = (genus, "")
                    display = genus
                    taxid = ""  # the id belongs to the binomial, not to the genus
                else:
                    key = (genus, species)
                    display = r["scientific_name"]
                    taxid = clean_str(r.get("ncbi_tax_id"))
                if key in merged:
                    merged[key]["value"] += value
                else:
                    merged[key] = {
                        "value": value,
                        "scientific_name": display,
                        "ncbi_tax_id": taxid,
                    }
            for entry in merged.values():
                export_rows.append({
                    "run_id": run_accession,
                    "scientific_name": entry["scientific_name"],
                    "taxonomic_rank": wanted_rank,
                    "ncbi_tax_id": entry["ncbi_tax_id"],
                    "relative_abundance": round(entry["value"] / total, 8),
                    "source_database": "GMrepo",
                })
                counters["rows_exported"] += 1

    if counters["runs_failed"] and not counters["runs_with_profile"]:
        raise RuntimeError(
            f"All {counters['runs_failed']} GMrepo profile requests failed; "
            f"first errors: {'; '.join(failures[:5])}"
        )
    if failures:
        # Surfaced even on a partial failure. Previously these were counted but
        # never shown unless every request failed, so a run could report success
        # while hundreds of samples were silently skipped for an unknown reason.
        print(
            f"  WARNING: {counters['runs_failed']:,} request(s) failed; "
            f"sample errors: {' | '.join(failures[:3])}"
        )

    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "run_id", "scientific_name", "taxonomic_rank", "ncbi_tax_id",
        "relative_abundance", "source_database",
    ]
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(export_rows)

    if unclassified_by_run:
        cur = connection.cursor()
        for run_accession, fraction in unclassified_by_run.items():
            cur.execute(
                "UPDATE samples SET unclassified_fraction=%s WHERE run_accession=%s",
                (fraction, run_accession),
            )
        cur.close()
    stats = load_abundances(connection, destination)
    return stats, counters


def enrich_taxonomy(
    connection: MySQLConnection,
    client: TaxonomyClient,
    limit: int | None = None,
    batch_size: int = 150,
    sleep_seconds: float = 0.4,
    resolve_missing_ids: bool = True,
) -> tuple[LoadStats, dict[str, int]]:
    """Fill taxonomic lineage (and missing ncbi_tax_id) from NCBI Taxonomy.

    Existing values are never overwritten; only NULL columns are filled. Phylum
    names are folded to the classic vocabulary already used in this table, so a
    phylum stays one groupable value rather than splitting across the pre- and
    post-2021 names.
    """
    stats = LoadStats()
    run_id = _start_run(
        connection, "NCBI_Taxonomy",
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy",
    )
    counters = {
        "with_taxid": 0, "names_resolved": 0, "names_unresolved": 0,
        "lineage_filled": 0, "taxid_filled": 0, "no_lineage_returned": 0,
    }
    try:
        cursor = connection.cursor()
        query = """
            SELECT id, scientific_name, ncbi_tax_id, superkingdom
            FROM taxa
            WHERE phylum IS NULL OR class_name IS NULL OR order_name IS NULL
               OR family IS NULL OR superkingdom IS NULL OR ncbi_tax_id IS NULL
            ORDER BY id
        """
        if limit:
            query += f" LIMIT {int(limit)}"
        cursor.execute(query)
        pending = cursor.fetchall()
        cursor.close()
        stats.read = len(pending)

        # resolve names for rows lacking a taxid, verifying each match
        resolved: dict[int, str] = {}
        for taxon_id, name, taxid, known_domain in pending:
            if taxid:
                resolved[taxon_id] = str(taxid)
                counters["with_taxid"] += 1
            elif resolve_missing_ids:
                found = client.resolve_name(name, known_domain or "")
                time.sleep(sleep_seconds)
                if found:
                    resolved[taxon_id] = found
                    counters["names_resolved"] += 1
                else:
                    counters["names_unresolved"] += 1

        by_taxid: dict[str, list[int]] = defaultdict(list)
        for taxon_id, taxid in resolved.items():
            by_taxid[taxid].append(taxon_id)
        original_taxid = {t_id: tx for t_id, _, tx, _sk in pending}

        all_taxids = list(by_taxid)
        for batch in batched(all_taxids, batch_size):
            lineages = client.fetch_lineages(batch)
            time.sleep(sleep_seconds)
            for taxid in batch:
                record = lineages.get(taxid)
                if not record:
                    counters["no_lineage_returned"] += 1
                    continue
                superkingdom = record.get("superkingdom") or record.get("domain") or None
                phylum = normalize_phylum(record.get("phylum")) or None
                for taxon_id in by_taxid[taxid]:
                    cur = connection.cursor()
                    cur.execute(
                        """
                        UPDATE taxa SET
                            superkingdom = COALESCE(superkingdom, %s),
                            phylum       = COALESCE(phylum, %s),
                            class_name   = COALESCE(class_name, %s),
                            order_name   = COALESCE(order_name, %s),
                            family       = COALESCE(family, %s),
                            ncbi_tax_id  = COALESCE(ncbi_tax_id, %s)
                        WHERE id = %s
                        """,
                        (
                            superkingdom, phylum,
                            record.get("class") or None,
                            record.get("order") or None,
                            record.get("family") or None,
                            int(taxid),
                            taxon_id,
                        ),
                    )
                    cur.close()
                    counters["lineage_filled"] += 1
                    if not original_taxid.get(taxon_id):
                        counters["taxid_filled"] += 1
                    stats.updated += 1

        _finish_run(connection, run_id, stats)
        return stats, counters
    except Exception as error:
        _finish_run(connection, run_id, stats, error)
        raise


def clear_contradictory_lineages(connection: MySQLConnection) -> int:
    """Blank lineages that contradict the organism's own superkingdom.

    A prokaryote cannot sit in Chordata or Streptophyta. Such rows come from a
    cross-kingdom homonym ("Proteus" the bacterium vs the salamander) being
    resolved to the wrong branch, and because each column is filled independently
    the result is a chimera: correct bacterial superkingdom, animal phylum. The
    lineage is cleared so a domain-constrained pass can refill it correctly.
    """
    placeholders = ", ".join(["%s"] * len(NON_PROKARYOTIC_PHYLA))
    cursor = connection.cursor()
    cursor.execute(
        f"""
        UPDATE taxa
        SET phylum = NULL, class_name = NULL, order_name = NULL, family = NULL,
            ncbi_tax_id = NULL
        WHERE superkingdom IN ('Bacteria', 'Archaea')
          AND LOWER(phylum) IN ({placeholders})
        """,
        tuple(NON_PROKARYOTIC_PHYLA),
    )
    changed = cursor.rowcount
    cursor.close()
    return changed


def normalize_existing_phyla(connection: MySQLConnection) -> int:
    """Fold already-stored post-2021 phylum names onto the classic name."""
    cursor = connection.cursor()
    cursor.execute("SELECT DISTINCT phylum FROM taxa WHERE phylum IS NOT NULL")
    changed = 0
    for (phylum,) in cursor.fetchall():
        canonical = normalize_phylum(phylum)
        if canonical and canonical != phylum:
            update = connection.cursor()
            update.execute("UPDATE taxa SET phylum = %s WHERE phylum = %s", (canonical, phylum))
            changed += update.rowcount
            update.close()
    cursor.close()
    return changed


def enrich_cmd_demographics(
    connection: MySQLConnection,
    client: CMDCurationClient,
    audit_csv: str | Path | None = None,
    overwrite: bool = False,
) -> tuple[LoadStats, dict[str, int]]:
    """Fill sample sex/age/BMI from curatedMetagenomicData's curated metadata.

    Only NULL fields are filled. Where cMD and an existing value disagree the
    stored value wins and the disagreement is counted and written to the audit
    file: cMD is a second-hand transcription from the paper, so a conflict is
    something to look at rather than something to silently resolve.
    """
    stats = LoadStats()
    run_id = _start_run(connection, "cMD_curation", CMD_RAW)
    counters = {
        "cmd_studies": 0, "rows_seen": 0, "matched": 0,
        "sex_filled": 0, "age_filled": 0, "bmi_filled": 0,
        "sex_conflicts": 0, "age_conflicts": 0, "age_rejected": 0,
    }
    audit_rows: list[dict[str, Any]] = []
    try:
        cursor = connection.cursor()
        cursor.execute("SELECT run_accession, id, sex, age_years, bmi FROM samples")
        ours = {
            r[0]: {"id": r[1], "sex": r[2], "age": r[3], "bmi": r[4]}
            for r in cursor.fetchall()
        }
        cursor.close()

        for study in client.studies():
            rows = client.study_samples(study)
            if not rows:
                continue
            counters["cmd_studies"] += 1
            for row in rows:
                counters["rows_seen"] += 1
                for accession in run_accessions(row):
                    current = ours.get(accession)
                    if not current:
                        continue
                    counters["matched"] += 1
                    stats.read += 1

                    sex = parse_sex(row)
                    age, age_note = parse_age_years(row)
                    bmi = parse_cmd_bmi(row)
                    if age is None and age_note.startswith(("out_of_range", "unknown_unit")):
                        counters["age_rejected"] += 1

                    set_sex = set_age = set_bmi = None
                    if sex:
                        if current["sex"] is None:
                            set_sex = sex
                            counters["sex_filled"] += 1
                        elif sex != str(current["sex"]).casefold():
                            counters["sex_conflicts"] += 1
                    if age is not None:
                        if current["age"] is None:
                            set_age = age
                            counters["age_filled"] += 1
                        elif abs(float(current["age"]) - age) > 1.0:
                            counters["age_conflicts"] += 1
                    if bmi is not None and current["bmi"] is None:
                        set_bmi = bmi
                        counters["bmi_filled"] += 1

                    if set_sex is None and set_age is None and set_bmi is None:
                        stats.skipped += 1
                    else:
                        update = connection.cursor()
                        update.execute(
                            """UPDATE samples SET
                                 sex = COALESCE(sex, %s),
                                 age_years = COALESCE(age_years, %s),
                                 bmi = COALESCE(bmi, %s)
                               WHERE id = %s""",
                            (set_sex, set_age, set_bmi, current["id"]),
                        )
                        update.close()
                        stats.updated += 1

                    if audit_csv:
                        audit_rows.append({
                            "run_accession": accession, "cmd_study": study,
                            "cmd_sex": sex or "", "cmd_age_years": age if age is not None else "",
                            "age_note": age_note, "cmd_bmi": bmi if bmi is not None else "",
                            "filled_sex": set_sex or "", "filled_age": set_age or "",
                            "filled_bmi": set_bmi or "",
                        })

        if audit_csv and audit_rows:
            destination = Path(audit_csv)
            destination.parent.mkdir(parents=True, exist_ok=True)
            with destination.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(audit_rows[0]))
                writer.writeheader()
                writer.writerows(audit_rows)

        _finish_run(connection, run_id, stats)
        return stats, counters
    except Exception as error:
        _finish_run(connection, run_id, stats, error)
        raise


def validate_database(connection: MySQLConnection) -> list[tuple[str, int]]:
    checks = {
        "taxa": "SELECT COUNT(*) FROM taxa",
        "diseases": "SELECT COUNT(*) FROM diseases",
        "studies": "SELECT COUNT(*) FROM studies",
        "associations": "SELECT COUNT(*) FROM taxon_disease_associations",
        "samples": "SELECT COUNT(*) FROM samples",
        "abundance_rows": "SELECT COUNT(*) FROM sample_taxon_abundances",
        "taxa_missing_genus": "SELECT COUNT(*) FROM taxa WHERE genus = '' OR genus IS NULL",
        "orphan_associations": """
            SELECT COUNT(*) FROM taxon_disease_associations a
            LEFT JOIN taxa t ON t.id = a.taxon_id
            LEFT JOIN diseases d ON d.id = a.disease_id
            WHERE t.id IS NULL OR d.id IS NULL
        """,
        "orphan_abundances": """
            SELECT COUNT(*) FROM sample_taxon_abundances a
            LEFT JOIN taxa t ON t.id = a.taxon_id
            LEFT JOIN samples s ON s.id = a.sample_id
            WHERE t.id IS NULL OR s.id IS NULL
        """,
        # A sample may carry more than one taxonomic rank (genus for every sample,
        # species additionally where the sequencing supports it), so abundances sum
        # to 1 within a rank rather than across the whole sample. Summing a sample
        # without filtering rank returns the number of ranks present, not 1.
        "abundance_rank_groups_not_summing_to_1": """
            SELECT COUNT(*) FROM (
                SELECT a.sample_id, t.taxonomic_rank, SUM(a.relative_abundance) AS total
                FROM sample_taxon_abundances a
                JOIN taxa t ON t.id = a.taxon_id
                GROUP BY a.sample_id, t.taxonomic_rank
                HAVING ABS(total - 1.0) > 0.01
            ) bad
        """,
    }
    cursor = connection.cursor()
    results = []
    for label, sql in checks.items():
        cursor.execute(sql)
        results.append((label, int(cursor.fetchone()[0])))
    cursor.close()
    return results
