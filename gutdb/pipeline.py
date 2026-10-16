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
    normalize_platform,
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
    resolve_effect_type,
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
        # Trust the extension before guessing. Sniffer infers the delimiter from
        # the first 8KB and picks a character that merely appears often, so a
        # free-text notes column containing semicolons or pipes can make it read
        # a comma-separated file with the wrong delimiter. That is how all 15
        # rows of literature_preeclampsia_markers.csv loaded with their columns
        # shifted, putting a sentence from `notes` into source_database -- and it
        # fails silently, because every field still has a value.
        suffix = source.suffix.lower()
        if suffix in {".tsv", ".tab"}:
            dialect: Any = csv.excel_tab
        elif suffix == ".csv":
            dialect = csv.excel
        else:
            try:
                dialect = csv.Sniffer().sniff(sample, delimiters=",\t;|")
            except csv.Error:
                dialect = csv.excel
        reader = csv.DictReader(handle, dialect=dialect)
        # A single-field header means the chosen delimiter is absent from the
        # file, so fall back to sniffing rather than loading one giant column.
        if reader.fieldnames is not None and len(reader.fieldnames) == 1:
            handle.seek(0)
            try:
                reader = csv.DictReader(
                    handle, dialect=csv.Sniffer().sniff(sample, delimiters=",\t;|"))
            except csv.Error:
                handle.seek(0)
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
    if source == "microbiomehd":
        return "MicrobiomeHD"
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


def _platform_for(
    connection: MySQLConnection, project: str, incoming: str
) -> str | None:
    """Reconcile a study's platform with what is already recorded for it."""
    if not incoming:
        return None
    cursor = connection.cursor()
    cursor.execute(
        "SELECT data_type FROM studies WHERE project_accession = %s", (project,))
    existing = cursor.fetchone()
    cursor.close()
    current = clean_str(existing[0]) if existing and existing[0] else ""
    if not current or current == incoming:
        return incoming
    if current == "mixed":
        return "mixed"
    return "mixed"


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
            data_quality = IF(%s, VALUES(data_quality), data_quality)
        """,
        (
            project,
            first_value(row, "project_title", "title") or None,
            first_value(row, "project_description", "description") or None,
            # Platform is normalized, and a study whose runs disagree becomes
            # 'mixed' rather than taking whichever row happened to be written
            # last. Ten GMrepo projects hold both amplicon and shotgun runs, and
            # last-write-wins had been labelling each of them by one half: the
            # amplicon runs of a study stamped mNGS then looked like shotgun
            # samples with no taxonomic profile.
            _platform_for(connection, project,
                          normalize_platform(first_value(row, "data_type", "experiment_type"))),
            first_value(row, "source_database") or "GMrepo",
            first_value(row, "data_quality") if first_value(row, "data_quality") in {"curated", "qualified"} else "unknown",
            bool(first_value(row, "data_quality")),
        ),
    )
    cursor.execute("SELECT id FROM studies WHERE project_accession = %s", (project,))
    study_id = cursor.fetchone()[0]
    cursor.close()
    return study_id


def _stable_method(
    connection: MySQLConnection,
    project: str,
    phenotype_a: str,
    phenotype_b: str,
    rank: str,
    declared: str,
) -> str:
    """Resolve the method for a comparison without inventing a new identity.

    `method` is part of comparison_key, and it has to be, because one study can
    legitimately hold two analyses of the same arms: MicrobiomeHD's standardized
    re-analysis of PMID 25432777 sits beside that paper's own reported result,
    and collapsing them would destroy the comparison between re-analysis and
    publication that the two were loaded to enable.

    What it must not do is change value for the same data. It used to fall back
    to `effect_type` when a file declared no method, so when the PRJNA705217
    marker export lost its method column it resolved to "LDA" where the first
    load had recorded "LEfSe" -- a second comparison key for one comparison, and
    all 34 of its IBS associations duplicated beneath it.

    A file that declares a method is taken at its word. One that does not reuses
    the method already recorded for this study, arms and rank, so a reload
    attaches to the comparison it attached to before.
    """
    declared = clean_str(declared)
    if declared:
        return declared
    cursor = connection.cursor()
    cursor.execute(
        """SELECT c.method FROM phenotype_comparisons c
           JOIN studies s ON s.id = c.study_id
           JOIN diseases da ON da.id = c.phenotype_a_id
           JOIN diseases db ON db.id = c.phenotype_b_id
           WHERE s.project_accession = %s AND da.name = %s AND db.name = %s
             AND c.taxonomic_level = %s AND c.method IS NOT NULL
           ORDER BY c.id LIMIT 1""",
        (project, phenotype_a, phenotype_b, rank),
    )
    existing = cursor.fetchone()
    cursor.close()
    if existing and clean_str(existing[0]):
        return clean_str(existing[0])
    return "LEfSe"


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
            method = _stable_method(
                connection, project_id, phenotype_a, phenotype_b, rank,
                first_value(row, "method"),
            )
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
                (taxon_id, disease_id, comparison_id,
                 resolve_effect_type(first_value(row, "effect_type"), lda)),
            )
            existed = cursor.fetchone() is not None
            cursor.execute(
                """
                INSERT INTO taxon_disease_associations
                    (taxon_id, disease_id, comparison_id, direction, effect_type, effect_size,
                     p_value, q_value, source_database)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
                ON DUPLICATE KEY UPDATE
                    direction = VALUES(direction),
                    effect_size = IF(%s, VALUES(effect_size), COALESCE(effect_size, VALUES(effect_size))),
                    p_value = COALESCE(VALUES(p_value), p_value),
                    q_value = COALESCE(VALUES(q_value), q_value),
                    source_database = VALUES(source_database)
                """,
                (
                    taxon_id, disease_id, comparison_id, direction,
                    resolve_effect_type(first_value(row, "effect_type"), lda), lda,
                    # The hand-curated marker files carry these and they were
                    # being dropped: 59 p values and 42 q values across the 15
                    # files. The columns existed only after p_value/q_value were
                    # added for MicrobiomeHD, so wiring them here is what makes
                    # that schema addition reach the literature set too.
                    nullable_float(first_value(row, "p_value", "pvalue", "p")),
                    nullable_float(first_value(row, "q_value", "qvalue", "fdr", "q")),
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
    "project_id", "project_title", "data_type", "data_quality",
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
    "project_id", "data_type", "run_id", "sample_id", "disease_name",
    "mesh_id", "sex", "age_years", "bmi", "country", "qc_status", "instrument_model",
    "nr_reads_sequenced",
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
    """Pick the phenotype a run belongs to, or nothing when GMrepo is ambiguous.

    Preferring the non-Health term is right while GMrepo lists the one phenotype
    that applies to a run. It now lists every phenotype associated with the run,
    so a control in a low-birth-weight cohort arrives as
    [Health, Infant Low Birth Weight] and this would call it a case. A re-sync on
    that basis relabelled roughly 14,000 controls as cases, which silently
    destroys the case/control structure every comparison in this database rests
    on -- and it fails open, because a wrong label looks exactly like a right one.

    There is no field in the runs endpoint that says which arm a run is in, so
    an ambiguous list yields nothing rather than a guess. load_samples upserts
    disease_id through COALESCE, so an existing label survives and a genuinely
    new sample is left unlabelled instead of mislabelled.
    """
    phenotypes = _parse_gmrepo_phenotypes(value)
    if not phenotypes:
        return "", ""
    terms = [clean_str(row.get("term")) for row in phenotypes]
    has_health = any(term.casefold() == "health" for term in terms)
    if len(phenotypes) > 1 and has_health:
        return "", ""
    chosen = next(
        (
            row for row in phenotypes
            if clean_str(row.get("term")).casefold() != "health"
        ),
        phenotypes[0],
    )
    return clean_str(chosen.get("term")), clean_str(chosen.get("disease"))


def prepare_gmrepo_sample_export(
    connection: MySQLConnection, rows: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    """Reconcile an export with stored sample metadata using exact accessions.

    Existing database demographics/labels include cMD/NCBI enrichment and
    curation that the raw GMrepo API does not carry. Preserve that work in the
    CSV as well. Missing measurements stay blank; no statistical imputation.
    Raw phenotype JSON is parsed transiently, never retained as an export field.
    """
    cursor = connection.cursor(dictionary=True)
    try:
        cursor.execute("""SELECT s.run_accession, s.gmrepo_sample_id,
            st.project_accession, s.sex, s.age_years, s.bmi,
            d.name AS disease_name, d.mesh_id
            FROM samples s JOIN studies st ON st.id = s.study_id
            LEFT JOIN diseases d ON d.id = s.disease_id""")
        stored = {r["run_accession"]: r for r in cursor.fetchall()}
    finally:
        cursor.close()
    output, audit = [], []
    for raw in rows:
        run = clean_str(raw.get("run_id"))
        current = stored.get(run, {})
        if current and (
            clean_str(current["project_accession"]) != clean_str(raw.get("project_id"))
            or (clean_str(current["gmrepo_sample_id"]) and clean_str(raw.get("sample_id"))
                and clean_str(current["gmrepo_sample_id"]) != clean_str(raw["sample_id"]))
        ):
            raise ValueError(f"Sample identity mismatch for {run}; export was not written")
        row = {field: raw.get(field, "") for field in GMREPO_SAMPLE_FIELDS}
        sources = {}
        for field, parser in (("sex", parse_host_sex),
                              ("age_years", lambda x: _plausible_age(nullable_float(x))),
                              ("bmi", lambda x: _plausible_bmi(nullable_float(x)))):
            original = parser(raw.get(field))
            recovered = parser(current.get(field))
            value = recovered if recovered is not None else original
            row[field] = "" if value is None else str(value)
            sources[field] = ("mysql.samples" if recovered is not None
                              else "source_validation")
        parsed_name, parsed_mesh = _choose_sample_phenotype(raw.get("phenotypes"))
        if not clean_str(row["disease_name"]) and not clean_str(row["mesh_id"]):
            row["disease_name"], row["mesh_id"] = parsed_name, parsed_mesh
        for field in ("disease_name", "mesh_id"):
            if current.get(field) is not None:
                row[field] = current[field]
                sources[field] = "mysql.diseases"
            else:
                row[field] = clean_str(row[field])
                sources[field] = "phenotype_json_or_existing_label"
        for field in ("sex", "age_years", "bmi", "disease_name", "mesh_id"):
            before, after = clean_str(raw.get(field)), clean_str(row[field])
            if before != after:
                audit.append({"run_id": run, "field": field, "before": before,
                              "after": after, "source": sources[field]})
        output.append(row)
    return output, audit


def write_sample_export_audit(path: str | Path, rows: list[dict[str, str]]) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["run_id", "field", "before", "after", "source"])
        writer.writeheader()
        writer.writerows(rows)


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
    # GMrepo returns each run twice from getAllRunsByProjectIDAsync, the two
    # records differing only in an internal id, so the export is keyed on run_id
    # rather than appended to. Without this every run appears twice: a re-sync
    # produced 67,160 rows for 33,580 runs.
    by_run: dict[str, dict[str, Any]] = {}
    for project_id in projects:
        skip = 0
        while True:
            page = client.project_runs(project_id, limit=limit, skip=skip)
            if not page:
                break
            for raw in page:
                disease_name, mesh_id = _choose_sample_phenotype(raw.get("phenotypes"))
                run_id = clean_str(raw.get("run_id"))
                if not run_id or run_id in by_run:
                    continue
                by_run[run_id] = (
                    {
                        "project_id": clean_str(raw.get("project_id")) or project_id,
                        "data_type": clean_str(raw.get("experiment_type")),
                        "run_id": clean_str(raw.get("run_id")),
                        "sample_id": clean_str(raw.get("sample_id")),
                        "disease_name": disease_name,
                        "mesh_id": mesh_id,
                        "sex": clean_str(raw.get("sex")),
                        "age_years": raw.get("host_age"),
                        # GMrepo exposes BMI and the export used to drop it,
                        # leaving values on the server for studies where this
                        # database had none.
                        "bmi": raw.get("BMI"),
                        "country": clean_str(raw.get("country")),
                        "qc_status": raw.get("QCStatus"),
                        "instrument_model": clean_str(raw.get("instrument_model")),
                        "nr_reads_sequenced": raw.get("nr_reads_sequenced"),
                    }
                )
            if len(page) < limit:
                break
            skip += limit

    rows, audit = prepare_gmrepo_sample_export(connection, list(by_run.values()))
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GMREPO_SAMPLE_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    write_sample_export_audit(
        destination.with_name(destination.stem + "_enrichment_audit.csv"), audit
    )

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


def _plausible_bmi(value: float | None) -> float | None:
    """Reject a BMI that cannot be a measurement.

    GMrepo records absent BMI as 0 for some cohorts -- 52 rows of PRJEB6172,
    whose real BMIs run 16.3 to 31.8. Loading 0 as a value makes a missing
    measurement look like an impossible one, and the sentinel came back every
    time that export was reloaded, so it is rejected here rather than cleaned up
    afterwards.
    """
    if value is None or value <= 0 or value > 200:
        return None
    return value


def _plausible_age(value: float | None) -> float | None:
    """Reject an age that cannot be a measurement.

    Zero is kept: it is a real birth-day sample in PRJNA716780, which holds 47
    samples under two years old. Negative and absurd values are not -- GMrepo
    serves -88 for two runs of PRJEB42155 and 119 for one of PRJNA521587.
    """
    if value is None or value < 0 or value > 110:
        return None
    return value


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
                    _plausible_age(nullable_float(first_value(row, "age_years", "age"))),
                    _plausible_bmi(nullable_float(first_value(row, "bmi"))),
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
                    (sample_id, taxon_id, relative_abundance)
                VALUES (%s, %s, %s)
                ON DUPLICATE KEY UPDATE relative_abundance = VALUES(relative_abundance)
                """,
                (sample[0], taxon_id, abundance),
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
        # has_biosample_attrs was dropped: a row is only written when NCBI
        # returned attributes for that run, so the flag was hardcoded True at
        # both write sites and said only "this row exists".
        "run_accession", "sex", "age_years", "age_note", "bmi",
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
                            "age_note": age_note, "bmi": "",
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
                        # Record the note only when it says something the row does
                        # not already show. "unit=year" and "empty" restate
                        # cmd_age_years being filled or blank, and filled every row
                        # of a year-only cohort with noise. A rejection
                        # ("out_of_range:140.0", "unknown_unit:decades") or a real
                        # unit conversion ("unit=month") is not recoverable from the
                        # value, so those are kept.
                        informative_note = (
                            "" if age_note in ("empty", "unit=year") else age_note
                        )
                        audit_rows.append({
                            "run_accession": accession, "cmd_study": study,
                            "cmd_sex": sex or "", "cmd_age_years": age if age is not None else "",
                            "age_note": informative_note,
                            "cmd_bmi": bmi if bmi is not None else "",
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
    checks = _validate_checks()
    cursor = connection.cursor()
    results = []
    for label, sql in checks.items():
        cursor.execute(sql)
        row = cursor.fetchone()
        # A tuple key means one query returning one row of several counts, so an
        # expensive view is materialised once while each invariant keeps a name.
        if isinstance(label, tuple):
            if len(row) != len(label):
                raise ValueError(
                    f"check {label[0]!r} returned {len(row)} columns "
                    f"for {len(label)} labels")
            results.extend((name, int(value or 0))
                           for name, value in zip(label, row))
        else:
            results.append((label, int(row[0] or 0)))
    cursor.close()
    return results


def _validate_checks() -> dict[str, str]:
    """The integrity checks, shared by both stores so neither drifts."""
    return {
        "taxa": "SELECT COUNT(*) FROM taxa",
        # A quarter of the taxa table is MiMeDB's trait reference: organisms
        # loaded for their gram stain, oxygen requirement and pathogenicity so
        # those are already there if a study ever reports them. They are not
        # observations, and reading "taxa" as the number of organisms this
        # database has seen overstates it by 1,416. Both numbers are reported so
        # neither has to be inferred.
        "taxa_observed": """
            SELECT COUNT(*) FROM taxa t
            WHERE EXISTS (SELECT 1 FROM sample_taxon_abundances ab
                          WHERE ab.taxon_id = t.id)
               OR EXISTS (SELECT 1 FROM taxon_disease_associations a
                          WHERE a.taxon_id = t.id)
        """,
        "taxa_reference_only": """
            SELECT COUNT(*) FROM taxa t
            WHERE NOT EXISTS (SELECT 1 FROM sample_taxon_abundances ab
                              WHERE ab.taxon_id = t.id)
              AND NOT EXISTS (SELECT 1 FROM taxon_disease_associations a
                              WHERE a.taxon_id = t.id)
        """,
        "diseases": "SELECT COUNT(*) FROM diseases",
        "studies": "SELECT COUNT(*) FROM studies",
        "associations": "SELECT COUNT(*) FROM taxon_disease_associations",
        "samples": "SELECT COUNT(*) FROM samples",
        "abundance_rows": "SELECT COUNT(*) FROM sample_taxon_abundances",
        "taxa_missing_genus": "SELECT COUNT(*) FROM taxa WHERE genus = '' OR genus IS NULL",
        # Gram stain is a property of the cell envelope, so it is near-constant
        # within a phylum, and a lone dissenter is almost always a bad source
        # value rather than a discovery. MiMeDB called Bacteroides vulgatus
        # gram-positive, and the trait backfill moved that onto Phocaeicola
        # vulgatus -- 10,263 abundance rows -- where it was the only positive
        # among 109 Bacteroidetes. On an orphan row nobody would have seen it.
        #
        # Deliberately conservative, because some phyla genuinely hold both:
        # Firmicutes contains the gram-negative Negativicutes, so it is well
        # under the 90% threshold and never flagged. Only 'positive'/'negative'
        # vote and are judged -- 'variable' and 'uncharacterized' are real
        # answers about an organism, not disagreements with its phylum.
        # Written with CTEs and no SUBSTRING_INDEX/GROUP_CONCAT so it runs on
        # the SQLite mirror too. The MySQL-only form left this check unable to
        # run there -- and the mirror turned out to hold six of the seven
        # dissenters, including Parabacteroides distasonis with 9,291 abundance
        # rows, because the repair had reached MySQL only.
        "gram_stain_against_phylum": """
            WITH votes AS (
                SELECT phylum, LOWER(TRIM(gram_stain)) AS g, COUNT(*) AS n
                FROM taxa
                WHERE phylum IS NOT NULL AND TRIM(phylum) <> ''
                  AND LOWER(TRIM(gram_stain)) IN ('positive', 'negative')
                GROUP BY phylum, LOWER(TRIM(gram_stain))
            ),
            totals AS (
                SELECT phylum, SUM(n) AS total, MAX(n) AS top
                FROM votes GROUP BY phylum
            ),
            majority AS (
                SELECT v.phylum, v.g
                FROM votes v
                JOIN totals t ON t.phylum = v.phylum AND v.n = t.top
                WHERE t.total >= 10
                  AND CAST(t.top AS REAL) / t.total >= 0.90
            )
            SELECT COUNT(*) FROM taxa t
            JOIN majority m ON m.phylum = t.phylum
            WHERE LOWER(TRIM(t.gram_stain)) IN ('positive', 'negative')
              AND LOWER(TRIM(t.gram_stain)) <> m.g
        """,
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
        # --- view cardinality ------------------------------------------------
        #
        # Nothing in the test suite reaches a SQL view: the tests are pure
        # functions so they need no credentials, and the views carry some of the
        # most load-bearing logic here. v_taxon_disease_evidence reversed 33
        # directions when it was wrong, v_ml_ready_samples overclaimed what its
        # flag guaranteed, and two defects in v_attribute_concept_coverage were
        # caught by reading its output rather than by any check.
        #
        # These counters close the gap for the failure a view hides best. A view
        # that LEFT JOINs a table on a non-unique key FANS OUT -- one base row
        # becomes several -- and every count, mean and fraction computed from it
        # is then silently wrong while each individual row still looks correct.
        # v_sample_attributes_canonical joins two alias tables and would be the
        # first to go; a duplicate (attribute, raw_value) is barred by its
        # primary key today, and this is what notices if that stops being true.
        #
        # Each counter is a DIFFERENCE from the base table's own count, so zero
        # is the only passing value and the pre-push hook already gates on it.
        "v_ml_ready_samples_row_drift": """
            SELECT ABS((SELECT COUNT(*) FROM v_ml_ready_samples)
                       - (SELECT COUNT(*) FROM samples))
        """,
        "v_sample_attributes_canonical_row_drift": """
            SELECT ABS((SELECT COUNT(*) FROM v_sample_attributes_canonical)
                       - (SELECT COUNT(*) FROM sample_attributes))
        """,
        # Every abundance row belongs to exactly one rank view, so the two must
        # partition the table: a gap means a rank nothing exposes, an overlap
        # means a row counted twice.
        "abundance_rank_views_not_partitioning": """
            SELECT ABS((SELECT COUNT(*) FROM v_abundance_genus)
                       + (SELECT COUNT(*) FROM v_abundance_species)
                       - (SELECT COUNT(*) FROM sample_taxon_abundances))
        """,
        "v_sample_attribute_keys_row_drift": """
            SELECT ABS((SELECT COUNT(*) FROM v_sample_attribute_keys)
                       - (SELECT COUNT(DISTINCT attribute) FROM sample_attributes))
        """,
        "v_taxon_disease_summary_row_drift": """
            SELECT ABS((SELECT COUNT(*) FROM v_taxon_disease_summary)
                       - (SELECT COUNT(*) FROM taxon_disease_associations))
        """,
        # n_taxa_profiled counts rows across every rank at once; the per-rank
        # columns must account for all of them, or a filter written against one
        # of them silently misses a rank.
        "ml_ready_per_rank_counts_not_adding_up": """
            SELECT COUNT(*) FROM v_ml_ready_samples
            WHERE n_genera_profiled + n_species_profiled <> n_taxa_profiled
        """,
        # A boolean concept whose value the alias map does not cover would let a
        # modeller pool a key with an inverted encoding -- the error the alias
        # tables exist to prevent.
        "attribute_booleans_undecoded": """
            SELECT COUNT(*) FROM v_sample_attributes_canonical
            WHERE boolean_undecoded = 1
        """,
        # --- v_taxon_disease_evidence ----------------------------------------
        #
        # The view that reported Faecalibacterium prausnitzii as ENRICHED in
        # ulcerative colitis, reversing one of the field's most replicated
        # findings, and 32 others with it. It had no invariant until now, which
        # is backwards: it is the most consequential view here and the only one
        # known to have been catastrophically wrong.
        #
        # ONLY THE FIRST OF THESE REPRODUCES THAT BUG. The old view pooled
        # disease-vs-disease comparisons with case/control ones, so all 753 pure
        # disease-vs-disease pairs carried a direction built from the wrong
        # question; direction_on_disease_vs_disease fires on every one of them.
        # The second rule it broke -- two independent COUNT(DISTINCT study_id),
        # letting one study be counted on both sides -- cannot be caught on
        # today's data, because no (taxon, disease, study) group has vs-Health
        # associations that disagree with each other, so n_studies_split is 0
        # everywhere and the old form would double-count nothing. The partition
        # counter is therefore a guard against a future regression rather than a
        # test of the historical defect, and worth keeping on exactly those
        # terms: two comparison methods inside one study could disagree tomorrow.
        #
        # ONE QUERY, EIGHT COUNTERS. A tuple key returns one row of counts, so
        # this CTE-heavy view is materialised once per store while each invariant
        # keeps its own name; a single combined total would say a direction was
        # wrong without saying which rule it broke.
        #
        # Measured, because the first version of this comment guessed and was
        # wrong: eight separate queries cost 0.89s against MySQL and the combined
        # one 0.16s, so the collapse saves 0.73s -- real, and worth having, but
        # NOT the reason validate takes a minute. That is the mirror pass, which
        # is 56% of the total; see validate_sqlite_mirror.
        (
            "evidence_direction_on_disease_vs_disease",
            "evidence_study_votes_not_partitioning",
            "evidence_votes_exceed_studies",
            "evidence_direction_without_a_case_control_vote",
            "evidence_missing_direction_despite_votes",
            "evidence_direction_contradicts_its_counts",
            "evidence_agreement_ratio_out_of_range",
            "evidence_row_drift",
        ): """
            SELECT
              SUM(CASE WHEN contrast_scope = 'disease_vs_disease'
                        AND consensus_direction IS NOT NULL
                       THEN 1 ELSE 0 END),
              -- One study, one vote: the three buckets must partition the
              -- case/control studies exactly.
              SUM(CASE WHEN n_studies_enriched + n_studies_depleted
                            + n_studies_split <> n_studies_case_control
                       THEN 1 ELSE 0 END),
              SUM(CASE WHEN n_studies_case_control + n_studies_other_contrast
                            <> n_studies
                       THEN 1 ELSE 0 END),
              -- A direction exists exactly when a case/control study voted,
              -- checked both ways round so neither a spurious direction nor a
              -- lost one passes.
              SUM(CASE WHEN consensus_direction IS NOT NULL
                        AND n_studies_enriched + n_studies_depleted = 0
                       THEN 1 ELSE 0 END),
              SUM(CASE WHEN consensus_direction IS NULL
                        AND n_studies_enriched + n_studies_depleted > 0
                       THEN 1 ELSE 0 END),
              -- The direction must agree with the counts it summarises. This is
              -- what catches a CASE arm written the wrong way round, the single
              -- edit that would silently reverse every direction in the table.
              SUM(CASE WHEN (consensus_direction = 'enriched'
                             AND n_studies_enriched <= n_studies_depleted)
                         OR (consensus_direction = 'depleted'
                             AND n_studies_depleted <= n_studies_enriched)
                         OR (consensus_direction = 'tied'
                             AND n_studies_enriched <> n_studies_depleted)
                       THEN 1 ELSE 0 END),
              -- agreement_ratio is a majority over a non-empty vote, so it
              -- cannot fall below half or above one.
              SUM(CASE WHEN agreement_ratio IS NOT NULL
                        AND (agreement_ratio < 0.5 OR agreement_ratio > 1.0)
                       THEN 1 ELSE 0 END),
              ABS(COUNT(*) - (SELECT COUNT(*) FROM (
                    SELECT DISTINCT taxon_id, disease_id
                    FROM taxon_disease_associations) pairs))
            FROM v_taxon_disease_evidence
        """,
        # --- the LDA sign order ----------------------------------------------
        #
        # An LDA score's sign is meaningless without its comparison order, which
        # is why phenotype_comparisons records which phenotype a POSITIVE score
        # refers to. These two counters are the only things that make that
        # design self-enforcing.
        #
        # One association of 16,400 had a magnitude and no sign order -- 
        # Bacteroides / Vitiligo at LDA 4.31 -- so the 4.31 could not be read.
        # Found by hand; this is what finds the next one.
        #
        # 56 comparisons carry no sign order and that is FINE: the other 55 hold
        # only unscored associations, a direction with no magnitude, where sign
        # order has nothing to disambiguate. So the counter asks about scored
        # associations, not about comparisons.
        "scored_associations_without_a_sign_order": """
            SELECT COUNT(*) FROM taxon_disease_associations a
            JOIN phenotype_comparisons c ON c.id = a.comparison_id
            WHERE a.effect_size IS NOT NULL
              AND c.positive_score_enriched_in_id IS NULL
              AND c.negative_score_enriched_in_id IS NULL
        """,
        # The arm a score's sign points at must be the disease exactly when the
        # direction says enriched. Written as an XOR rather than two cases
        # because the obvious form -- "does the signed arm equal the disease" --
        # is TRUE for enriched and FALSE for depleted, and reading a 0 there as a
        # failure is a mistake easily made: 9,183 correct rows look like
        # disagreements under it.
        "direction_disagrees_with_the_signed_arm": """
            SELECT COUNT(*) FROM taxon_disease_associations a
            JOIN phenotype_comparisons c ON c.id = a.comparison_id
            WHERE a.effect_size IS NOT NULL AND a.effect_size <> 0
              AND a.direction IN ('enriched', 'depleted')
              AND CASE WHEN a.effect_size > 0
                       THEN c.positive_score_enriched_in_id
                       ELSE c.negative_score_enriched_in_id END IS NOT NULL
              AND ((CASE WHEN a.effect_size > 0
                         THEN c.positive_score_enriched_in_id
                         ELSE c.negative_score_enriched_in_id END = a.disease_id)
                   <> (a.direction = 'enriched'))
        """,
        # ONE ORGANISM COUNTED TWICE IN ONE SAMPLE. Two taxa rows sharing a tax
        # ID are the same organism, so a sample holding both counts it twice --
        # it takes double share while every other taxon there is squeezed. The
        # sum-to-1 counter CANNOT see this, because renormalisation happened
        # after the duplication, which is why the mirror carried 582 such cells
        # across 12,392 of its 23,494 samples while every one of validate's other
        # counters read zero.
        #
        # This gates on the harm rather than on duplicate tax IDs, because five
        # duplicate groups legitimately remain in the mirror and are listed in
        # data/taxa_merge_mirror_audit.csv: one whose OVERRIDES action names
        # MySQL row ids, one where two different Vibrio species share a tax ID
        # (a mis-assignment, not a duplicate), and three bracket-notation
        # artifacts that sit at a different rank from their twin and so must
        # never be merged. None of them is counted twice in any sample.
        # GROUPED PER RANK, which the first version of this counter was not. The
        # genus and species layers are separate compositions, each summing to 1,
        # so one organism appearing in both is normal -- that is how every
        # organism works. Only two rows at the SAME rank sharing a tax ID in one
        # sample are the same measurement counted twice. Grouping across ranks
        # read 0 by luck and would have fired the moment a genus-layer feature
        # was given the tax ID its species-layer twin already had.
        "taxon_counted_twice_in_a_sample": """
            SELECT COUNT(*) FROM (
                SELECT a.sample_id, t.taxonomic_rank, t.ncbi_tax_id
                FROM sample_taxon_abundances a
                JOIN taxa t ON t.id = a.taxon_id
                WHERE t.ncbi_tax_id IS NOT NULL
                GROUP BY a.sample_id, t.taxonomic_rank, t.ncbi_tax_id
                HAVING COUNT(DISTINCT a.taxon_id) > 1
            ) d
        """,
        # The three rank flags on v_abundance_genus must stay consistent:
        # coarser and finer are mutually exclusive, and not_a_genus is exactly
        # their union. They are three hand-written rank lists over one fact, so
        # an edit to one that misses another would silently disagree -- a feature
        # reported as neither coarser nor finer yet not a genus, or as both.
        "genus_rank_flags_inconsistent": """
            SELECT COUNT(*) FROM v_abundance_genus
            WHERE (coarser_than_genus = 1 AND finer_than_genus = 1)
               OR (not_a_genus <> (coarser_than_genus OR finer_than_genus))
        """,
        # effect_type records what the stored number IS, not how it was
        # produced, and resolve_effect_type() exists to make that conversion.
        # The mirror held 123 rows reading 'LEfSe' -- a method -- because its
        # build never called the resolver. MySQL had none, so nothing noticed.
        # The method names are listed in gutdb/transform.py; this is the
        # case-folded form of that set.
        "effect_type_names_a_method": """
            SELECT COUNT(*) FROM taxon_disease_associations
            WHERE LOWER(TRIM(effect_type)) IN
                  ('lefse', 'lefse analysis', 'linear discriminant analysis')
        """,
        # A comparison with one arm recorded and not the other is half a fact:
        # the sign order is only usable as a pair.
        "comparisons_with_only_one_arm": """
            SELECT COUNT(*) FROM phenotype_comparisons
            WHERE (positive_score_enriched_in_id IS NULL)
               <> (negative_score_enriched_in_id IS NULL)
        """,
        # v_taxon_specificity.depleted_fraction had the same contrast-mixing
        # defect and was fixed the same way; a fraction outside [0, 1] means the
        # denominator has drifted from the numerator again.
        "specificity_depleted_fraction_out_of_range": """
            SELECT COUNT(*) FROM v_taxon_specificity
            WHERE depleted_fraction IS NOT NULL
              AND (depleted_fraction < 0 OR depleted_fraction > 1)
        """,
    }


# Checks that cannot run against the SQLite mirror, because it has no
# sample_attributes table and so none of the alias machinery built on it.
_MYSQL_ONLY_CHECKS = frozenset({
    "v_sample_attributes_canonical_row_drift",
    "v_sample_attribute_keys_row_drift",
    "attribute_booleans_undecoded",
})


def validate_sqlite_mirror(sqlite_path: str | None = None) -> list[tuple[str, int]]:
    """Run the portable checks against the SQLite mirror.

    Nothing validated the mirror at all until now. validate_database is typed to
    MySQL, so the mirror's 8 views had no invariants while MySQL's had twelve --
    which is the defect shape this whole schema keeps producing: a guard that
    reached one store only. The (SeqCode) annotations, body_site, the trait
    vocabularies and studies.data_type were all the same mistake in data; this
    was it in the checking.

    The queries are reused rather than rewritten, so the two stores cannot be
    checked against different definitions. Labels are prefixed `sqlite_`, which
    means the pre-push hook picks them up with no change: its gate now treats
    every counter as a defect unless it is named a census figure, so a mirror
    counter is enforced the moment it appears.

    Returns an empty list when the mirror is absent -- it is a build artifact, so
    a checkout without one is normal and must not fail the push.

    THIS PASS IS THE BULK OF validate's RUNTIME: 36.8s of 65.2s, against 28.4s
    for MySQL, mostly the sum-to-1 check over the mirror's 2.05M abundance rows
    and its own copy of v_taxon_disease_evidence. Validate went from roughly 32s
    to 65s when this was added -- the cost of checking the second store at all,
    not of any one counter. Worth knowing before adding more: the pre-push hook
    pays it on every push.
    """
    import os
    import sqlite3
    from pathlib import Path

    if sqlite_path is None:
        sqlite_path = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")
    if not Path(sqlite_path).exists():
        return []
    checks = {label: sql for label, sql in _validate_checks().items()
              if isinstance(label, tuple) or label not in _MYSQL_ONLY_CHECKS}
    connection = sqlite3.connect(sqlite_path)
    try:
        cursor = connection.cursor()
        results = []
        for label, sql in checks.items():
            names = label if isinstance(label, tuple) else (label,)
            try:
                cursor.execute(sql)
                row = cursor.fetchone()
                if len(row) != len(names):
                    raise sqlite3.DataError(
                        f"{len(row)} columns for {len(names)} labels")
                results.extend((f"sqlite_{name}", int(value or 0))
                               for name, value in zip(names, row))
            except sqlite3.Error as error:
                # A check that cannot run is reported as a failure rather than
                # skipped: a silently absent counter is how the mirror went
                # unchecked in the first place.
                results.append((f"sqlite_{names[0]}_FAILED_TO_RUN", 1))
                results.append(
                    (f"sqlite_{names[0]}_error_{type(error).__name__}", 1))
        return results
    finally:
        connection.close()
