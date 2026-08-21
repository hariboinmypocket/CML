from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from .config import Settings
from .db import initialize_database, transaction
from .gmrepo import GMrepoClient
from .cmd_curation import CMDCurationClient
from .mgnify import MGnifyClient
from .ncbi import SRAClient, TaxonomyClient
from .pipeline import (
    LoadStats,
    enrich_sample_demographics,
    clear_contradictory_lineages,
    enrich_cmd_demographics,
    enrich_taxonomy,
    normalize_existing_phyla,
    load_associations,
    load_abundances,
    load_baseline_taxa,
    load_mimedb,
    load_samples,
    sync_gmrepo_phenotypes,
    sync_gmrepo_comparison,
    sync_all_gmrepo_comparisons,
    sync_gmrepo_samples,
    sync_gmrepo_abundances,
    sync_legacy_microbiome_tables,
    sync_mgnify_study,
    validate_database,
)


ROOT = Path(__file__).resolve().parent.parent


def print_stats(label: str, stats: LoadStats) -> None:
    print(
        f"{label}: read={stats.read:,}, inserted={stats.inserted:,}, "
        f"updated={stats.updated:,}, skipped={stats.skipped:,}"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Populate and maintain the normalized gut microbiome database."
    )
    parser.add_argument("--env-file", default=".env", help="Database settings file")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("init", help="Create the database and schema")

    baseline = sub.add_parser("load-baseline", help="Load the existing gut taxa CSV")
    baseline.add_argument("path", nargs="?", default="gut_microbiome_new_full.csv")
    baseline.add_argument("--overwrite", action="store_true")

    mimedb = sub.add_parser("load-mimedb", help="Enrich taxa from MiMeDB")
    mimedb.add_argument("path", nargs="?", default="mimedb_microbes_v2.csv")
    mimedb.add_argument("--overwrite", action="store_true")

    assoc = sub.add_parser("load-associations", help="Load a GMrepo marker TSV/CSV")
    assoc.add_argument("path")
    assoc.add_argument("--fill-only", action="store_true", help="Do not replace existing statistics")

    samples = sub.add_parser("load-samples", help="Load GMrepo clinical/run metadata CSV or TSV")
    samples.add_argument("path")

    abundances = sub.add_parser("load-abundances", help="Load long-format sample taxon abundances")
    abundances.add_argument("path")

    sub.add_parser("sync-phenotypes", help="Load GMrepo phenotype names and MeSH IDs via API")

    comparison = sub.add_parser(
        "sync-comparison", help="Fetch and load a curated GMrepo phenotype comparison"
    )
    comparison.add_argument("--mesh-a", required=True)
    comparison.add_argument("--mesh-b", required=True)
    comparison.add_argument("--output", required=True)
    comparison.add_argument("--exclude-project", action="append", default=[])
    comparison.add_argument("--fill-only", action="store_true")

    all_comparisons = sub.add_parser(
        "sync-all-comparisons",
        help="Fetch every curated GMrepo disease comparison and synchronize legacy tables",
    )
    all_comparisons.add_argument("--output", default="data/gmrepo_all_disease_associations.csv")
    all_comparisons.add_argument("--workers", type=int, default=6)
    all_comparisons.add_argument("--include-health", action="store_true")
    all_comparisons.add_argument("--fill-only", action="store_true")

    sample_sync = sub.add_parser(
        "sync-gmrepo-samples",
        help="Fetch GMrepo run/sample metadata for every imported study project",
    )
    sample_sync.add_argument("--output", default="data/gmrepo_project_samples.csv")
    sample_sync.add_argument("--limit", type=int, default=1000)

    sub.add_parser(
        "sync-legacy", help="Synchronize normalized taxa and disease evidence into microbiome_dataset"
    )

    gm_abund = sub.add_parser(
        "sync-gmrepo-abundances",
        help="Fetch per-run taxonomic profiles from GMrepo into sample_taxon_abundances",
    )
    gm_abund.add_argument(
        "--rank", default="genus", choices=["species", "genus"],
        help="Single taxonomic rank to load per run (loading several would break the sum-to-1 invariant)",
    )
    gm_abund.add_argument("--limit", type=int, default=None, help="Only process this many un-profiled runs")
    gm_abund.add_argument("--workers", type=int, default=4)
    gm_abund.add_argument(
        "--only-complete-demographics", action="store_true",
        help="Restrict to runs that already have sex and age (the immediately ML-usable cohort)",
    )
    gm_abund.add_argument(
        "--data-type", default="",
        help="Restrict to studies of this data_type (e.g. mNGS); GMrepo returns no species for 16S",
    )
    gm_abund.add_argument("--output", default="data/gmrepo_run_abundances.csv")

    mgnify = sub.add_parser(
        "sync-mgnify-study",
        help="Pull a completed MGnify study's sample metadata and precomputed genus/species taxonomy",
    )
    mgnify.add_argument("--study", required=True, help="MGnify study accession, e.g. MGYS00001188")
    mgnify.add_argument("--title", required=True, help="Project title to store on the study row")
    mgnify.add_argument(
        "--stage-map", required=True,
        help="Path to a JSON file mapping a sample metadata disease-stage value to [disease_name, mesh_id]",
    )
    mgnify.add_argument(
        "--body-site", default="",
        help="Fallback body site if a sample's own 'environment (feature)' metadata is missing",
    )
    mgnify.add_argument("--samples-output", default="data/mgnify_samples.csv")
    mgnify.add_argument("--abundances-output", default="data/mgnify_abundances.csv")
    mgnify.add_argument("--sleep", type=float, default=0.3, help="Seconds to sleep between MGnify requests")

    tax = sub.add_parser(
        "enrich-taxonomy",
        help="Fill taxonomic lineage and missing ncbi_tax_id from NCBI Taxonomy (fills NULLs only)",
    )
    tax.add_argument("--limit", type=int, default=None)
    tax.add_argument("--batch-size", type=int, default=150)
    tax.add_argument("--sleep", type=float, default=0.4)
    tax.add_argument("--no-name-resolution", action="store_true",
                     help="Only use taxa that already have an ncbi_tax_id; skip name lookups")
    tax.add_argument("--email", default=None)
    tax.add_argument("--api-key", default=None)

    cmd_dem = sub.add_parser(
        "enrich-cmd-demographics",
        help="Fill sample sex/age/BMI from curatedMetagenomicData's curated metadata (fills NULLs only)",
    )
    cmd_dem.add_argument("--output", default="data/cmd_demographics_audit.csv")

    enrich = sub.add_parser(
        "enrich-demographics",
        help="Backfill sample sex/age/BMI from NCBI BioSample attributes (fills NULLs only)",
    )
    enrich.add_argument(
        "--limit", type=int, default=None,
        help="Only process this many samples missing sex/age/bmi (omit to run against all of them)",
    )
    enrich.add_argument("--batch-size", type=int, default=100, help="Run accessions per NCBI request")
    enrich.add_argument("--sleep", type=float, default=0.4, help="Seconds to sleep between NCBI requests")
    enrich.add_argument(
        "--output", default="data/ncbi_biosample_demographics.csv",
        help="Audit CSV of resolved values, including rejected/unparseable attempts",
    )
    enrich.add_argument("--email", default=None, help="Contact email sent to NCBI E-utilities (recommended)")
    enrich.add_argument("--api-key", default=None, help="NCBI API key for a higher rate limit (10 req/sec)")

    all_cmd = sub.add_parser("all", help="Initialize and run the standard local-data pipeline")
    all_cmd.add_argument("--baseline", default="gut_microbiome_new_full.csv")
    all_cmd.add_argument("--mimedb", default="mimedb_microbes_v2.csv")
    all_cmd.add_argument("--associations", action="append", default=[])
    all_cmd.add_argument("--sync-phenotypes", action="store_true")
    all_cmd.add_argument("--overwrite", action="store_true")

    sub.add_parser("validate", help="Print row counts and referential-integrity checks")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    settings = Settings.from_env(args.env_file)
    schema_path = Path(__file__).with_name("schema.sql")

    try:
        if args.command in {"init", "all"}:
            initialize_database(settings, schema_path)
            print(f"Initialized MySQL database: {settings.database}")
            if args.command == "init":
                return 0

        with transaction(settings) as connection:
            if args.command == "load-baseline":
                print_stats("Baseline", load_baseline_taxa(connection, args.path, args.overwrite))
            elif args.command == "load-mimedb":
                print_stats("MiMeDB", load_mimedb(connection, args.path, args.overwrite))
            elif args.command == "load-associations":
                print_stats("Associations", load_associations(connection, args.path, not args.fill_only))
            elif args.command == "load-samples":
                print_stats("Samples", load_samples(connection, args.path))
            elif args.command == "load-abundances":
                print_stats("Abundances", load_abundances(connection, args.path))
            elif args.command == "sync-phenotypes":
                client = GMrepoClient(settings.gmrepo_base_url)
                print_stats("GMrepo phenotypes", sync_gmrepo_phenotypes(connection, client))
            elif args.command == "sync-comparison":
                client = GMrepoClient(settings.gmrepo_base_url)
                stats, projects = sync_gmrepo_comparison(
                    connection,
                    client,
                    args.mesh_a,
                    args.mesh_b,
                    args.output,
                    args.exclude_project,
                    not args.fill_only,
                )
                print(f"GMrepo projects: {', '.join(projects)}")
                print_stats("GMrepo comparison", stats)
            elif args.command == "sync-all-comparisons":
                client = GMrepoClient(settings.gmrepo_base_url)
                stats, metadata = sync_all_gmrepo_comparisons(
                    connection,
                    client,
                    args.output,
                    args.workers,
                    args.include_health,
                    not args.fill_only,
                )
                print(
                    "GMrepo catalog: "
                    f"comparisons={metadata['comparisons']:,}, "
                    f"projects={metadata['projects']:,}, "
                    f"exported_rows={metadata['exported_rows']:,}, "
                    f"pruned_rows={metadata['pruned_rows']:,}"
                )
                print_stats("GMrepo all comparisons", stats)
                legacy = sync_legacy_microbiome_tables(connection)
                print(
                    "Legacy sync: "
                    f"taxa_before={legacy['legacy_taxa_before']:,}, "
                    f"taxa_inserted={legacy['legacy_taxa_inserted']:,}, "
                    f"taxa_after={legacy['legacy_taxa_after']:,}, "
                    f"associations={legacy['legacy_associations']:,}"
                )
            elif args.command == "sync-gmrepo-samples":
                client = GMrepoClient(settings.gmrepo_base_url)
                stats, metadata = sync_gmrepo_samples(
                    connection,
                    client,
                    args.output,
                    args.limit,
                )
                print(
                    "GMrepo samples: "
                    f"projects={metadata['projects']:,}, "
                    f"exported_rows={metadata['exported_rows']:,}"
                )
                print_stats("Samples", stats)
            elif args.command == "sync-gmrepo-abundances":
                client = GMrepoClient(settings.gmrepo_base_url)
                stats, meta = sync_gmrepo_abundances(
                    connection,
                    client,
                    args.output,
                    args.rank,
                    args.limit,
                    args.workers,
                    args.only_complete_demographics,
                    args.data_type,
                )
                print(
                    "GMrepo profiles: "
                    f"requested={meta['runs_requested']:,}, "
                    f"with_profile={meta['runs_with_profile']:,}, "
                    f"empty={meta['runs_empty']:,}, "
                    f"failed={meta['runs_failed']:,}, "
                    f"genus_aggregated={meta['runs_aggregated_to_genus']:,}, "
                    f"rows_exported={meta['rows_exported']:,}"
                )
                print_stats("GMrepo abundances", stats)
            elif args.command == "sync-mgnify-study":
                with open(args.stage_map, encoding="utf-8") as handle:
                    raw_map = json.load(handle)
                stage_map = {k.casefold(): tuple(v) for k, v in raw_map.items()}
                client = MGnifyClient()
                sample_stats, abundance_stats, metadata = sync_mgnify_study(
                    connection,
                    client,
                    args.study,
                    args.title,
                    stage_map,
                    args.body_site,
                    args.samples_output,
                    args.abundances_output,
                    args.sleep,
                )
                print(
                    "MGnify study: "
                    f"samples_seen={metadata['samples_seen']:,}, "
                    f"samples_no_analysis={metadata['samples_no_analysis']:,}, "
                    f"taxa_rows={metadata['taxa_rows']:,}"
                )
                print_stats("MGnify samples", sample_stats)
                print_stats("MGnify abundances", abundance_stats)
            elif args.command == "enrich-taxonomy":
                folded = normalize_existing_phyla(connection)
                print(f"Folded {folded:,} existing rows onto classic phylum names")
                cleared = clear_contradictory_lineages(connection)
                print(f"Cleared {cleared:,} lineages contradicting their own superkingdom")
                client = TaxonomyClient(email=args.email, api_key=args.api_key)
                stats, counters = enrich_taxonomy(
                    connection, client, args.limit, args.batch_size,
                    args.sleep, not args.no_name_resolution,
                )
                print_stats("NCBI Taxonomy", stats)
                print(
                    f"  had_taxid={counters['with_taxid']:,}, "
                    f"names_resolved={counters['names_resolved']:,}, "
                    f"names_unresolved={counters['names_unresolved']:,}, "
                    f"lineage_filled={counters['lineage_filled']:,}, "
                    f"taxid_filled={counters['taxid_filled']:,}, "
                    f"no_lineage={counters['no_lineage_returned']:,}"
                )
            elif args.command == "enrich-cmd-demographics":
                stats, counters = enrich_cmd_demographics(
                    connection, CMDCurationClient(), args.output)
                print_stats("cMD demographics", stats)
                print(
                    f"  studies={counters['cmd_studies']}, matched={counters['matched']:,} | "
                    f"filled sex={counters['sex_filled']:,} age={counters['age_filled']:,} "
                    f"bmi={counters['bmi_filled']:,} | conflicts sex={counters['sex_conflicts']:,} "
                    f"age={counters['age_conflicts']:,} | age_rejected={counters['age_rejected']:,}"
                )
            elif args.command == "enrich-demographics":
                client = SRAClient(email=args.email, api_key=args.api_key)
                stats, counters = enrich_sample_demographics(
                    connection,
                    client,
                    limit=args.limit,
                    batch_size=args.batch_size,
                    sleep_seconds=args.sleep,
                    audit_csv=args.output,
                )
                print_stats("NCBI demographics", stats)
                print(
                    "Filled: "
                    f"sex={counters['sex_filled']:,}, age={counters['age_filled']:,}, "
                    f"bmi={counters['bmi_filled']:,}; "
                    f"no_biosample_attrs={counters['no_biosample_attrs']:,}"
                )
                print(f"Audit CSV: {args.output}")
            elif args.command == "sync-legacy":
                legacy = sync_legacy_microbiome_tables(connection)
                print(
                    "Legacy sync: "
                    f"taxa_before={legacy['legacy_taxa_before']:,}, "
                    f"taxa_inserted={legacy['legacy_taxa_inserted']:,}, "
                    f"taxa_after={legacy['legacy_taxa_after']:,}, "
                    f"associations={legacy['legacy_associations']:,}"
                )
            elif args.command == "all":
                print_stats("Baseline", load_baseline_taxa(connection, args.baseline, args.overwrite))
                print_stats("MiMeDB", load_mimedb(connection, args.mimedb, args.overwrite))
                if args.sync_phenotypes:
                    client = GMrepoClient(settings.gmrepo_base_url)
                    print_stats("GMrepo phenotypes", sync_gmrepo_phenotypes(connection, client))
                for path in args.associations:
                    print_stats(path, load_associations(connection, path, args.overwrite))
            elif args.command == "validate":
                for label, count in validate_database(connection):
                    print(f"{label}: {count:,}")
        return 0
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
