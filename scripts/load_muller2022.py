#!/usr/bin/env python3
"""Load the Muller 2022 microbiome-metabolome collection.

Muller, Algavi & Borenstein 2022, npj Biofilms and Microbiomes 8:79
(doi:10.1038/s41522-022-00345-5). Fetch the files first with
scripts/fetch_muller2022.py and review data/gtdb_ncbi_genus_map.csv, which
scripts/map_gtdb_genera.py writes and this loader consumes rather than
re-deciding.

This is the database's second sample-level source. 98.3% of samples came from
GMrepo and MicrobiomeHD contributes 20 studies but no samples, so the point is
less the sample count than having a second cohort-level opinion -- plus a
metabolome, which nothing here has carried before.

Four kinds of double-counting were found in the data and are all refused:

  iHMP_IBDMDB_2019 is already loaded as PRJNA398089 from GMrepo. The whole
  cohort is skipped; its 382 samples are not new.

  53 samples appear in BOTH YACHIDA_CRC_2019 and ERAWIJANTARI_GASTRIC_CANCER_2020
  and the collection flags them itself. All 53 are Healthy controls. They are
  dropped from YACHIDA, which keeps 74 controls and a working contrast; dropping
  them from ERAWIJANTARI instead would leave it with one. Their subject ids
  differ only by a suffix ('10025' against '10025.Healthy'), so a dedup keyed on
  subject_id would not have caught this -- the suffix is stripped on load so the
  two cohorts agree about who these people are.

  Several GTDB labels map to one NCBI genus, so abundances are SUMMED per
  target. Clostridium receives 40 GTDB fragments, Enterococcus 11, Ruminococcus
  7, Bacteroides 6, and 52% of all abundance mass sits in such targets. An
  insert-per-label would have let the last write win and silently discarded half
  the data, which is what INSERT OR IGNORE did to the SQLite mirror.

  Samples carry one row per (sample, taxon), so repeated runs of one subject stay
  distinguishable through subject_id rather than being merged.

Two cohort arms are loaded WITHOUT a disease label, pending a judgement this
script should not make:

  YACHIDA_CRC_2019 'HS' (30 samples) is a group code the metadata does not
  decode and the paper is needed for.

  ERAWIJANTARI 'Gastrectomy' (42 samples) is a post-resection state, not an
  active tumour -- Surgery_Type confirms subtotal and total gastrectomies -- so
  labelling it Stomach Neoplasms would assert something false. The 54 controls
  load normally.

Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect

csv.field_size_limit(10 ** 9)

DRY_RUN = "--apply" not in sys.argv
SOURCE = Path("data/muller2022")
GENUS_MAP = Path("data/gtdb_ncbi_genus_map.csv")
SOURCE_DB = "MullerMM"

# Already here as PRJNA398089 from GMrepo.
SKIP_COHORTS = {"iHMP_IBDMDB_2019"}

# Drop rows this cohort shares with another cohort in the same collection.
SHARED_DROP = {"YACHIDA_CRC_2019": "Shared.w.ERAWIJANTARI_2020"}

# Genus labels below this share of total collection mass are not worth a taxa
# row: mapping everything would add 1,853 genera beginning Abyssibacter,
# Acaryochloris and Acetohalobium, which are marine and environmental kraken2
# noise. This floor keeps 395 genera and 99.81% of mappable mass for 126 new
# rows.
MASS_FLOOR = 0.00001

# 16S resolves to genus at best and shotgun sprays reads across far more
# lineages, so genus richness separates the two cleanly here: the 16S cohorts
# hold 82-584 genus labels and the shotgun cohorts 2,895-11,942. Taken from the
# data because the collection does not ship a per-cohort platform column, and
# YACHIDA (11,942 labels, unmistakably shotgun) ships no species table, so file
# presence is not a usable signal.
MNGS_GENUS_THRESHOLD = 1000

SKIP = "__skip__"           # do not load the sample at all
UNLABELLED = None           # load the sample, leave disease_id NULL

# Study.Group -> disease name in `diseases`. Hand-written per cohort: these are
# free-text codes, several undecodable without the original paper, and guessing
# one mislabels every sample in an arm.
PHENOTYPES: dict[str, dict[str, str | None]] = {
    "SINHA_CRC_2016": {"1": "Colorectal Neoplasms", "0": "Health"},
    "MARS_IBS_2020": {
        # IBS-D and IBS-C are subtypes of one MeSH disease; the schema has no
        # subtype field, so both map to IBS and the distinction is lost here.
        "D": "Irritable Bowel Syndrome",
        "C": "Irritable Bowel Syndrome",
        "H": "Health",
    },
    "YACHIDA_CRC_2019": {
        "Healthy": "Health",
        "Stage_0": "Colorectal Neoplasms",
        "Stage_I_II": "Colorectal Neoplasms",
        "Stage_III_IV": "Colorectal Neoplasms",
        "MP": "Adenoma",
        "HS": UNLABELLED,
    },
    "ERAWIJANTARI_GASTRIC_CANCER_2020": {
        "Healthy": "Health",
        "Gastrectomy": UNLABELLED,
    },
    "FRANZOSA_IBD_2019": {
        "CD": "Crohn Disease", "UC": "Colitis, Ulcerative", "Control": "Health"},
    "JACOBS_IBD_FAMILIES_2016": {
        "CD": "Crohn Disease", "UC": "Colitis, Ulcerative", "Normal": "Health"},
    "KANG_AUTISM_2017": {
        "Autistic": "Autism Spectrum Disorder", "Neurotypical": "Health"},
    "KIM_ADENOMAS_2020": {
        "Adenoma": "Adenoma", "Carcinoma": "Colorectal Neoplasms", "Control": "Health"},
    "KOSTIC_INFANTS_DIABETES_2015": {
        "case": "Diabetes Mellitus, Type 1", "control": "Health"},
    "WANDRO_PRETERMS_2018": {
        "control": "Health", "septic": "Sepsis",
        "nec": "Enterocolitis, Necrotizing",
        "GI": SKIP,                      # n=1, no usable arm
    },
    "WANG_ESRD_2020": {"ESRD": "Kidney Failure, Chronic", "Control": "Health"},
    # Healthy infants on different diets, and healthy stool-bank donors. The
    # groups are diet or timepoint rather than disease, so every sample is a
    # control. Age is recorded, which is what makes infant cohorts usable
    # alongside adult ones rather than a hidden batch effect.
    "HE_INFANTS_MFGM_2019": {
        "Baseline": "Health", "Month12": "Health",
        "With.comp.food": "Health", "Without.comp.food": "Health"},
    "POYET_BIO_ML_2019": {"": "Health"},
}

# MeSH ids for rows this loader may have to create.
NEW_DISEASES = {"Adenoma": "D000236"}

HMDB_RE = re.compile(r"^HMDB\d{5,7}$")
KEGG_RE = re.compile(r"^C\d{5}$")
RUN_ACC_RE = re.compile(r"^(SRR|ERR|DRR)\w+$")
NA = {"", "NA", "na", "NaN", "None", "-", "nan"}


def norm_hmdb(value: str) -> str:
    """HMDB00002 and HMDB0000002 are one compound; store the 7-digit form."""
    return "HMDB" + value[4:].zfill(7)


def clean(value: str | None) -> str:
    return (value or "").strip()


def load_genus_map() -> dict[str, str]:
    """GTDB lineage label -> NCBI genus name, for labels above the mass floor."""
    if not GENUS_MAP.exists():
        raise SystemExit(f"{GENUS_MAP} missing -- run scripts/map_gtdb_genera.py first")
    keep: dict[str, str] = {}
    by_target: dict[str, float] = defaultdict(float)
    rows = list(csv.DictReader(GENUS_MAP.open(newline="", encoding="utf-8")))
    for row in rows:
        if row["decision"].startswith("map"):
            by_target[row["base_genus"]] += float(row["mass_share"])
    for row in rows:
        if row["decision"].startswith("map") and by_target[row["base_genus"]] >= MASS_FLOOR:
            keep[row["gtdb_label"]] = row["base_genus"]
    return keep


def parse_float(value: str | None) -> float | None:
    v = clean(value)
    if v in NA:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def parse_sex(value: str | None) -> str | None:
    v = clean(value).lower()
    if v in ("male", "m"):
        return "male"
    if v in ("female", "f"):
        return "female"
    return None


def age_years(row: dict) -> float | None:
    """Age in years, honouring Age.Units -- infant cohorts report months or days."""
    raw = parse_float(row.get("Age"))
    if raw is None:
        return None
    unit = clean(row.get("Age.Units")).lower()
    if unit in ("years", "year", "y", ""):
        return raw
    if unit in ("months", "month"):
        return raw / 12.0
    if unit in ("weeks", "week"):
        return raw / 52.0
    if unit in ("days", "day"):
        return raw / 365.25
    return None


def subject_id(study: str, row: dict) -> str:
    """Subject, with the cohort-local suffix stripped.

    ERAWIJANTARI writes the 53 samples it shares with YACHIDA as '10025.Healthy'
    where YACHIDA writes '10025'. Keeping both spellings would hide that they
    are the same person from any later cross-cohort check.
    """
    value = clean(row.get("Subject"))
    if study == "ERAWIJANTARI_GASTRIC_CANCER_2020" and "." in value:
        return value.split(".", 1)[0]
    return value


def ensure_diseases(cur) -> dict[str, int]:
    """Resolve every disease name the mappings use, creating only known gaps."""
    wanted = {name for table in PHENOTYPES.values() for name in table.values()
              if name not in (SKIP, UNLABELLED)}
    cur.execute("SELECT name, id FROM diseases")
    have = {name: did for name, did in cur.fetchall()}
    for name in sorted(wanted - set(have)):
        mesh = NEW_DISEASES.get(name)
        if mesh is None:
            raise SystemExit(
                f"disease {name!r} is not in `diseases` and has no MeSH id in "
                f"NEW_DISEASES -- refusing to invent one")
        if not DRY_RUN:
            cur.execute(
                "INSERT INTO diseases (mesh_id, name, name_key) VALUES (%s, %s, %s)",
                (mesh, name, name.casefold()))
            have[name] = cur.lastrowid
        else:
            have[name] = -1
        print(f"   disease created: {name} ({mesh})")
    return {n: have[n] for n in wanted}


def ensure_genera(cur, needed: set[str], taxdump) -> dict[str, int]:
    """Resolve each NCBI genus name to a taxa row, creating the ones missing.

    New rows get their lineage from the local NCBI dump, not from GTDB: GTDB
    writes Firmicutes_A where NCBI writes Firmicutes, and normalize_phylum and
    the phylum columns elsewhere in this schema are NCBI-shaped. The domain the
    GTDB label declares is used only to pick between homonyms, which is the
    guard Bacillus-the-stick-insect needs.
    """
    cur.execute("SELECT LOWER(scientific_name), id FROM taxa WHERE taxonomic_rank = 'genus'")
    have = {name: tid for name, tid in cur.fetchall()}
    resolved = {g: have[g.lower()] for g in needed if g.lower() in have}
    missing = sorted(g for g in needed if g.lower() not in have)
    if not missing:
        return resolved

    name_to_taxid, lineages = taxdump
    created = skipped = 0
    for genus in missing:
        candidates = name_to_taxid.get(genus.casefold(), [])
        pick = None
        for taxid in candidates:
            lineage = lineages.get(taxid)
            top = lineage.get("superkingdom") or lineage.get("domain") or ""
            if top in ("Bacteria", "Archaea"):
                pick = (taxid, lineage)
                break
        if pick is None:
            skipped += 1
            continue
        taxid, lineage = pick
        if not DRY_RUN:
            cur.execute(
                """INSERT INTO taxa (superkingdom, phylum, class_name, order_name,
                                     family, genus, species, genus_key, species_key,
                                     scientific_name, taxonomic_rank, ncbi_tax_id,
                                     source_database)
                   VALUES (%s,%s,%s,%s,%s,%s,'',%s,'',%s,'genus',%s,%s)""",
                (lineage.get("superkingdom") or lineage.get("domain"),
                 lineage.get("phylum"), lineage.get("class"), lineage.get("order"),
                 lineage.get("family"), genus, genus.casefold(), genus, taxid, SOURCE_DB))
            resolved[genus] = cur.lastrowid
        else:
            resolved[genus] = -1
        created += 1
    print(f"   genus rows: {len(resolved) - created} existing, {created} created"
          + (f", {skipped} unresolvable in the dump and skipped" if skipped else ""))
    return resolved


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_abundances(path: Path, keep: dict[str, str]) -> tuple[dict[str, dict[str, float]], float]:
    """sample -> {NCBI genus: relative abundance}, summed per target then renormalized.

    Returns the profiles and the mean share of each sample's mass that survived
    the map, so the caller can report what renormalizing had to absorb.
    """
    profiles: dict[str, dict[str, float]] = {}
    retained: list[float] = []
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)[1:]
        targets = [keep.get(col) for col in header]
        for row in reader:
            sample = row[0]
            summed: dict[str, float] = defaultdict(float)
            total = 0.0
            for target, value in zip(targets, row[1:]):
                if value in ("", "NA"):
                    continue
                x = float(value)
                total += x
                if target:
                    summed[target] += x
            kept = sum(summed.values())
            if kept <= 0:
                continue
            retained.append(kept / total if total else 0.0)
            # Restore the per-rank sum-to-1 invariant after the rejected labels
            # (MAG bins, above-genus lineages) are dropped.
            #
            # Zeros are dropped, not stored. The GTDB tables are dense and write
            # 0.0 for every genus a sample lacks, which would have added 233,965
            # rows meaning "measured as absent" to a table where absence is
            # represented by the ABSENCE of a row. That distinction is load
            # bearing: scripts/adjudicate_conflicts.py counts a missing row as
            # zero and divides by the number of samples profiled at the rank, so
            # explicit zeros would not change an arithmetic mean but would make
            # every prevalence or detection count wrong. Dropping them cannot
            # disturb the sum.
            profiles[sample] = {g: v / kept for g, v in summed.items() if v > 0}
    return profiles, (sum(retained) / len(retained) if retained else 0.0)


def load_metabolites(study_dir: Path) -> tuple[dict, dict, dict]:
    """sample -> {compound column: level}, and compound column -> identifiers.

    Only identifier-bearing compounds are kept: the rest are unannotated m/z
    peaks with no identity shared across cohorts, so they cannot be compared
    with anything. Identifier columns are validated rather than trusted --
    upstream carries MetaboAnalyst METPA ids and KEGG C##### values in the HMDB
    column, and KEGG DRUG D##### ids in the KEGG column.
    """
    mtb = study_dir / "mtb.tsv"
    mapping = study_dir / "mtb.map.tsv"
    if not (mtb.exists() and mapping.exists()):
        return {}, {}, {}
    compounds: dict[str, dict] = {}
    for row in read_rows(mapping):
        h, k = clean(row.get("HMDB")), clean(row.get("KEGG"))
        hmdb = norm_hmdb(h) if HMDB_RE.match(h) else None
        kegg = k if KEGG_RE.match(k) else None
        if not (hmdb or kegg):
            continue
        compounds[row["Compound"]] = {
            "hmdb_id": hmdb,
            "kegg_id": kegg,
            "identity": hmdb or f"kegg:{kegg}",
            "name": clean(row.get("Compound.Name")) or row["Compound"],
            "high_confidence": clean(row.get("High.Confidence.Annotation")).upper() != "FALSE",
        }

    # Several measured columns in one cohort can carry the same identifier, and
    # (sample_id, metabolite_id) is a primary key, so something has to give. The
    # two causes are indistinguishable from the data:
    #
    #   the same analyte measured twice -- HILIC_NEG_alanine beside
    #   HILIC_POS_alanine, or 'cholate' beside 'Cholic acid'
    #
    #   an upstream annotation error -- ERAWIJANTARI gives both C00197_3PG and
    #   C00661_G3P the identifier HMDB0000807, but 3-phosphoglycerate and
    #   glycerol-3-phosphate are different compounds
    #
    # Averaging would invent a number across two different molecules in the
    # second case, and letting the last insert win would discard a measurement
    # silently in both. The colliding identities are therefore dropped from this
    # cohort and reported. 124 of 17,066 columns are affected, and the compounds
    # remain available from the cohorts where they do not collide.
    seen: dict[str, list[str]] = defaultdict(list)
    for col, meta in compounds.items():
        seen[meta["identity"]].append(col)
    collided = {ident: cols for ident, cols in seen.items() if len(cols) > 1}
    for cols in collided.values():
        for col in cols:
            compounds.pop(col, None)
    levels: dict[str, dict[str, float]] = {}
    with mtb.open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)[1:]
        for row in reader:
            per_sample = {}
            for col, value in zip(header, row[1:]):
                if col in compounds:
                    x = parse_float(value)
                    if x is not None:
                        per_sample[col] = x
            if per_sample:
                levels[row[0]] = per_sample
    return levels, compounds, collided


def ensure_metabolite_rows(cur, compounds: dict[str, dict], cache: dict) -> dict[str, int]:
    """Resolve each compound to a metabolites row, keyed on HMDB then KEGG."""
    out: dict[str, int] = {}
    for col, meta in compounds.items():
        key = meta["hmdb_id"] or f"kegg:{meta['kegg_id']}"
        if key in cache:
            out[col] = cache[key]
            continue
        row = None
        if meta["hmdb_id"]:
            cur.execute("SELECT id FROM metabolites WHERE hmdb_id = %s", (meta["hmdb_id"],))
            row = cur.fetchone()
        if row is None and meta["kegg_id"]:
            cur.execute("SELECT id FROM metabolites WHERE kegg_id = %s", (meta["kegg_id"],))
            row = cur.fetchone()
        if row is not None:
            cache[key] = out[col] = row[0]
            continue
        if DRY_RUN:
            cache[key] = out[col] = -1
            continue
        cur.execute(
            "INSERT INTO metabolites (hmdb_id, kegg_id, name) VALUES (%s, %s, %s)",
            (meta["hmdb_id"], meta["kegg_id"], meta["name"][:255]))
        cache[key] = out[col] = cur.lastrowid
    return out


def main() -> int:
    keep = load_genus_map()
    targets = sorted(set(keep.values()))
    print(f"genus map: {len(keep):,} GTDB labels -> {len(targets)} NCBI genera "
          f"(mass floor {100*MASS_FLOOR:g}%)")

    cohorts = sorted(d for d in SOURCE.glob("*") if d.is_dir())
    cohorts = [d for d in cohorts if d.name not in SKIP_COHORTS]
    print(f"cohorts to load: {len(cohorts)} "
          f"({', '.join(sorted(SKIP_COHORTS))} skipped as already present)\n")

    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()

    diseases = ensure_diseases(cur)

    taxdump = None
    missing_genera = set()
    cur.execute("SELECT LOWER(scientific_name) FROM taxa WHERE taxonomic_rank = 'genus'")
    present = {r[0] for r in cur.fetchall()}
    missing_genera = {g for g in targets if g.lower() not in present}
    if missing_genera:
        print(f"   {len(missing_genera)} genera need taxa rows; reading the NCBI dump")
        from fix_curated_csv_lineage import load_taxdump
        taxdump = load_taxdump()
    genera = ensure_genera(cur, set(targets), taxdump or ({}, None))

    # Drop labels whose target never got a taxa row BEFORE any renormalizing.
    # Copromorpha is the live case: GTDB uses it for a bacterial genus, and the
    # only NCBI taxid of that name (1181387) is a moth, so ensure_genera's
    # domain guard refuses it -- the same guard the curated-CSV lineage repair
    # exists to enforce. Renormalizing first and filtering afterwards would
    # leave its 0.05% of mass in the denominator with no row to carry it, and
    # every affected sample would then fail the per-rank sum-to-1 invariant.
    unresolved = sorted({t for t in targets if t not in genera})
    if unresolved:
        keep = {label: target for label, target in keep.items() if target in genera}
        print(f"   dropped before renormalizing (no NCBI genus row): "
              f"{', '.join(unresolved)}")

    met_cache: dict[str, int] = {}
    totals = defaultdict(int)
    print(f"\n{'cohort':34} {'samples':>8} {'abund':>9} {'mtb lvls':>9} {'kept%':>7}  platform")
    for study_dir in cohorts:
        study = study_dir.name
        table = PHENOTYPES.get(study)
        if table is None:
            print(f"{study:34} -- no phenotype mapping, skipped")
            continue

        meta = read_rows(study_dir / "metadata.tsv")
        shared_col = SHARED_DROP.get(study)
        profiles, kept_share = load_abundances(study_dir / "genera.tsv", keep)
        levels, compounds, collided = load_metabolites(study_dir)
        n_genera = len(next(csv.reader(
            (study_dir / "genera.tsv").open(newline="", encoding="utf-8"),
            delimiter="\t"))) - 1
        platform = "mNGS" if n_genera >= MNGS_GENUS_THRESHOLD else "16S"

        if not DRY_RUN:
            cur.execute(
                """INSERT INTO studies (project_accession, title, data_type,
                                        source_database, metabolomics_method, data_quality)
                   VALUES (%s,%s,%s,%s,%s,'curated')
                   ON DUPLICATE KEY UPDATE data_type = VALUES(data_type),
                       source_database = VALUES(source_database)""",
                (f"{SOURCE_DB}_{study}", meta[0].get("Publication.Name", study)[:500],
                 platform, SOURCE_DB, "see mtb.map.tsv" if compounds else None))
            cur.execute("SELECT id FROM studies WHERE project_accession = %s",
                        (f"{SOURCE_DB}_{study}",))
            study_id = cur.fetchone()[0]
        else:
            study_id = -1

        # Resolve compound columns to metabolites rows once per cohort, then
        # drop any that land on the SAME row. Checking identifier strings is not
        # enough: one column carrying only an HMDB id and another carrying only
        # a KEGG id can both resolve to a row that holds both, which collides on
        # (sample_id, metabolite_id) just the same and was quietly costing 700
        # levels.
        met_ids: dict[str, int] = {}
        if compounds and not DRY_RUN:
            met_ids = ensure_metabolite_rows(cur, compounds, met_cache)
            per_row: dict[int, list[str]] = defaultdict(list)
            for col, mid in met_ids.items():
                per_row[mid].append(col)
            for mid, cols in per_row.items():
                if len(cols) > 1:
                    totals["collided_metabolites"] += len(cols)
                    for col in cols:
                        met_ids.pop(col, None)

        n_samples = n_abund = n_levels = 0
        skipped_group = dropped_shared = unlabelled = 0
        for row in meta:
            sample = clean(row.get("Sample"))
            if shared_col and clean(row.get(shared_col)).upper() == "TRUE":
                dropped_shared += 1
                continue
            group = clean(row.get("Study.Group"))
            if group not in table:
                skipped_group += 1
                continue
            label = table[group]
            if label == SKIP:
                skipped_group += 1
                continue
            disease_id = diseases.get(label) if label is not UNLABELLED else None
            if label is UNLABELLED:
                unlabelled += 1
            profile = profiles.get(sample)
            if not profile:
                continue

            run = clean(row.get("Run"))
            accession = run if RUN_ACC_RE.match(run) else f"{SOURCE_DB}_{study}_{sample}"
            if not DRY_RUN:
                cur.execute(
                    """INSERT INTO samples (study_id, disease_id, run_accession,
                                            subject_id, sex, age_years, bmi, qc_status)
                       VALUES (%s,%s,%s,%s,%s,%s,%s,'curated')
                       ON DUPLICATE KEY UPDATE disease_id = VALUES(disease_id),
                           subject_id = VALUES(subject_id), sex = VALUES(sex),
                           age_years = VALUES(age_years), bmi = VALUES(bmi)""",
                    (study_id, disease_id, accession[:128], subject_id(study, row)[:128],
                     parse_sex(row.get("Gender")), age_years(row), parse_float(row.get("BMI"))))
                cur.execute("SELECT id FROM samples WHERE run_accession = %s", (accession[:128],))
                sample_id = cur.fetchone()[0]
                rows_ab = [(sample_id, genera[g], v) for g, v in profile.items()
                           if g in genera and genera[g] > 0]
                cur.executemany(
                    """INSERT INTO sample_taxon_abundances (sample_id, taxon_id, relative_abundance)
                       VALUES (%s,%s,%s)
                       ON DUPLICATE KEY UPDATE relative_abundance = VALUES(relative_abundance)""",
                    rows_ab)
                n_abund += len(rows_ab)
                if sample in levels and met_ids:
                    rows_mt = [(sample_id, met_ids[c], v, compounds[c]["high_confidence"])
                               for c, v in levels[sample].items() if met_ids.get(c, 0) > 0]
                    cur.executemany(
                        """INSERT INTO sample_metabolite_levels
                             (sample_id, metabolite_id, level, high_confidence_annotation)
                           VALUES (%s,%s,%s,%s)
                           ON DUPLICATE KEY UPDATE level = VALUES(level)""",
                        rows_mt)
                    n_levels += len(rows_mt)
            else:
                n_abund += sum(1 for g in profile if g in genera)
                n_levels += len(levels.get(sample, {}))   # upper bound: see note
            n_samples += 1

        totals["samples"] += n_samples
        totals["abundances"] += n_abund
        totals["levels"] += n_levels
        totals["unlabelled"] += unlabelled
        totals["dropped_shared"] += dropped_shared
        totals["skipped_group"] += skipped_group
        notes = []
        if dropped_shared:
            notes.append(f"-{dropped_shared} shared")
        if skipped_group:
            notes.append(f"-{skipped_group} group")
        if unlabelled:
            notes.append(f"{unlabelled} unlabelled")
        if collided:
            dropped_cols = sum(len(v) for v in collided.values())
            totals["collided_metabolites"] += dropped_cols
            notes.append(f"-{dropped_cols} mtb id clashes")
        print(f"{study:34} {n_samples:>8} {n_abund:>9} {n_levels:>9} "
              f"{100*kept_share:>6.1f}%  {platform}"
              + (f"  ({', '.join(notes)})" if notes else ""))

    if not DRY_RUN:
        conn.commit()
    print(f"\nsamples {totals['samples']:,}   abundance rows {totals['abundances']:,}   "
          f"metabolite levels {totals['levels']:,}")
    print(f"unlabelled (disease_id NULL): {totals['unlabelled']}   "
          f"dropped as shared: {totals['dropped_shared']}   "
          f"skipped by group: {totals['skipped_group']}")
    print(f"metabolite columns dropped for sharing an identifier within a cohort: "
          f"{totals['collided_metabolites']}")
    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")

    cur.close()
    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
