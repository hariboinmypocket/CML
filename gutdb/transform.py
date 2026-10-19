from __future__ import annotations

from collections import Counter
from decimal import Decimal, ROUND_FLOOR, ROUND_HALF_UP
import math
import re
from typing import Any, Iterable, Mapping


NA_VALUES = {"", "n/a", "na", "none", "null", "nan", "unknown", "not available"}


def verified_demographic_link(record: Mapping[str, Any]) -> bool:
    """A plausible library prefix is not a clinical-to-run identifier match."""
    return all(record.get(key) is True for key in (
        "published_values_match_audit", "ena_library_to_run_verified",
        "exact_source_identifier_in_ena",
    ))


def same_demographic_value(field: str, left: Any, right: Any) -> bool:
    """Compare CSV and SQL values without treating missing values as evidence."""
    a, b = clean_str(left), clean_str(right)
    if not a or not b:
        return False
    if field == "sex":
        return a.casefold() == b.casefold()
    try:
        return abs(Decimal(a) - Decimal(b)) < Decimal("0.000001")
    except (ValueError, ArithmeticError):
        return False


def reconcile_demographic_cells(current: Mapping[str, Any], source: Mapping[str, Any],
                                verified: bool) -> tuple[dict[str, Any], list[str]]:
    """Withdraw matching unsafe imports, fill verified gaps, preserve conflicts.

    Shared by both database stores and CSV repair. Never overwrite an independent
    nonmissing measurement, and never manufacture a source-missing measurement.
    """
    changes, conflicts = {}, []
    for field in ("sex", "age_years", "bmi"):
        present, imported = current.get(field), source.get(field)
        if not clean_str(imported):
            continue
        if not clean_str(present):
            if verified:
                changes[field] = imported
        elif same_demographic_value(field, present, imported):
            if not verified:
                changes[field] = None
        else:
            conflicts.append(field)
    return changes, conflicts


def clean_str(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    text = re.sub(r"\s+", " ", str(value)).strip()
    return "" if text.casefold() in NA_VALUES else text


def first_value(row: Mapping[str, Any], *names: str) -> str:
    folded = {str(k).strip().casefold(): v for k, v in row.items()}
    for name in names:
        value = clean_str(folded.get(name.casefold()))
        if value:
            return value
    return ""


def normalize_genus(value: Any) -> str:
    genus = clean_str(value).strip("[]")
    if not genus:
        return ""
    parts = genus.split()
    if len(parts) >= 2 and parts[0].casefold() == "candidatus":
        return "Candidatus " + parts[1][0].upper() + parts[1][1:].lower()
    return genus[0].upper() + genus[1:].lower()


def extract_species_epithet(genus: Any, species_field: Any) -> str:
    """Return the species portion while retaining `sp.` and strain identifiers."""
    g = normalize_genus(genus)
    species = clean_str(species_field).replace("[", "").replace("]", "")
    if not species:
        return ""
    parts = species.split()
    if g and parts and parts[0].casefold() == g.casefold():
        parts = parts[1:]
    if not parts:
        return ""
    return " ".join(parts).lower()


def parse_scientific_name(value: Any) -> tuple[str, str]:
    name = clean_str(value).replace("[", "").replace("]", "")
    parts = name.split()
    if len(parts) < 2:
        return normalize_genus(parts[0]) if parts else "", ""
    if parts[0].casefold() == "candidatus":
        # "Candidatus" is a nomenclatural status prefix for uncultured organisms,
        # not a genus. The genus is the prefix plus the next word, so a two-word
        # name like "Candidatus Soleaferrea" is a bare genus with no epithet.
        # Treating it as genus="Candidatus" invents a genus that does not exist
        # and collapses unrelated organisms under it.
        if len(parts) == 2:
            return normalize_genus(" ".join(parts[:2])), ""
        if len(parts) >= 3:
            genus = normalize_genus(" ".join(parts[:2]))
            return genus, " ".join(parts[2:]).lower()
    genus = normalize_genus(parts[0])
    return genus, " ".join(parts[1:]).lower()


# A source field holding "Genus epithet" names its own genus, so the `genus`
# column is not needed to read it. These guards keep that reading to fields that
# really are binomials: a rank name above genus is not a genus (NCBI calls an
# unnamed member of a family "Lachnospiraceae bacterium", which must not become
# genus Lachnospiraceae), and neither is the first word of a virus name.
RANK_NAME_SUFFIXES = (
    "aceae", "ales", "ineae", "idae", "aria", "eae",
    "mycetes", "mycota", "bacteria", "bacteriota",
)
NON_EPITHETS = {
    "bacterium", "archaeon", "symbiont", "endosymbiont", "oral", "group",
    "complex", "sp", "spp", "subsp", "str", "var", "cf", "aff",
}
# "Parabacteroides sp. 20_3" names its genus as surely as a binomial does; the
# unnamed species is still that genus's. Only these stand in for an epithet --
# "subsp."/"var." qualify an epithet rather than replace one, so a field that
# opens with them is not naming a species at all.
EPITHET_PLACEHOLDERS = {"sp.", "sp", "spp.", "spp"}


def names_a_genus(token: Any) -> bool:
    """True when a token can be read as a genus name on its own."""
    text = clean_str(token).strip("[]")
    if not re.fullmatch(r"[A-Za-z]{4,}", text):
        return False
    return not text.casefold().endswith(RANK_NAME_SUFFIXES)


def is_species_epithet(token: Any) -> bool:
    """True when a token can be read as a species epithet."""
    text = clean_str(token).strip("[]").casefold()
    if not re.fullmatch(r"[a-z]{3,}", text):
        return False
    return text not in NON_EPITHETS and "virus" not in text and "phage" not in text


def resolve_source_binomial(genus_field: Any, species_field: Any) -> tuple[str, str]:
    """Read (genus, epithet) from a source row's `genus` and `species` fields.

    MiMeDB stores a full binomial in `species` while `genus` can hold a
    different genus for the same organism -- the two fields disagree because one
    of them was updated for a reclassification and the other was not. Keeping
    both, which is what `extract_species_epithet` alone does (it only strips the
    leading token when it equals the `genus` column), concatenates two genera
    into a name no organism bears: "Lactobacillus limosilactobacillus reuteri",
    "Clavispora candida auris".

    Which field is current cannot be told from the fields themselves: MiMeDB has
    the new name in `species` for Limosilactobacillus reuteri and the old one
    there for Bacteroides vulgatus (where `genus` says Phocaeicola), and
    MMDBm0000292 is the reverse again. So this does not try to pick the current
    name -- it picks the field that names one organism. Taking the binomial as
    written always yields a real name, at worst a synonym, which
    `merge_duplicate_taxa.py` folds onto the current name by tax ID. Keeping the
    `genus` column and dropping the foreign token instead would invent names:
    "Clavispora auris", "Alfipia birgiae", "Brochotrix thermosphacta".

    An empty `genus` column is left empty rather than filled from the species
    field: `_collapse_mimedb` skips those rows, and inventing a genus for them
    would ingest ~40 organisms the loader has never loaded.
    """
    genus = normalize_genus(genus_field)
    field = clean_str(species_field)
    if not genus:
        return genus, extract_species_epithet(genus, field)

    parts = field.replace("[", "").replace("]", "").split()
    def follows_a_genus(token: str) -> bool:
        return is_species_epithet(token) or token.casefold() in EPITHET_PLACEHOLDERS

    if len(parts) >= 3 and parts[0].casefold() == "candidatus":
        # "Candidatus" is a status prefix, so the genus is the next word too.
        takes_field = (
            " ".join(parts[:2]).casefold() != genus.casefold()
            and names_a_genus(parts[1])
            and follows_a_genus(parts[2])
        )
    elif len(parts) >= 2:
        takes_field = (
            parts[0].casefold() != genus.casefold()
            and names_a_genus(parts[0])
            and follows_a_genus(parts[1])
        )
    else:
        takes_field = False

    if takes_field:
        field_genus, field_species = parse_scientific_name(field)
        if field_genus and field_species:
            return field_genus, field_species
    return genus, extract_species_epithet(genus, field)


def taxon_key(genus: Any, species: Any) -> tuple[str, str]:
    return normalize_genus(genus).casefold(), extract_species_epithet(genus, species).casefold()


def derive_energy_mode(
    oxygen_requirement: Any = "",
    metabolism: Any = "",
    energy_source: Any = "",
    activity: Any = "",
) -> str:
    oxygen = clean_str(oxygen_requirement).casefold()
    evidence = " ".join(
        clean_str(x).casefold() for x in (metabolism, energy_source, activity)
    )

    fermenter_terms = (
        "ferment",
        "saccharolytic",
        "cellulolytic",
        "ligninolytic",
    )
    respirer_terms = (
        "respirat",
        "non-ferment",
        "nonferment",
        "lithotroph",
        "chemolith",
        "phototroph",
        "autotroph",
    )
    fermenter = any(term in evidence for term in fermenter_terms)
    respirer = any(term in evidence for term in respirer_terms)

    if "facultative" in oxygen or (fermenter and respirer):
        return "mixed"
    if fermenter:
        return "fermenter"
    if respirer:
        return "respirator"
    if "anaerob" in oxygen and "aerotolerant" not in oxygen:
        return "fermenter"
    if any(term in oxygen for term in ("aerobe", "aerobic", "microaeroph")):
        return "respirator"
    return "N/A"


def derive_primary_food_source(
    metabolism: Any = "", activity: Any = "", energy_source: Any = ""
) -> str:
    evidence = " ".join(
        clean_str(x).casefold() for x in (metabolism, activity, energy_source)
    )
    rules = (
        (("cellulolytic", "cellulose", "ligninolytic"), "complex_carbs_fiber"),
        (("saccharolytic", "carbohydrate", "glycolytic"), "carbohydrates"),
        (("proteolytic", "protein", "amino acid", "peptidolytic"), "amino_acids_proteins"),
        (("methylotroph", "methanotroph", "c1 compound"), "c1_compounds"),
        (("lithotroph", "chemolith", "inorganic"), "inorganic"),
    )
    for terms, category in rules:
        if any(term in evidence for term in terms):
            return category
    if any(term in evidence for term in ("heterotroph", "organotroph", "osmotroph")):
        return "organic_unspecified"
    return "N/A"


def collapse_energy_modes(values: Iterable[Any]) -> str:
    modes = {clean_str(v).casefold() for v in values}
    modes -= NA_VALUES
    if not modes:
        return "N/A"
    if "mixed" in modes or len(modes) > 1:
        return "mixed"
    return next(iter(modes))


FOOD_PRIORITY = {
    "complex_carbs_fiber": 0,
    "carbohydrates": 1,
    "amino_acids_proteins": 2,
    "c1_compounds": 3,
    "inorganic": 4,
    "organic_unspecified": 5,
}


def collapse_food_sources(values: Iterable[Any]) -> str:
    usable = [clean_str(v).casefold() for v in values]
    usable = [v for v in usable if v not in NA_VALUES]
    if not usable:
        return "N/A"
    counts = Counter(usable)
    return min(counts, key=lambda item: (-counts[item], FOOD_PRIORITY.get(item, 99), item))


# The 2021 ICNP revision renamed most bacterial phyla (Firmicutes -> Bacillota and
# so on). NCBI now returns the new names, MiMeDB carries a mix of both, and the
# microbiome literature still overwhelmingly uses the classic ones (the
# Firmicutes/Bacteroidetes ratio being the obvious example). Everything is folded
# to the classic name so a phylum groups as one value instead of two.
PHYLUM_SYNONYMS = {
    "acidobacteriota": "Acidobacteria",
    "actinomycetota": "Actinobacteria",
    "aquificota": "Aquificae",
    "armatimonadota": "Armatimonadetes",
    "bacillota": "Firmicutes",
    "bacteroidota": "Bacteroidetes",
    "caldisericota": "Caldiserica",
    # Same NCBI node (28889) under two names, so a rename like the others. The
    # thermophilic archaea filed here -- Aeropyrum, Metallosphaera, Ignicoccus,
    # Pyrobaculum -- are not gut organisms, but they come in through MiMeDB's
    # reference set and a split phylum is a split phylum.
    "thermoproteota": "Crenarchaeota",
    "chlamydiota": "Chlamydiae",
    "chlorobiota": "Chlorobi",
    "chloroflexota": "Chloroflexi",
    "chrysiogenota": "Chrysiogenetes",
    "cyanobacteriota": "Cyanobacteria",
    "deferribacterota": "Deferribacteres",
    "deinococcota": "Deinococcus-Thermus",
    "dictyoglomota": "Dictyoglomi",
    "elusimicrobiota": "Elusimicrobia",
    "fibrobacterota": "Fibrobacteres",
    "fusobacteriota": "Fusobacteria",
    "gemmatimonadota": "Gemmatimonadetes",
    "ignavibacteriota": "Ignavibacteriae",
    "lentisphaerota": "Lentisphaerae",
    # Euryarchaeota and Methanobacteriota are NCBI taxid 28890 under two names,
    # so this is a rename and folds safely. Its neighbours are not: GTDB split
    # the old Euryarchaeota, and Thermoplasmatota (2283796) and Halobacteriota
    # (3100166) are separate taxids, i.e. real reclassifications. Those belong in
    # scripts/normalize_reclassified_phyla.py's APPROVED set, where a human signs
    # off on moving taxa between phyla, not in a synonym table.
    "methanobacteriota": "Euryarchaeota",
    "mycoplasmatota": "Tenericutes",
    "nitrospinota": "Nitrospinae",
    "nitrospirota": "Nitrospirae",
    "planctomycetota": "Planctomycetes",
    "pseudomonadota": "Proteobacteria",
    # Not a name anyone publishes: a one-character slip for Pseudomonadota that
    # reached the taxa table and would otherwise sit as an 87th phylum of one.
    "pseudomonadata": "Proteobacteria",
    "spirochaetota": "Spirochaetes",
    # A Candidatus prefix marks an uncultured lineage, not a different phylum.
    "candidatus thermoplasmatota": "Thermoplasmatota",
    "synergistota": "Synergistetes",
    # NOT "Thermodesulfobacteria". The ICNP did rename that classic phylum to
    # Thermodesulfobacteriota, but NCBI subsequently folded Desulfobacterota --
    # the former Deltaproteobacteria sulfate reducers -- into the same name, so
    # NCBI's Thermodesulfobacteriota is now a broad phylum, not the small
    # thermophile one. Folding it onto the classic name filed Desulfovibrio,
    # Bilophila and Desulfobulbus (gut sulfate reducers studied in IBD and CRC)
    # as hot-spring thermophiles. Left unmapped so the source's own name stands.
    "thermotogota": "Thermotogae",
    "verrucomicrobiota": "Verrucomicrobia",
}


# NCBI appends the nomenclatural code to names published outside the ICNP, so
# the dump carries "Methylomirabilota (SeqCode)" and nine siblings across the
# class, order and family ranks. The code is metadata about how a name was
# published, not part of the name, and leaving it attached creates a rank value
# of one that would split from the bare name the moment anything else arrives
# under it. Only these two are stripped, and only in trailing position: a
# parenthetical is not always an annotation, and "Neobittarella (ex bilen et al.
# 2018)" is a scientific name this database legitimately held.
NOMENCLATURAL_CODE = re.compile(r"\s*\((?:SeqCode|ICNP)\)\s*$", re.IGNORECASE)


def strip_nomenclatural_code(value: Any) -> str:
    """Drop a trailing (SeqCode) or (ICNP) annotation from a lineage rank."""
    return NOMENCLATURAL_CODE.sub("", clean_str(value)).strip()


def normalize_phylum(value: Any) -> str:
    """Fold post-2021 phylum names onto the classic name already used here."""
    text = strip_nomenclatural_code(value)
    if not text:
        return ""
    return PHYLUM_SYNONYMS.get(text.casefold(), text)


LATIN_INFLECTIONS = ("us", "um", "a", "is", "e", "ii", "i", "ae", "os", "on", "s", "ys")


def latin_stem(word: str) -> str:
    """Strip a trailing Latin case/gender ending, keeping at least four letters."""
    text = clean_str(word).casefold()
    for suffix in sorted(LATIN_INFLECTIONS, key=len, reverse=True):
        if text.endswith(suffix) and len(text) - len(suffix) >= 4:
            return text[: -len(suffix)]
    return text


def same_epithet(first: Any, second: Any) -> bool:
    """True when two species epithets differ only by Latin gender agreement.

    A species moved to a new genus takes that genus's grammatical gender, so
    `Erysipelatoclostridium ramosum` becomes `Thomasclavelia ramosa` and
    `Butyrivibrio crossotus` becomes `Eshraghiella crossota` -- same organism.
    A genuinely different epithet (`prausnitzii` vs `duncaniae`) is a different
    organism, and conflating the two relabels real measurements, so this is
    deliberately strict: it only tolerates a shared stem of four or more
    letters, never a prefix overlap.
    """
    one, two = clean_str(first).casefold(), clean_str(second).casefold()
    if not one or not two:
        return False
    if one == two:
        return True
    stem_one, stem_two = latin_stem(one), latin_stem(two)
    return len(stem_one) >= 4 and stem_one == stem_two


PLATFORM_SYNONYMS = {
    "16s": "16S",
    "16s rrna": "16S",
    "amplicon": "16S",
    "mngs": "mNGS",
    "metagenomics": "mNGS",
    "shotgun metagenomics": "mNGS",
    "wgs": "mNGS",
    "meta-analysis": "meta-analysis",
}


def normalize_platform(value: Any) -> str:
    """Fold a study's sequencing platform onto one of 16S / mNGS.

    The distinction drives the whole species-vs-genus feature space, so it has
    to be one value per platform. GMrepo has reported this field as both
    "16S"/"mNGS" and "AMPLICON"/"Metagenomics" at different times, and a
    re-sync that wrote the second spelling silently split 229 16S studies into
    191 AMPLICON plus 48 16S -- enough to drop a third of them out of a
    `--data-type 16S` filter without any error. Normalizing on the way in means
    the vocabulary cannot drift again whatever the upstream calls it.
    """
    text = clean_str(value)
    if not text:
        return ""
    return PLATFORM_SYNONYMS.get(text.casefold(), text)


def normalize_yes_no(value: Any) -> str:
    """Normalize a MiMeDB yes/no trait, discarding uncertain calls.

    MiMeDB marks unsure observations with a trailing '?' ("Yes?", "No?"). Those
    are dropped rather than stripped to a bare Yes/No, since promoting a hedged
    call to a definite one invents certainty the source never claimed.
    """
    text = clean_str(value)
    if not text or "?" in text:
        return ""
    folded = text.casefold()
    if folded in {"yes", "y", "true", "1"}:
        return "Yes"
    if folded in {"no", "n", "false", "0"}:
        return "No"
    return ""


def normalize_pathogen_flag(value: Any) -> int | None:
    """Return 1 only where MiMeDB positively documents a human pathogen.

    The source records the flag exclusively as a positive assertion, so anything
    else is unknown and stays NULL. Never return 0: that would fabricate negative
    evidence for ~97% of taxa that simply were not annotated.
    """
    return 1 if clean_str(value) == "1" else None


def nullable_float(value: Any) -> float | None:
    text = clean_str(value)
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


# Methods the literature marker files name in their own `effect_type` column.
# That column answers "how was this obtained", which is not what
# taxon_disease_associations.effect_type records -- that says WHAT the number in
# effect_size is. Copying it verbatim put 125 rows under effect_type 'LEfSe', a
# method standing beside LDA and log2_fold_change, and 124 of them had no effect
# size at all. The method is not lost by fixing this: phenotype_comparisons.method
# already carries it, and scripts/cleanup_gutdb.py reads it from the same column
# for exactly that purpose.
EFFECT_TYPE_METHODS = frozenset({"lefse", "lefse analysis", "linear discriminant analysis"})


def resolve_effect_type(declared: Any, effect_size: Any = None) -> str:
    """What the stored number is, rather than how it was produced.

    A declared method becomes 'LDA' when a score came with it and 'reported'
    when none did -- which is what the rows loaded from these same files already
    use for a direction with no magnitude.
    """
    text = clean_str(declared).casefold()
    if not text:
        return "LDA"
    if text in EFFECT_TYPE_METHODS:
        return "LDA" if effect_size is not None else "reported"
    return clean_str(declared)


def quantize_effect_size(value: Any) -> float | None:
    """Quantize an effect size to two decimals: round positives, floor negatives.

    Negatives are floored rather than rounded so a score's magnitude never
    shrinks: an LDA score's absolute value is the strength of the effect, so
    -4.301 becomes -4.31, whereas rounding to -4.30 would understate it.
    Positives use round-half-up rather than Python's default round-half-even,
    which would send 2.345 to 2.34.
    """
    number = nullable_float(value)
    if number is None:
        return None
    quantized = Decimal(str(number)).quantize(
        Decimal("0.01"),
        rounding=ROUND_FLOOR if number < 0 else ROUND_HALF_UP,
    )
    return float(quantized)


def nullable_int(value: Any) -> int | None:
    number = nullable_float(value)
    return int(number) if number is not None else None


def association_direction(
    lda_score: float | None,
    disease_name: str,
    positive_enriched_in: str = "",
    negative_enriched_in: str = "",
    explicit_direction: str = "",
) -> str:
    explicit = clean_str(explicit_direction).casefold()
    if explicit in {"enriched", "depleted", "marker", "no_difference"}:
        return explicit
    if lda_score is None or lda_score == 0:
        return "marker"
    enriched_group = positive_enriched_in if lda_score > 0 else negative_enriched_in
    if not clean_str(enriched_group):
        return "marker"
    return (
        "enriched"
        if clean_str(enriched_group).casefold() == clean_str(disease_name).casefold()
        else "depleted"
    )


# --- controlled vocabularies for the taxa trait columns ----------------------
#
# These lived in scripts/cleanup_gutdb.py, which runs only against the SQLite
# mirror, so for three columns the two stores held different vocabularies for the
# same rows. MySQL carried 30 distinct oxygen_requirement values where the mirror
# carried 10 -- 'anaerobe' (576), 'Anaerobic' (50) and 'Obligate anaerobe' (35)
# as three groups rather than one -- and `GROUP BY oxygen_requirement` therefore
# answered differently depending on which store was asked. Oxygen tolerance is
# one of the few genuinely mechanistic traits for a gut community, so a
# fragmented categorical is a real loss, not untidiness.
#
# They are here rather than in a script because this is the third repair that
# existed for one store only: the (SeqCode) annotations and body_site were the
# first two. One implementation, imported by both sides.
#
# Case and synonym variants collapse onto one lowercase canonical form. The
# right-hand sides are the vocabulary; a value the vocabulary does not list is
# kept, lowercased, so an unanticipated category is never silently dropped.
OXYGEN = {
    "aerobe": "aerobe", "aerobic": "aerobe", "obligate aerobe": "aerobe",
    "obligate aerobic": "aerobe",
    "anaerobe": "anaerobe", "anaerobic": "anaerobe",
    "obligate anaerobe": "anaerobe",
    "facultative": "facultative anaerobe",
    "facultative anaerobe": "facultative anaerobe",
    "facultative aerobe": "facultative anaerobe",
    "facultatively anaerobe": "facultative anaerobe",
    "aerobe/facultative anaerobe": "facultative anaerobe",
    "aerobe; facultative anaerobe": "facultative anaerobe",
    "facultative anaerobe, air+co2": "facultative anaerobe",
    "microaerophile": "microaerophile", "microaerophilic": "microaerophile",
    # Organisms reported as tolerating both: 'variable' rather than a pick, since
    # collapsing them onto either end would assert something the source does not.
    "aerobe; anaerobe": "variable",
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
    "mesophilic": "mesophilic", "mesophil": "mesophilic",
    "thermophilic": "thermophilic", "hyperthermophilic": "hyperthermophilic",
    "psychrophilic": "psychrophilic", "psychrotolerant": "psychrotolerant",
}
YES_NO = {"yes": "yes", "no": "no", "yes?": "yes", "no?": "no",
          "1": "yes", "0": "no", "true": "yes", "false": "no", "not known": ""}
PH = {"alkaliphile": "alkaliphile", "acidophile": "acidophile",
      "neutrophile": "neutrophile"}
# `superkingdom` collects three different vocabularies: MiMeDB's domains, the
# curated table's `kingdom` (Eubacteria, Fungi), and NCBI's post-2024 kingdom
# names (Bacillati, Pseudomonadati). All fold onto the three domains + Viruses.
SUPERKINGDOM = {
    "bacteria": "Bacteria", "eubacteria": "Bacteria", "bacillati": "Bacteria",
    "pseudomonadati": "Bacteria", "fusobacteriati": "Bacteria",
    "thermotogati": "Bacteria",
    "archaea": "Archaea", "methanobacteriati": "Archaea",
    "thermoproteati": "Archaea",
    "eukaryota": "Eukaryota", "fungi": "Eukaryota", "metazoa": "Eukaryota",
    "viridiplantae": "Eukaryota", "viruses": "Viruses",
}
# Terms that describe how cells sit next to each other, which is a different axis
# from whether they form spores. A `sporulation` value made only of these is
# cell-arrangement text in the wrong column.
ARRANGEMENT_TERMS = frozenset({"singles", "pairs", "chains", "clusters",
                               "tetrads", "filaments", "v-shaped forms"})

TRAIT_VOCABULARIES = {
    "oxygen_requirement": OXYGEN,
    "gram_stain": GRAM,
    "ph_preference": PH,
    "superkingdom": SUPERKINGDOM,
    "temperature_range": TEMPERATURE,
    "mobility": YES_NO,
    "flagella_presence": YES_NO,
}
# Columns that get whitespace and case tidied but have no closed vocabulary --
# free text whose categories nobody has enumerated.
FREE_TEXT_TRAITS = ("human_pathogen", "shape", "metabolism", "habitat",
                    "energy_source", "biotic_relationship",
                    "optimal_temperature", "cell_arrangement")


def canon_trait(value, mapping: dict[str, str] | None = None) -> str:
    """Clean a trait value; map it through a vocabulary when one is given.

    Under a vocabulary, an unlisted value is kept but lowercased, so
    'Aerotolerant' and 'aerotolerant' stop being two values while a category
    nobody anticipated survives. SUPERKINGDOM is the exception: its canonical
    forms are capitalised because they are taxon names.
    """
    text = clean_str(value).replace("\xa0", " ")
    text = re.sub(r"\s+", " ", text).strip()
    if not text:
        return ""
    if mapping is None:
        return text
    return mapping.get(text.casefold(),
                       text if mapping is SUPERKINGDOM else text.casefold())


def split_sporulation(value) -> tuple[str, str]:
    """Split one mixed `Spore`-style value into (sporulation, cell_arrangement).

    The column holds three kinds of value at once: a yes/no flag, prose about
    spore type ('Sporulating (ascospores); reproduction also by budding'), and
    cell-arrangement text that belongs to another column entirely ('Chains -
    Pairs - Singles'). Anything the three rules below do not recognise returns
    ('', ''), which leaves the cell alone rather than guessing at it.
    """
    text = canon_trait(value)
    if not text:
        return "", ""
    folded = text.casefold()
    if folded in YES_NO:
        return YES_NO[folded], ""
    if "sporulat" in folded or "spotulating" in folded:
        return ("no" if folded.startswith("non") else "yes"), ""
    if "budding" in folded:
        return ("no" if folded.startswith(("no;", "no ")) else "yes"), ""
    # Two delimiter styles are in use and a single split handles only one of
    # them. 'Chains - Pairs - Singles' and 'Singles - V-shaped forms' separate on
    # ' - ', while 'Pairs-Singles' separates on a bare hyphen -- and splitting
    # 'Singles - V-shaped forms' on the bare hyphen yields 'v' and 'shaped
    # forms', which matches nothing. Try the spaced form first, then the bare one.
    for parts in (folded.replace("/", " - ").split(" - "),
                  folded.replace("/", "-").split("-")):
        tokens = {p.strip() for p in parts if p.strip()}
        if tokens and tokens <= ARRANGEMENT_TERMS:
            return "", text
    return "", ""


# Values that name a study design or summarise several platforms rather than
# naming one. A SAMPLE must never inherit any of these: 'mixed' is not a
# resolution and not a batch.
NOT_A_PLATFORM = frozenset({"mixed", "meta-analysis", ""})


def derive_study_platform(platforms) -> str:
    """Summarise a study's per-sample platforms into its data_type.

    One platform across the samples means that platform; more than one means
    'mixed'. Keeping this as a function rather than inline SQL is the point --
    the SQLite mirror held FOUR vocabularies for two platforms ('AMPLICON'
    beside '16S rRNA', 'Metagenomics' beside 'shotgun metagenomics') where MySQL
    held the normalized pair, so `WHERE data_type = '16S'` matched 236 studies
    in one store and none in the other. Deriving both stores' value from the same
    rule is what stops that recurring.

    Returns '' for a study with nothing to derive from, which the callers read as
    "leave this study's data_type alone" -- an association-only study has no
    samples, and 'meta-analysis' is a true statement about those.
    """
    distinct = {normalize_platform(p) for p in platforms if p}
    distinct -= {""}
    if not distinct:
        return ""
    if len(distinct) == 1:
        return distinct.pop()
    return "mixed"


# --- sample_attributes.value_numeric ----------------------------------------
#
# The column's contract is "a number you may do arithmetic on". Two classes of
# value break it while looking fine:
#
#   The number was FABRICATED from a string that is not one. 434 rows under
#   `unique_id` hold '10007539_0' with value_numeric = 100075390 -- the
#   underscore was dropped and the digits run together, which can collide two
#   distinct ids ('10007539_0' and '1000753_90' both give 100075390).
#
#   The number is a CODE, not a quantity, and the key name says so:
#   'Antibiotic_in_the_past_month._1_no_2_yes' stores 1 and 2. Averaging those
#   is averaging labels, and worse, 1 means NO there while 1 means YES under the
#   neighbouring `Antibiotics` key -- so pooling the two inverts 343 rows. Same
#   family of error as the LDA sign ambiguity this schema already guards.
#
# A key whose values MIX numbers and words is a third case and is deliberately
# NOT cleared: the number is correct for the rows that have one. What is missing
# there is the denominator -- `Stage` holds 0 beside I, IIA and IIIB, so
# AVG(value_numeric) is the mean of 27 zeros out of 150 rows, which is 0 -- and
# the fix for a missing denominator is to expose it, exactly as with absence in
# the abundance matrix. v_sample_attribute_keys reports it per key.
PLAIN_NUMBER = re.compile(r"^\s*-?\d+(\.\d+)?\s*$")
# A key name that spells out its own code mapping: '..._1_no_2_yes', 'yes_1_no_0'.
CODED_IN_NAME = re.compile(r"\d\s*_\s*(?:no|yes)|(?:no|yes)\s*_\s*\d", re.IGNORECASE)


def numeric_value_is_trustworthy(attribute: str, value: str) -> bool:
    """Whether value_numeric may hold a number for this (attribute, value).

    False means the cell should be NULL: either the string is not the number
    stored against it, or the number is a category code whose meaning lives in
    the key name rather than on the number line.
    """
    if not PLAIN_NUMBER.match(clean_str(value)):
        return False
    if CODED_IN_NAME.search(clean_str(attribute)):
        return False
    return True


# --- metabolite level scales -------------------------------------------------
#
# sample_metabolite_levels mixes scales that differ in KIND, not magnitude, and
# nothing recorded which. Two of the 13 cohorts arrive already log- or
# z-transformed and hold 70,067 negative values; the other 11 are non-negative
# intensities whose maxima span nine orders of magnitude, from 461 in
# KANG_AUTISM_2017 to 3.07e10 in KOSTIC_INFANTS_DIABETES_2015.
#
# THE DANGEROUS PART IS THE ZERO. On a non-negative scale a zero means the
# compound was not detected -- 78,400 of YACHIDA_CRC_2019's 130,242 levels, 60%
# of them. On a centred scale a zero is mid-distribution: SINHA_CRC_2016 has
# 23,553 negatives against 22,540 positives with a mean of -0.04, so its 281
# zeros mean "at this cohort's average", the opposite of absence. The same stored
# number means opposite things depending on the study, and the standing advice to
# "normalize within study" does not cover it -- a uniform log1p over this table
# yields NaN on 70,067 rows and takes the log of a log on 136,134 more.
METABOLITE_SCALE_CENTRED = "centred"
METABOLITE_SCALE_NON_NEGATIVE = "non_negative"


def derive_metabolite_scale(min_level, max_level) -> str:
    """Classify a study's level scale from its own range.

    A negative level cannot be a concentration, so it is positive evidence that
    the cohort was transformed upstream. The absence of negatives is NOT evidence
    of a particular unit, only that the values are non-negative -- which is why
    the other label says exactly that and claims no unit. The source does not
    record units anywhere; studies.metabolomics_method carries the protocol DOI.
    """
    if min_level is None:
        return ""
    if float(min_level) < 0:
        return METABOLITE_SCALE_CENTRED
    return METABOLITE_SCALE_NON_NEGATIVE


def zero_means_absence(scale: str) -> bool:
    """Whether a stored level of 0 may be read as "not detected" on this scale."""
    return scale == METABOLITE_SCALE_NON_NEGATIVE


# --- NCBI rank vs the abundance layer ---------------------------------------
#
# taxonomic_rank is the LAYER a row was loaded into -- GMrepo publishes a genus
# array and a species array, and the per-(sample, rank) sum to 1 is defined over
# that split. ncbi_rank is what the organism actually is. 408 rows differ.
#
# Ranks at or below genus: a genus-layer row carrying one of these is not coarser
# than its layer. 'no rank' and 'clade' assert nothing either way, so they are
# included here rather than counted as a disagreement.
RANKS_AT_OR_BELOW_GENUS = frozenset({
    "genus", "subgenus", "species", "subspecies", "strain", "species group",
    "species subgroup", "isolate", "serotype", "no rank", "clade",
})


def ncbi_rank_differs(ncbi_rank, layer_rank: str) -> bool:
    """Whether NCBI puts this taxon at a rank other than the layer it sits in.

    Covers both directions, because both are the same statement -- "this feature
    is not what its layer says" -- and only one of them was exposed:

      COARSER. 123 genus-layer features are families, orders, classes or phyla,
      carrying 16,098 abundance rows. coarser_than_genus() has flagged these
      since they were found.

      FINER. 4 genus-layer features are species, carrying 1,233 rows: the three
      bracket organisms NCBI leaves without a genus, plus Phenoliferia
      psychrophila. Nothing flagged them, and three of them exist because
      unwind_bracket_genus_repoint.py created them -- so the gap was opened by
      the fix for a different one.

      And 181 SPECIES-layer features are strains, carrying 189,390 rows, 14.8%
      of that layer, with v_abundance_species exposing neither the rank nor a
      flag.

    'no rank' and 'clade' assert nothing, so they are not a difference. An empty
    rank is not either -- 39 taxa carry a tax ID the local dump does not list.
    """
    rank = clean_str(ncbi_rank).casefold()
    if not rank or rank in ("no rank", "clade"):
        return False
    return rank != clean_str(layer_rank).casefold()


def coarser_than_genus(ncbi_rank) -> bool:
    """Whether a genus-layer feature is really a family, order, class or phylum.

    123 observed taxa are: 90 families (Muribaculaceae, Tannerellaceae), 19
    orders (Candidatus Gastranaerophilales, Bivalvulida), 7 classes, 5 phyla.
    They carry 16,098 abundance rows and about 1% of the genus layer's mass.

    They are NOT errors to move out of the layer. GMrepo labels them genus, they
    are the unclassified residual within their family rather than a duplicate of
    it -- mean shares of 0.2-1.6% beside real genera carrying far more -- and the
    parts are disjoint, so the sum to 1 holds. The only thing missing was a way
    to tell them apart from a real genus, which is what this is for.
    """
    rank = clean_str(ncbi_rank).casefold()
    if not rank:
        return False
    return rank not in RANKS_AT_OR_BELOW_GENUS


# --- class, order and family synonyms ---------------------------------------
#
# PHYLUM_SYNONYMS solved this for one rank. The other three never got it:
# scripts/backfill_lineage.py passes ONLY `phylum` through normalize_phylum and
# writes class, order and family verbatim from the dump, so a row filled from
# NCBI carries the current scientific name while one loaded earlier keeps the
# old one, and the two sit side by side. Measured by the tax ID identity test --
# two stored names resolving to the SAME NCBI node are synonyms of one taxon,
# never a reclassification:
#
#     class_name    1 node  split over   409 rows
#     order_name   10 nodes split over 1,500 rows
#     family       13 nodes split over   557 rows
#
# So `GROUP BY order_name` reported Clostridiales as two orders -- 413 rows under
# Eubacteriales and 127 under Clostridiales -- which is the Firmicutes/Bacillota
# defect again one rank down.
#
# THE TARGET IS THE MAJORITY STORED NAME, NOT NCBI'S SCIENTIFIC NAME, and that
# is deliberate. Folding onto the scientific name would have:
#
#   re-introduced the nomenclatural code this schema strips --
#   Saccharimonadales (SeqCode), Saccharimonadaceae (SeqCode) and
#   Obscuribacteraceae (SeqCode) are NCBI's scientific names for three of these
#   nodes, and writing them back would fail the pre-push phylum check; and
#
#   moved 107 Streptomycetales rows onto Kitasatosporales, which NCBI prefers and
#   almost nobody uses -- the same reason PHYLUM_SYNONYMS keeps Firmicutes over
#   Bacillota.
#
# Every target below is a name already present in this database. The map is
# explicit rather than computed at run time so the fold cannot drift as row
# counts change; the counts in the comments are from 2026-10-14.
LINEAGE_SYNONYMS = {
    # order: the 2021-era renamings, each confirmed by a shared tax ID
    "clostridiales": "Eubacteriales",                 # 127 -> 413
    "bacillales": "Caryophanales",                    # 122 -> 167
    "rhizobiales": "Hyphomicrobiales",                #  36 -> 113
    "corynebacteriales": "Mycobacteriales",           #  26 ->  79
    "enterobacteriales": "Enterobacterales",          #  15 -> 230
    "xanthomonadales": "Lysobacterales",              #   4 ->  37
    "kitasatosporales": "Streptomycetales",           #   3 -> 107, against NCBI
    "erysiphales": "Helotiales",                      #   1 ->   6
    "haloferacales": "Halobacteriales",               #   1 ->   8
    "saccharimonadales": "Candidatus Saccharimonadales",   # 1 -> 4
    # family
    "ruminococcaceae": "Oscillospiraceae",            #  32 -> 161
    "leuconostocaceae": "Lactobacillaceae",           #   8 -> 189
    "planococcaceae": "Caryophanaceae",               #   8 ->  13
    "bradyrhizobiaceae": "Nitrobacteraceae",          #   6 ->   7
    "rhodobacteraceae": "Paracoccaceae",              #   6 ->  32
    "xanthomonadaceae": "Lysobacteraceae",            #   4 ->  26
    "hungateiclostridiaceae": "Oscillospiraceae",     #   3 -> 161
    "nocardiopsidaceae": "Nocardiopsaceae",           #   2 ->   3, against NCBI
    "clostridiales family xiii. incertae sedis": "Anaerovoracaceae",  # 1 -> 36
    "methanosaetaceae": "Methanotrichaceae",          #   1 ->   2
    "sinobacteraceae": "Nevskiaceae",                 #   1 ->   7
    "thiovulaceae": "Sulfurimonadaceae",              #   1 ->   2
    "saccharimonadaceae": "Candidatus Saccharimonadaceae",  # 1 -> 3
    "obscuribacteraceae": "Candidatus Obscuribacteraceae",  # 1 -> 1, tie
    # class
    "actinomycetia": "Actinomycetes",                 #   1 -> 408
    # The Candidatus status marker, which NCBI writes as part of the designation
    # for an uncultivated taxon. Neither store held both forms internally, but
    # they disagreed with each other across 5 values -- MySQL keeping the prefix,
    # the mirror dropping it -- so a lineage string read from one store did not
    # match the other. Folding onto the prefixed form matches NCBI and the
    # Saccharimonadales/Obscuribacteraceae entries above.
    "cloacimonadia": "Candidatus Cloacimonadia",
    "cloacimonadales": "Candidatus Cloacimonadales",
    "cloacimonadaceae": "Candidatus Cloacimonadaceae",
    "hepatoplasmataceae": "Candidatus Hepatoplasmataceae",
    "nucleicultricaceae": "Candidatus Nucleicultricaceae",
}


def normalize_lineage_name(value) -> str:
    """Fold a class, order or family name onto one name per NCBI node.

    Applied AFTER strip_nomenclatural_code, so a stored '(SeqCode)' value folds
    on its bare name. Returns the value unchanged when it is not a known synonym,
    so a name nobody has mapped is never discarded.
    """
    text = strip_nomenclatural_code(value)
    if not text:
        return ""
    return LINEAGE_SYNONYMS.get(text.casefold(), text)
