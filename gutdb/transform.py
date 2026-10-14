from __future__ import annotations

from collections import Counter
from decimal import Decimal, ROUND_FLOOR, ROUND_HALF_UP
import math
import re
from typing import Any, Iterable, Mapping


NA_VALUES = {"", "n/a", "na", "none", "null", "nan", "unknown", "not available"}


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
