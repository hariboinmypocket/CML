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
    "mycoplasmatota": "Tenericutes",
    "nitrospinota": "Nitrospinae",
    "nitrospirota": "Nitrospirae",
    "planctomycetota": "Planctomycetes",
    "pseudomonadota": "Proteobacteria",
    "spirochaetota": "Spirochaetes",
    "synergistota": "Synergistetes",
    "thermodesulfobacteriota": "Thermodesulfobacteria",
    "thermotogota": "Thermotogae",
    "verrucomicrobiota": "Verrucomicrobia",
}


def normalize_phylum(value: Any) -> str:
    """Fold post-2021 phylum names onto the classic name already used here."""
    text = clean_str(value)
    if not text:
        return ""
    return PHYLUM_SYNONYMS.get(text.casefold(), text)


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
