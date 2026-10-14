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
