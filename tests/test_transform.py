import unittest
import csv
from pathlib import Path

from gutdb.transform import (
    GRAM,
    NOT_A_PLATFORM,
    derive_metabolite_scale,
    derive_study_platform,
    numeric_value_is_trustworthy,
    zero_means_absence,
    OXYGEN,
    SUPERKINGDOM,
    canon_trait,
    coarser_than_genus,
    split_sporulation,
    quantize_effect_size,
    association_direction,
    collapse_energy_modes,
    collapse_food_sources,
    derive_energy_mode,
    derive_primary_food_source,
    extract_species_epithet,
    normalize_genus,
    normalize_phylum,
    resolve_effect_type,
    strip_nomenclatural_code,
    normalize_platform,
    parse_scientific_name,
    same_epithet,
)


class TransformTests(unittest.TestCase):
    def test_taxonomy_normalization(self):
        self.assertEqual(normalize_genus("  bAcTeRoIdEs  "), "Bacteroides")
        self.assertEqual(
            extract_species_epithet("Bacteroides", "Bacteroides ovatus"), "ovatus"
        )
        self.assertEqual(
            parse_scientific_name("[Clostridium] scindens"), ("Clostridium", "scindens")
        )
        self.assertEqual(
            parse_scientific_name("Oscillibacter sp. MSJ-31"),
            ("Oscillibacter", "sp. msj-31"),
        )
        self.assertEqual(
            parse_scientific_name("Candidatus Borkfalkia ceftriaxoniphila"),
            ("Candidatus Borkfalkia", "ceftriaxoniphila"),
        )

    def test_energy_mode(self):
        self.assertEqual(derive_energy_mode("Anaerobe"), "fermenter")
        self.assertEqual(derive_energy_mode("Aerobe"), "respirator")
        self.assertEqual(derive_energy_mode("Facultative anaerobe"), "mixed")
        self.assertEqual(derive_energy_mode("", "fermentation and respiration"), "mixed")
        self.assertEqual(collapse_energy_modes(["fermenter", "respirator"]), "mixed")

    def test_food_source(self):
        self.assertEqual(
            derive_primary_food_source("cellulolytic and saccharolytic"),
            "complex_carbs_fiber",
        )
        self.assertEqual(
            collapse_food_sources(["carbohydrates", "inorganic", "carbohydrates"]),
            "carbohydrates",
        )

    def test_direction_is_relative_to_disease(self):
        self.assertEqual(
            association_direction(-3.2, "IBS", "IBS", "Health"), "depleted"
        )
        self.assertEqual(
            association_direction(3.2, "IBS", "IBS", "Health"), "enriched"
        )

    def test_seed_marker_file_has_both_directions(self):
        path = Path(__file__).parents[1] / "data" / "gmrepo_PRJNA705217_ibs_markers.csv"
        with path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        directions = {
            association_direction(
                float(row["lda_score"]),
                row["disease_name"],
                row["positive_enriched_in"],
                row["negative_enriched_in"],
            )
            for row in rows
        }
        self.assertEqual(len(rows), 34)
        self.assertEqual(directions, {"enriched", "depleted"})


class QuantizeEffectSizeTests(unittest.TestCase):
    def test_positive_rounds_half_up(self):
        self.assertEqual(quantize_effect_size("2.3576382479"), 2.36)
        self.assertEqual(quantize_effect_size("2.345"), 2.35)  # not 2.34 (half-even)

    def test_negative_floors_away_from_zero(self):
        # a score's magnitude is the effect strength, so it must not shrink
        self.assertEqual(quantize_effect_size("-4.3012"), -4.31)
        self.assertEqual(quantize_effect_size("-2.031"), -2.04)

    def test_already_two_decimals_unchanged(self):
        self.assertEqual(quantize_effect_size("-4.31"), -4.31)
        self.assertEqual(quantize_effect_size(3.5), 3.5)

    def test_missing_and_unparseable_stay_none(self):
        for value in ("", None, "NA", "abc"):
            self.assertIsNone(quantize_effect_size(value))


class SameEpithetTests(unittest.TestCase):
    def test_gender_agreement_is_the_same_organism(self):
        # A species moved to a new genus takes that genus's gender.
        self.assertTrue(same_epithet("ramosum", "ramosa"))
        self.assertTrue(same_epithet("crossotus", "crossota"))
        self.assertTrue(same_epithet("firmus", "firmus"))

    def test_different_epithets_are_different_organisms(self):
        # The case that must never merge: tax ID 411483 is F. prausnitzii
        # A2-165, renamed F. duncaniae, while the species is tax ID 853.
        self.assertFalse(same_epithet("prausnitzii", "duncaniae"))
        self.assertFalse(same_epithet("hominis", "pacaense"))

    def test_short_or_prefix_overlap_is_not_enough(self):
        self.assertFalse(same_epithet("coli", "colitis"))
        self.assertFalse(same_epithet("avium", "avidum"))

    def test_empty_never_matches(self):
        for pair in (("", "ramosa"), ("ramosa", ""), ("", "")):
            self.assertFalse(same_epithet(*pair))


class NormalizePlatformTests(unittest.TestCase):
    def test_amplicon_spellings_fold_to_16S(self):
        for value in ("16S", "16S rRNA", "AMPLICON", "amplicon"):
            self.assertEqual(normalize_platform(value), "16S")

    def test_shotgun_spellings_fold_to_mNGS(self):
        for value in ("mNGS", "Metagenomics", "shotgun metagenomics", "WGS"):
            self.assertEqual(normalize_platform(value), "mNGS")

    def test_unknown_platform_survives(self):
        # An unanticipated platform is kept rather than coerced, so it shows up
        # as itself instead of silently joining 16S or mNGS.
        self.assertEqual(normalize_platform("Nanopore"), "Nanopore")
        self.assertEqual(normalize_platform("meta-analysis"), "meta-analysis")

    def test_blank(self):
        self.assertEqual(normalize_platform(""), "")
        self.assertEqual(normalize_platform(None), "")


class NormalizePhylumTests(unittest.TestCase):
    def test_classic_renames_fold_to_literature_name(self):
        self.assertEqual(normalize_phylum("Bacillota"), "Firmicutes")
        self.assertEqual(normalize_phylum("Bacteroidota"), "Bacteroidetes")
        self.assertEqual(normalize_phylum("Pseudomonadota"), "Proteobacteria")

    def test_thermodesulfobacteriota_is_not_folded(self):
        """NCBI's Thermodesulfobacteriota absorbed Desulfobacterota.

        Mapping it onto the classic small Thermodesulfobacteria phylum filed
        the gut sulfate reducers (Desulfovibrio, Bilophila) as thermophiles.
        """
        self.assertEqual(
            normalize_phylum("Thermodesulfobacteriota"), "Thermodesulfobacteriota"
        )
        self.assertEqual(
            normalize_phylum("Thermodesulfobacteria"), "Thermodesulfobacteria"
        )

    def test_euryarchaeota_rename_folds(self):
        """Same NCBI tax ID (28890) under two names, so it is a rename."""
        self.assertEqual(normalize_phylum("Methanobacteriota"), "Euryarchaeota")

    def test_nomenclatural_code_annotation_is_stripped(self):
        """NCBI appends the code to names published outside the ICNP.

        It describes how a name was published, not the name, and leaving it on
        makes a rank value of one that splits from the bare name later.
        """
        self.assertEqual(
            strip_nomenclatural_code("Methylomirabilota (SeqCode)"), "Methylomirabilota")
        self.assertEqual(strip_nomenclatural_code("Bathyarchaeia (ICNP)"), "Bathyarchaeia")
        self.assertEqual(normalize_phylum("Methylomirabilota (SeqCode)"), "Methylomirabilota")

    def test_other_parentheticals_survive(self):
        """A parenthetical is not always an annotation.

        "Neobittarella (ex bilen et al. 2018)" is a scientific name this database
        legitimately held, so the strip is deliberately limited to the two codes
        and to trailing position.
        """
        self.assertEqual(
            strip_nomenclatural_code("Neobittarella (ex bilen et al. 2018)"),
            "Neobittarella (ex bilen et al. 2018)")

    def test_crenarchaeota_rename_folds(self):
        """Thermoproteota and Crenarchaeota are both NCBI taxid 28889."""
        self.assertEqual(normalize_phylum("Thermoproteota"), "Crenarchaeota")

    def test_spelling_and_candidatus_variants_fold(self):
        """Neither is a phylum anyone publishes.

        'Pseudomonadata' is a one-character slip for Pseudomonadota that reached
        the taxa table, and a Candidatus prefix marks an uncultured lineage
        rather than a different phylum. Left alone, each sits as a phylum of one
        and splits a group.
        """
        self.assertEqual(normalize_phylum("Pseudomonadata"), "Proteobacteria")
        self.assertEqual(
            normalize_phylum("Candidatus thermoplasmatota"), "Thermoplasmatota")

    def test_euryarchaeota_splits_are_not_folded(self):
        """GTDB split the old Euryarchaeota, and NCBI gave the pieces their own
        nodes: Thermoplasmatota is 2283796 and Halobacteriota is 3100166, not
        28890. Folding either onto Euryarchaeota would move real organisms
        between phyla, which is a taxonomic decision rather than a spelling one
        and belongs in normalize_reclassified_phyla.py's APPROVED set.
        """
        self.assertEqual(normalize_phylum("Thermoplasmatota"), "Thermoplasmatota")
        self.assertEqual(normalize_phylum("Halobacteriota"), "Halobacteriota")
        self.assertEqual(normalize_phylum("Myxococcota"), "Myxococcota")

    def test_every_synonym_target_is_itself_stable(self):
        """Folding must reach a fixed point in one step.

        If a target were itself a key, the stored value would depend on how many
        times the normalizer happened to run -- and these names are written by
        several loaders.
        """
        from gutdb.transform import PHYLUM_SYNONYMS

        for source, target in PHYLUM_SYNONYMS.items():
            self.assertEqual(
                normalize_phylum(target), target,
                f"{source} -> {target}, but {target} folds again")

    def test_unknown_and_empty(self):
        self.assertEqual(normalize_phylum("Firmicutes"), "Firmicutes")
        self.assertEqual(normalize_phylum(""), "")
        self.assertEqual(normalize_phylum(None), "")


class ResolveEffectTypeTests(unittest.TestCase):
    """effect_type says WHAT the number is, not how it was produced."""

    def test_a_declared_method_becomes_the_statistic_it_produces(self):
        self.assertEqual(resolve_effect_type("LEfSe", 4.31), "LDA")

    def test_a_declared_method_with_no_score_is_reported(self):
        """124 of the 125 rows that carried effect_type 'LEfSe' had no score.

        A direction with no magnitude is exactly what 'reported' already means
        for the other rows loaded from these same literature files.
        """
        self.assertEqual(resolve_effect_type("LEfSe", None), "reported")

    def test_real_effect_types_pass_through(self):
        self.assertEqual(resolve_effect_type("log2_fold_change", 1.2), "log2_fold_change")
        self.assertEqual(resolve_effect_type("q_value", None), "q_value")
        self.assertEqual(resolve_effect_type("reported", None), "reported")

    def test_nothing_declared_defaults_to_lda(self):
        """GMrepo supplies LDA scores and declares no effect_type."""
        self.assertEqual(resolve_effect_type("", None), "LDA")
        self.assertEqual(resolve_effect_type(None, 2.0), "LDA")


if __name__ == "__main__":
    unittest.main()


class TraitVocabularyTests(unittest.TestCase):
    """The vocabularies moved here from scripts/cleanup_gutdb.py.

    They lived in a script that runs only against the SQLite mirror, so MySQL
    held 21 distinct oxygen_requirement values where the mirror held 10. These
    tests pin the folding itself, which is now the single implementation both
    stores use.
    """

    def test_oxygen_case_and_synonym_variants_are_one_value(self):
        for variant in ("anaerobe", "Anaerobic", "ANAEROBE", "Obligate anaerobe"):
            self.assertEqual(canon_trait(variant, OXYGEN), "anaerobe", variant)
        for variant in ("aerobe", "Aerobic", "obligate aerobe", "Obligate aerobic"):
            self.assertEqual(canon_trait(variant, OXYGEN), "aerobe", variant)

    def test_facultative_forms_collapse_onto_the_anaerobe_spelling(self):
        for variant in ("Facultative", "Facultative aerobe", "Facultatively anaerobe",
                        "Aerobe/facultative anaerobe", "Facultative anaerobe, air+co2"):
            self.assertEqual(canon_trait(variant, OXYGEN),
                             "facultative anaerobe", variant)

    def test_both_tolerances_reported_becomes_variable_not_a_pick(self):
        """Collapsing onto either end would assert what the source does not say."""
        self.assertEqual(canon_trait("Aerobe; anaerobe", OXYGEN), "variable")

    def test_an_unlisted_category_survives_lowercased(self):
        """Never silently dropped: 'Aerotolerant' is a real category nobody mapped."""
        self.assertEqual(canon_trait("Aerotolerant", OXYGEN), "aerotolerant")
        self.assertEqual(canon_trait("Nanaerobe", OXYGEN), "nanaerobe")

    def test_a_structurally_positive_wall_that_does_not_stain_is_variable(self):
        for variant in ("Structurally positive but stains negative",
                        "Structurally positive but may stain negative"):
            self.assertEqual(canon_trait(variant, GRAM), "variable", variant)

    def test_superkingdom_keeps_its_capitals_because_they_are_taxon_names(self):
        self.assertEqual(canon_trait("bacillati", SUPERKINGDOM), "Bacteria")
        self.assertEqual(canon_trait("Eubacteria", SUPERKINGDOM), "Bacteria")
        self.assertEqual(canon_trait("methanobacteriati", SUPERKINGDOM), "Archaea")

    def test_whitespace_and_nbsp_are_normalized(self):
        self.assertEqual(canon_trait("  Aerobic  "), "Aerobic")

    def test_empty_stays_empty_rather_than_becoming_a_value(self):
        for value in (None, "", "   ", " "):
            self.assertEqual(canon_trait(value, OXYGEN), "")


class SporulationSplitTests(unittest.TestCase):
    """sporulation held a flag, prose about spore type, and cell arrangement."""

    def test_integer_and_word_encodings_of_the_flag_agree(self):
        self.assertEqual(split_sporulation("0"), ("no", ""))
        self.assertEqual(split_sporulation("1"), ("yes", ""))
        self.assertEqual(split_sporulation("Nonsporulating"), ("no", ""))
        self.assertEqual(split_sporulation("Sporulating"), ("yes", ""))

    def test_prose_about_spore_type_reduces_to_the_flag(self):
        self.assertEqual(
            split_sporulation("Sporulating (ascospores); reproduction also by budding"),
            ("yes", ""))
        self.assertEqual(split_sporulation("No; Reproduction by unipolar budding"),
                         ("no", ""))

    def test_a_known_typo_is_still_recognised(self):
        self.assertEqual(split_sporulation("Spotulating"), ("yes", ""))

    def test_cell_arrangement_text_is_routed_to_the_other_column(self):
        for value in ("Singles", "Pairs - Singles", "Chains - Pairs - Singles",
                      "Clusters - Singles"):
            spore, arrangement = split_sporulation(value)
            self.assertEqual(spore, "", value)
            self.assertEqual(arrangement, value, value)

    def test_both_delimiter_styles_are_handled(self):
        """'Pairs-Singles' splits on a bare hyphen; 'Singles - V-shaped forms' does not.

        Splitting the second on a bare hyphen yields 'v' and 'shaped forms',
        which matches nothing, so it was being left in the wrong column.
        """
        self.assertEqual(split_sporulation("Pairs-Singles")[1], "Pairs-Singles")
        self.assertEqual(split_sporulation("Singles - V-shaped forms")[1],
                         "Singles - V-shaped forms")

    def test_an_unrecognised_value_is_left_alone_rather_than_guessed_at(self):
        """('', '') means 'do not touch this cell', not 'clear it'."""
        self.assertEqual(split_sporulation("Septated hyphae"), ("", ""))
        self.assertEqual(split_sporulation("something nobody anticipated"), ("", ""))


class StudyPlatformDerivationTests(unittest.TestCase):
    """studies.data_type is derived from the samples so the stores cannot differ.

    The mirror held four vocabularies for two platforms -- 'AMPLICON' beside
    '16S rRNA', 'Metagenomics' beside 'shotgun metagenomics' -- where MySQL held
    '16S'/'mNGS', so `WHERE data_type = '16S'` matched 236 studies in one store
    and none in the other.
    """

    def test_one_platform_however_it_is_spelled(self):
        for spellings in (("AMPLICON",), ("16S",), ("16S rRNA",),
                          ("AMPLICON", "16S rRNA", "16S")):
            self.assertEqual(derive_study_platform(spellings), "16S", spellings)
        for spellings in (("Metagenomics",), ("mNGS",),
                          ("shotgun metagenomics", "Metagenomics")):
            self.assertEqual(derive_study_platform(spellings), "mNGS", spellings)

    def test_two_real_platforms_become_mixed(self):
        self.assertEqual(derive_study_platform(["AMPLICON", "Metagenomics"]), "mixed")
        self.assertEqual(derive_study_platform(["16S", "mNGS"]), "mixed")

    def test_nothing_to_derive_from_returns_empty_meaning_leave_it_alone(self):
        """An association-only study has no samples; 'meta-analysis' is true there."""
        for empty in ([], [None], [""], [None, ""]):
            self.assertEqual(derive_study_platform(empty), "", empty)

    def test_mixed_does_not_survive_as_an_input_value(self):
        """'mixed' is a summary, so feeding it back must not make it a platform."""
        self.assertEqual(derive_study_platform(["mixed"]), "mixed")
        self.assertEqual(derive_study_platform(["mixed", "AMPLICON"]), "mixed")

    def test_a_sample_may_not_inherit_a_summary_or_a_study_design(self):
        self.assertIn("mixed", NOT_A_PLATFORM)
        self.assertIn("meta-analysis", NOT_A_PLATFORM)
        self.assertIn("", NOT_A_PLATFORM)
        self.assertNotIn("16S", NOT_A_PLATFORM)
        self.assertNotIn("mNGS", NOT_A_PLATFORM)


class AttributeNumericTrustTests(unittest.TestCase):
    """value_numeric must hold a number you can do arithmetic on.

    Two classes broke that while looking fine: a number fabricated from a string
    that is not one, and a category code whose meaning lives in the key name.
    """

    def test_a_number_is_kept(self):
        for attribute, value in (("Creatinine", "1.42"), ("eGFR", "-1.5"),
                                 ("Brinkman Index", "2400"), ("Timepoint", " 3 ")):
            self.assertTrue(numeric_value_is_trustworthy(attribute, value),
                            (attribute, value))

    def test_a_fabricated_number_is_refused(self):
        """'10007539_0' was stored as 100075390 -- and '1000753_90' collides with it."""
        self.assertFalse(numeric_value_is_trustworthy("unique_id", "10007539_0"))
        self.assertFalse(numeric_value_is_trustworthy("unique_id", "1000753_90"))

    def test_an_interval_is_not_a_number(self):
        for value in ("0-1", "2-4", "3+", "1-", "AADF", "IIIB", ""):
            self.assertFalse(numeric_value_is_trustworthy("anything", value), value)

    def test_a_code_declared_in_the_key_name_is_refused(self):
        """1 means NO here and YES under `Antibiotics`; pooling them inverts rows."""
        for attribute in ("Antibiotic_in_the_past_month._1_no_2_yes",
                          "X_probiotic_in_the_past_2_weeks_1_no_2_yes",
                          "smoker_yes_1_no_0"):
            self.assertFalse(numeric_value_is_trustworthy(attribute, "1"), attribute)

    def test_a_key_merely_mentioning_a_range_is_not_treated_as_coded(self):
        """'IBS_symptom_severity_1_500' is a scale, not a yes/no code."""
        self.assertTrue(numeric_value_is_trustworthy("IBS_symptom_severity_1_500",
                                                     "320"))
        self.assertTrue(numeric_value_is_trustworthy("Effect_on_life_in_general_0_100",
                                                     "45"))

    def test_a_plain_number_under_a_mixed_key_is_still_kept(self):
        """Stage = 0 is correct; what is missing is the denominator, not the value."""
        self.assertTrue(numeric_value_is_trustworthy("Stage", "0"))
        self.assertTrue(numeric_value_is_trustworthy("Antibiotics", "1"))


class MetaboliteScaleTests(unittest.TestCase):
    """sample_metabolite_levels mixes scales that differ in kind.

    Two of 13 cohorts arrive log- or z-transformed and hold 70,067 negative
    levels; the rest are non-negative with maxima spanning nine orders of
    magnitude. The hazard is the zero: absence on one scale, mid-distribution on
    the other.
    """

    def test_a_negative_level_proves_an_upstream_transform(self):
        """A concentration cannot be negative, so this is positive evidence."""
        self.assertEqual(derive_metabolite_scale(-6.725, 7.6), "centred")
        self.assertEqual(derive_metabolite_scale(-2.481, 3.3), "centred")

    def test_no_negatives_claims_only_non_negativity(self):
        """Not a unit claim: the source records units nowhere."""
        for lo, hi in ((0.0, 665180.6), (0.0, 266.7), (1.0, 3297391681.0)):
            self.assertEqual(derive_metabolite_scale(lo, hi), "non_negative")

    def test_a_zero_means_absence_only_on_the_non_negative_scale(self):
        self.assertTrue(zero_means_absence(derive_metabolite_scale(0.0, 665180.6)))
        self.assertFalse(zero_means_absence(derive_metabolite_scale(-6.725, 7.6)))

    def test_nothing_to_derive_from_claims_nothing(self):
        self.assertEqual(derive_metabolite_scale(None, None), "")
        self.assertFalse(zero_means_absence(""))

    def test_a_zero_minimum_is_non_negative_not_centred(self):
        """Boundary: 0.0 is not negative, so it must not be called centred."""
        self.assertEqual(derive_metabolite_scale(0.0, 1.0), "non_negative")


class NcbiRankTests(unittest.TestCase):
    """taxonomic_rank is the abundance layer; ncbi_rank is what the organism is."""

    def test_a_family_or_coarser_in_the_genus_layer_is_flagged(self):
        for rank in ("family", "order", "class", "phylum", "infraclass",
                     "infraorder", "superkingdom"):
            self.assertTrue(coarser_than_genus(rank), rank)

    def test_genus_and_anything_below_it_is_not_flagged(self):
        for rank in ("genus", "subgenus", "species", "subspecies", "strain",
                     "species group", "isolate"):
            self.assertFalse(coarser_than_genus(rank), rank)

    def test_an_uninformative_rank_asserts_nothing(self):
        """'no rank' and 'clade' must not be read as a disagreement."""
        for rank in ("no rank", "clade", "", None):
            self.assertFalse(coarser_than_genus(rank), rank)

    def test_case_and_whitespace_do_not_change_the_answer(self):
        self.assertTrue(coarser_than_genus("  Family "))
        self.assertFalse(coarser_than_genus(" GENUS "))
