import unittest
import csv
from pathlib import Path

from gutdb.transform import (
    association_direction,
    collapse_energy_modes,
    collapse_food_sources,
    derive_energy_mode,
    derive_primary_food_source,
    extract_species_epithet,
    normalize_genus,
    parse_scientific_name,
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


if __name__ == "__main__":
    unittest.main()
