#!/usr/bin/env python3
"""Map the source's 362 attribute names onto one concept each, with decoding.

sample_attributes keeps every key exactly as the source wrote it, which is right
-- a value can always be traced back to the file and column it came from -- and
useless for asking a question. Antibiotic exposure is spread over seven yes/no
keys, so `WHERE attribute = 'Antibiotics'` finds 518 of the 1,423 rows that
answer it.

THE VALUE MAP IS THE POINT, NOT THE NAME MAP. `1` means YES under `Antibiotics`
and NO under `Antibiotic_in_the_past_month._1_no_2_yes`. A name-only alias would
let those pool and invert 343 rows -- the same family of error as pooling an LDA
score without its comparison order. So every boolean key carries its own
raw -> canonical decoding and nothing is inferred from the raw value.

`concept` is the QUESTION, not the topic. Antibiotics carry three -- was there
exposure, which agent, and when -- and one `antibiotic` concept would put drug
names and booleans in the same column.

WHAT IS DELIBERATELY NOT MAPPED, because a shared word is not a shared concept.
Each of these was a false cluster a name heuristic produced and reading the
values rejected:

  `steroids` and `mesalamine` are distinct drugs, not members of an
  immunosuppressant concept. Only `immunosuppressant` and
  `Immunosuppressive_medication` are the same question.

  `Breast cancer` is a diagnosis. A '%breast%' match pulls it in beside the
  infant-feeding keys.

  `GADA_Level`, `IA2A_Level`, `IAA_Level`, `AAB_Level`, `ICA_Level` are five
  different autoantibodies, not one 'level'.

Unmapped keys are the normal case, not a backlog: most of the 362 are a single
cohort's own question and need no alias. The view returns them with concept NULL
rather than hiding them, and this script reports how many.

Idempotent. Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect

DRY_RUN = "--apply" not in sys.argv

# (attribute, concept, value_kind, notes)
ALIASES = [
    # --- antibiotics: three questions, not one topic ----------------------
    ("Antibiotics", "antibiotic_exposure", "boolean",
     "two encodings in one key: no/yes and 0/1"),
    ("antibiotic", "antibiotic_exposure", "boolean", None),
    ("Antibiotic_in_the_past_month._1_no_2_yes", "antibiotic_exposure", "boolean",
     "1 means NO here, the opposite of Antibiotics"),
    ("Antibiotics_Last3months", "antibiotic_exposure", "boolean",
     "the yes carries a timing detail the boolean drops"),
    ("AbxExposureAbsolute", "antibiotic_exposure", "boolean", None),
    ("PostAbxExposure", "antibiotic_exposure", "boolean",
     "exposure after collection, so not interchangeable with the others"),
    ("maternal_antibiotics", "antibiotic_exposure", "boolean",
     "the MOTHER's exposure, not the sampled infant's"),
    ("antibiotic_type", "antibiotic_agent", "agent", None),
    ("which_abx", "antibiotic_agent", "agent", None),
    ("Antibiotic_used", "antibiotic_agent", "agent", "agent plus dose and duration"),
    ("AbxAtCollection", "antibiotic_agent", "agent", "'no_abx' means no exposure"),
    ("AbxPreCollection", "antibiotic_agent", "agent", None),
    ("Antibiotic_start_date", "antibiotic_start_date", "date", None),
    # --- probiotics -------------------------------------------------------
    ("X_probiotic_in_the_past_2_weeks_1_no_2_yes", "probiotic_exposure", "boolean",
     "1 means NO"),
    ("Probiotics", "probiotic_product", "free_text", None),
    ("Probiotic_used", "probiotic_product", "free_text", None),
    # --- smoking: status and exposure are different measurements ----------
    ("SmokingStatus", "smoking_status", "category", None),
    ("Tobacco_Use", "smoking_status", "category", "mixes status with frequency"),
    ("Chem ID / Smoking history", "smoking_status", "category", None),
    ("Brinkman Index", "smoking_exposure_index", "quantity",
     "cigarettes per day times years smoked; 0 for a non-smoker"),
    # --- alcohol ----------------------------------------------------------
    ("AlcoholConsumption", "alcohol_status", "boolean", None),
    ("Drinking Status", "alcohol_status", "category",
     "three states plus Unknown, so not reducible to the boolean"),
    ("Alcohol", "alcohol_amount", "quantity", "cohort's own units, not comparable"),
    ("Alcohol.drink.week.", "alcohol_amount", "interval",
     "point answers beside ranges; the numeric subset is biased"),
    # --- infant feeding and delivery --------------------------------------
    ("feeding_type", "infant_feeding", "category", "BM / BM,F"),
    ("Infant_Formula", "infant_feeding", "boolean", None),
    ("Breast_Feeding", "infant_feeding", "boolean", None),
    ("Delivery_Route", "delivery_mode", "category", None),
    ("delivery_mode", "delivery_mode", "category", None),
    # --- immunosuppression: only these two are one question ---------------
    ("immunosuppressant", "immunosuppressant_use", "boolean", None),
    ("Immunosuppressive_medication", "immunosuppressant_use", "boolean", None),
]

# (attribute, raw_value, canonical_value). Every boolean key appears here; a
# boolean whose value is absent from this table surfaces as boolean_undecoded in
# v_sample_attributes_canonical rather than guessing.
VALUE_ALIASES = [
    ("Antibiotics", "yes", "yes"), ("Antibiotics", "1", "yes"),
    ("Antibiotics", "no", "no"), ("Antibiotics", "0", "no"),
    ("antibiotic", "Yes", "yes"), ("antibiotic", "No", "no"),
    # The inversion this table exists for.
    ("Antibiotic_in_the_past_month._1_no_2_yes", "1", "no"),
    ("Antibiotic_in_the_past_month._1_no_2_yes", "2", "yes"),
    ("Antibiotics_Last3months", "No", "no"),
    ("Antibiotics_Last3months", "Yes, 2 months ago", "yes"),
    ("AbxExposureAbsolute", "TRUE", "yes"), ("AbxExposureAbsolute", "FALSE", "no"),
    ("PostAbxExposure", "TRUE", "yes"), ("PostAbxExposure", "FALSE", "no"),
    ("maternal_antibiotics", "1", "yes"), ("maternal_antibiotics", "0", "no"),
    ("X_probiotic_in_the_past_2_weeks_1_no_2_yes", "1", "no"),
    ("X_probiotic_in_the_past_2_weeks_1_no_2_yes", "2", "yes"),
    ("AlcoholConsumption", "Yes", "yes"), ("AlcoholConsumption", "No", "no"),
    ("immunosuppressant", "Yes", "yes"), ("immunosuppressant", "No", "no"),
    ("Immunosuppressive_medication", "No", "no"),
    ("Infant_Formula", "TRUE", "yes"), ("Infant_Formula", "FALSE", "no"),
    ("Breast_Feeding", "No", "no"),
    # Categories, decoded where two cohorts spell one state differently.
    ("SmokingStatus", "Smoking", "current"),
    ("SmokingStatus", "NotSmoking", "never_or_not_current"),
    ("SmokingStatus", "StoppedSmoking", "former"),
    ("SmokingStatus", "Unknown", "unknown"),
    ("Chem ID / Smoking history", "Smoke", "current"),
    ("Chem ID / Smoking history", "Non-smoke", "never_or_not_current"),
    ("Chem ID / Smoking history", "Missing", "unknown"),
    ("Drinking Status", "Drink", "current"),
    ("Drinking Status", "Not Drinking", "never_or_not_current"),
    ("Drinking Status", "Stop Drinking", "former"),
    ("Drinking Status", "Unknown", "unknown"),
    ("Delivery_Route", "cesarian", "caesarean"),
    ("Delivery_Route", "vaginal", "vaginal"),
    ("delivery_mode", "c-section", "caesarean"),
    ("delivery_mode", "vaginal", "vaginal"),
]


def main() -> int:
    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()

    cur.execute("SELECT DISTINCT attribute FROM sample_attributes")
    present = {r[0] for r in cur.fetchall()}
    missing = [a for a, *_ in ALIASES if a not in present]
    if missing:
        print(f"keys in the seed that are NOT in sample_attributes ({len(missing)}):")
        for a in missing:
            print(f"   {a}")
        print("   -- a typo here would silently alias nothing, so this is fatal")
        conn.close()
        return 1
    print(f"all {len(ALIASES)} seeded keys exist in sample_attributes")

    orphan_values = [(a, v) for a, v, _ in VALUE_ALIASES
                     if a not in {x[0] for x in ALIASES}]
    if orphan_values:
        print(f"value aliases with no key alias: {orphan_values}")
        conn.close()
        return 1

    if not DRY_RUN:
        cur.executemany(
            """INSERT INTO sample_attribute_aliases (attribute, concept,
                                                     value_kind, notes)
               VALUES (%s,%s,%s,%s)
               ON DUPLICATE KEY UPDATE concept = VALUES(concept),
                   value_kind = VALUES(value_kind), notes = VALUES(notes)""",
            ALIASES)
        cur.executemany(
            """INSERT INTO sample_attribute_value_aliases
                 (attribute, raw_value, canonical_value)
               VALUES (%s,%s,%s)
               ON DUPLICATE KEY UPDATE canonical_value = VALUES(canonical_value)""",
            VALUE_ALIASES)
        conn.commit()

    cur.execute("""SELECT al.concept, COUNT(DISTINCT al.attribute) keys_,
                          COUNT(*) rows_
                   FROM sample_attribute_aliases al
                   JOIN sample_attributes a ON a.attribute = al.attribute
                   GROUP BY al.concept ORDER BY rows_ DESC""")
    print(f"\n{'concept':26} {'keys':>5} {'rows':>6}")
    for concept, keys_, rows_ in cur.fetchall():
        print(f"{concept:26} {keys_:>5} {rows_:>6}")

    cur.execute("""SELECT COUNT(*) FROM v_sample_attributes_canonical
                   WHERE boolean_undecoded = 1""")
    undecoded = cur.fetchone()[0]
    cur.execute("""SELECT COUNT(DISTINCT attribute) FROM sample_attributes
                   WHERE attribute NOT IN
                     (SELECT attribute FROM sample_attribute_aliases)""")
    unmapped = cur.fetchone()[0]
    print(f"\nboolean values with no decoding (must be 0): {undecoded}")
    print(f"keys left unmapped: {unmapped} of 362 -- the normal case, each a "
          f"single cohort's own question")
    conn.close()

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
