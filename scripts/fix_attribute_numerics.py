#!/usr/bin/env python3
"""Make sample_attributes.value_numeric mean what the schema says it means.

The column is documented as "populated only when `value` parses as a number ...
lets a numeric trait be filtered and aggregated without casting strings". Both
halves of that were false.

NOT A PARSE OF THE VALUE (434 rows, one key). `unique_id` holds '10007539_0'
with value_numeric = 100075390: the underscore was dropped and the digits run
together. That is a fabricated number, and it can collide two distinct ids --
'10007539_0' and '1000753_90' both reduce to 100075390.

A CODE, NOT A QUANTITY (687 rows, two keys). `Antibiotic_in_the_past_month._1_no_2_yes`
and `X_probiotic_in_the_past_2_weeks_1_no_2_yes` store 1 and 2 as labels. The
mean of those is the mean of labels, and the encoding INVERTS against its
neighbour: 1 means no here, while under `Antibiotics` 1 means yes. Pooling the
two silently flips 343 rows -- the same family of error as the LDA sign
ambiguity the schema already guards against.

Both classes are cleared, because in both the stored number does not mean the
quantity it appears to.

A THIRD CLASS IS DELIBERATELY LEFT ALONE: 13 keys mix numbers with words, and
there the number is correct for whichever rows have one. What is missing is the
denominator.

    Stage                 27/150 numeric   0 beside I, IIA, IIIB  -> AVG is 0
    Antibiotics           74/518           0, 1 beside no, yes
    Race                 131/371           0, 1 beside White, Black, Hispanic
    Alcohol.drink.week.   98/161           point values beside ranges '0-1', '2-4'
    BowelMovement_Day     63/162           same
    Timepoint            508/518           the 10 others say 'Flare'

`Stage` is the one to look at: AVG(value_numeric) over it is the mean of 27
zeros, which is 0, for a column whose real values run I to IVB. Deleting those
27 zeros would not help -- they are correct -- so this script exposes the
denominator instead, which is the same move the abundance matrix needed for
absence. v_sample_attribute_keys carries n_rows, n_numeric and the fraction, so
a mean can be read against what it was taken over.

The ranges are worth naming separately: under `Alcohol.drink.week.` and
`BowelMovement_Day` the non-numeric values are intervals, so the numeric subset
is not a random sample of the key -- it is exactly the respondents who gave a
point answer.

MySQL only; the SQLite mirror has no sample_attributes table.
Idempotent. Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.transform import CODED_IN_NAME, PLAIN_NUMBER, numeric_value_is_trustworthy

DRY_RUN = "--apply" not in sys.argv
AUDIT = Path("data/attribute_numeric_audit.csv")
AUDIT_FIELDS = ("attribute", "value", "value_numeric", "n_rows", "reason")


def main() -> int:
    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    rows: list[dict] = []
    counts: Counter = Counter()

    cur.execute("""SELECT attribute, value, value_numeric, COUNT(*)
                   FROM sample_attributes
                   WHERE value_numeric IS NOT NULL
                   GROUP BY attribute, value, value_numeric""")
    clear: list[tuple[str, str]] = []
    for attribute, value, numeric, n in cur.fetchall():
        if numeric_value_is_trustworthy(attribute, value):
            counts["kept"] += n
            continue
        reason = ("not_a_number" if not PLAIN_NUMBER.match(str(value).strip())
                  else "code_declared_in_key_name")
        counts[reason] += n
        rows.append({"attribute": attribute, "value": value,
                     "value_numeric": numeric, "n_rows": n, "reason": reason})
        clear.append((attribute, value))

    by_reason = Counter(r["reason"] for r in rows)
    print("value_numeric cells that do not mean the quantity they appear to:")
    for reason, n_keys in by_reason.items():
        print(f"   {reason:26} {counts[reason]:>5} rows across "
              f"{len({r['attribute'] for r in rows if r['reason'] == reason})} keys")
    print(f"   {'kept':26} {counts['kept']:>5} rows")

    for attribute in sorted({r["attribute"] for r in rows}):
        n = sum(r["n_rows"] for r in rows if r["attribute"] == attribute)
        why = {r["reason"] for r in rows if r["attribute"] == attribute}
        print(f"      {attribute[:52]:52} {n:>4}  {','.join(sorted(why))}")

    if not DRY_RUN and clear:
        cur.executemany("""UPDATE sample_attributes SET value_numeric = NULL
                           WHERE attribute = %s AND value = %s""", clear)
        conn.commit()
        print(f"\ncleared {cur.rowcount if cur.rowcount > 0 else len(clear)} "
              f"(attribute, value) pairs")

    print("\nkeys that mix numbers with words -- left alone, denominator exposed:")
    cur.execute("""SELECT attribute, COUNT(*) n, SUM(value_numeric IS NOT NULL) n_num
                   FROM sample_attributes GROUP BY attribute
                   HAVING n_num > 0 AND n_num < n ORDER BY n_num DESC""")
    for attribute, n, n_num in cur.fetchall():
        print(f"   {attribute[:40]:40} {n_num:>4}/{n:<4} "
              f"({n_num / n:.0%} carry a number)")

    AUDIT.parent.mkdir(parents=True, exist_ok=True)
    with AUDIT.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=AUDIT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {AUDIT} ({len(rows)} rows)")
    conn.close()

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
