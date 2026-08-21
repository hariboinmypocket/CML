#!/usr/bin/env python3
"""Harvest per-sample age/sex/BMI for PRJNA562327 (Ren et al. 2020, CKD).

Source: Supporting Excel 1 of https://doi.org/10.1002/advs.202001936
489 participants with age, gender and BMI and no missing values.

The supplement carries no SRR/SAMN column, so rows are joined to runs through
the SRA library_name. Three series exist and each is handled explicitly:

  CKD_n    -> Mobio_ZZ_CKD_n   only CKD library set; 159 distinct ids in 1-160
  HzCKD-n  -> hzCKD_n          only hzCKD library set; 57 ids in 1-58
  HC_n     -> zzHC_n           INFERRED, see below

The HC join is an inference, not a documented mapping. Two HC library sets
exist (zzHC_1-320 and Mobio_ZZ_HC_1-160), but the supplement's 273 HC ids form
one duplicate-free series spanning 1-320, and only zzHC spans that far, so the
series cannot be hosted by Mobio_ZZ_HC. The risk it carries is a permutation
error confined to HC rows; it cannot be detected by comparing distributions,
because any reordering leaves the distribution unchanged. Every written row is
tagged in the audit file with the series used so the HC subset can be reverted
on its own if the mapping is ever disproved.

Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import xlrd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from gutdb.config import Settings
from gutdb.db import connect

SUPP = Path("/private/tmp/claude-501/-Users-lucyshaa-Desktop-CML/"
            "0a6c9862-7ad6-4ee9-a30e-151a04129642/scratchpad/supp/ADVS-7-2001936-s003.xls")
ENA_TSV = Path("/tmp/ck_dates.tsv")
AUDIT = Path("data/ckd_study_demographics_audit.csv")
DRY_RUN = "--apply" not in sys.argv

SERIES = {"CKD": "Mobio_ZZ_CKD_{}", "HZCKD": "hzCKD_{}", "HC": "zzHC_{}"}


def norm_header(cell: str) -> str:
    """Sheets differ in punctuation width and column naming."""
    text = str(cell).strip().lower()
    text = text.replace("（", "(").replace("）", ")")
    if text.startswith("age"):
        return "age"
    if text.startswith("gender") or text.startswith("sex"):
        return "sex"
    if text.startswith("bmi"):
        return "bmi"
    if text in ("sample", "number"):
        return "sample"
    return text


def parse_supplement() -> list[dict[str, object]]:
    wb = xlrd.open_workbook(SUPP)
    out: list[dict[str, object]] = []
    for sheet in wb.sheets():
        header = [norm_header(c.value) for c in sheet.row(0)]
        idx = {name: header.index(name) for name in ("sample", "age", "sex", "bmi")
               if name in header}
        if "sample" not in idx:
            continue
        for r in range(1, sheet.nrows):
            raw_id = str(sheet.cell_value(r, idx["sample"])).strip()
            m = re.match(r"^(HzCKD|CKD|HC)[-_](\d+)$", raw_id, re.I)
            if not m:
                continue
            series = m.group(1).upper()
            number = int(m.group(2))

            def cell(name: str):
                if name not in idx:
                    return None
                value = str(sheet.cell_value(r, idx[name])).strip()
                return value or None

            age_raw, sex_raw, bmi_raw = cell("age"), cell("sex"), cell("bmi")
            try:
                age = float(age_raw) if age_raw not in (None, "NA") else None
            except ValueError:
                age = None
            sex = None
            if sex_raw:
                folded = sex_raw.strip().casefold()
                sex = "male" if folded == "male" else ("female" if folded == "female" else None)
            try:
                bmi = float(bmi_raw) if bmi_raw not in (None, "NA") else None
            except ValueError:
                bmi = None
            if age is not None and not (0 <= age <= 122):
                age = None
            if bmi is not None and not (10 <= bmi <= 80):
                bmi = None
            out.append({"raw_id": raw_id, "series": series, "number": number,
                        "age": age, "sex": sex, "bmi": bmi,
                        "library_name": SERIES[series].format(number)})
    return out


def main() -> int:
    rows = parse_supplement()
    print(f"supplement rows parsed: {len(rows)}")

    lib_to_run: dict[str, str] = {}
    with ENA_TSV.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            lib = (row.get("library_name") or "").strip()
            if lib:
                lib_to_run[lib] = row["run_accession"]
    print(f"SRA library_name -> run map: {len(lib_to_run)}")

    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    cur.execute("SELECT run_accession, id, sex, age_years, bmi FROM samples")
    ours = {r[0]: {"id": r[1], "sex": r[2], "age": r[3], "bmi": r[4]} for r in cur.fetchall()}
    cur.close()

    audit: list[dict[str, object]] = []
    stats = {"matched": 0, "no_library": 0, "not_in_db": 0,
             "sex_filled": 0, "age_filled": 0, "bmi_filled": 0,
             "sex_conflict": 0, "age_conflict": 0}

    for row in rows:
        run = lib_to_run.get(str(row["library_name"]))
        if not run:
            stats["no_library"] += 1
            continue
        current = ours.get(run)
        if not current:
            stats["not_in_db"] += 1
            continue
        stats["matched"] += 1

        set_sex = set_age = set_bmi = None
        if row["sex"]:
            if current["sex"] is None:
                set_sex = row["sex"]; stats["sex_filled"] += 1
            elif str(current["sex"]).casefold() != row["sex"]:
                stats["sex_conflict"] += 1
        if row["age"] is not None:
            if current["age"] is None:
                set_age = row["age"]; stats["age_filled"] += 1
            elif abs(float(current["age"]) - float(row["age"])) > 1.0:
                stats["age_conflict"] += 1
        if row["bmi"] is not None and current["bmi"] is None:
            set_bmi = row["bmi"]; stats["bmi_filled"] += 1

        if not DRY_RUN and (set_sex or set_age is not None or set_bmi is not None):
            upd = conn.cursor()
            upd.execute(
                """UPDATE samples SET sex = COALESCE(sex, %s),
                       age_years = COALESCE(age_years, %s), bmi = COALESCE(bmi, %s)
                   WHERE id = %s""",
                (set_sex, set_age, set_bmi, current["id"]),
            )
            upd.close()

        audit.append({"supplement_id": row["raw_id"], "series": row["series"],
                      "library_name": row["library_name"], "run_accession": run,
                      "join_basis": "documented" if row["series"] != "HC" else "inferred_single_series",
                      "age": row["age"], "sex": row["sex"], "bmi": row["bmi"],
                      "filled_age": set_age, "filled_sex": set_sex, "filled_bmi": set_bmi})

    if not DRY_RUN:
        conn.commit()
        AUDIT.parent.mkdir(parents=True, exist_ok=True)
        with AUDIT.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(audit[0]))
            writer.writeheader()
            writer.writerows(audit)
        print(f"audit written -> {AUDIT}")
    else:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")

    for key, value in stats.items():
        print(f"  {key:14} {value}")
    by_series: dict[str, int] = {}
    for a in audit:
        if a["filled_age"] is not None or a["filled_sex"] or a["filled_bmi"]:
            by_series[str(a["series"])] = by_series.get(str(a["series"]), 0) + 1
    print(f"  rows filled by series: {by_series}")

    cur = conn.cursor(); cur.close(); conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
