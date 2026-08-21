#!/usr/bin/env python3
"""Check how much of our sample set curatedMetagenomicData could annotate.

cMD's curation repo carries demographics extracted from the papers by human
curators, keyed to NCBI run accessions. That makes it independent of BioSample
and ENA, which is why it can hold age/sex for samples where those are empty.

Read-only: reports overlap, writes nothing to the database.
"""
from __future__ import annotations

import csv
import io
import sys
import time
from collections import defaultdict
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from gutdb.config import Settings
from gutdb.db import connect

API = "https://api.github.com/repos/waldronlab/curatedMetagenomicDataCuration/contents/inst/curated"
RAW = "https://raw.githubusercontent.com/waldronlab/curatedMetagenomicDataCuration/master/inst/curated"

session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=Retry(
    total=4, backoff_factor=1.5, status_forcelist=(429, 500, 502, 503, 504))))


def study_list() -> list[str]:
    r = session.get(API, timeout=60)
    r.raise_for_status()
    return [x["name"] for x in r.json() if x.get("type") == "dir"]


def fetch_samples(study: str) -> list[dict[str, str]]:
    for suffix in ("_sample.tsv", "_metadata.tsv"):
        url = f"{RAW}/{study}/{study}{suffix}"
        try:
            r = session.get(url, timeout=90)
        except Exception:
            continue
        if r.status_code == 200 and r.text.strip():
            return list(csv.DictReader(io.StringIO(r.text), delimiter="\t"))
    return []


def main() -> int:
    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    cur.execute("SELECT run_accession, sex, age_years, bmi FROM samples")
    ours = {r[0]: {"sex": r[1], "age": r[2], "bmi": r[3]} for r in cur.fetchall()}
    cur.close()
    conn.close()
    print(f"our samples: {len(ours):,}", flush=True)

    studies = study_list()
    print(f"cMD studies: {len(studies)}", flush=True)

    matched: dict[str, dict[str, str]] = {}
    per_study_hits: dict[str, int] = defaultdict(int)
    empty_studies = 0

    for i, study in enumerate(studies, 1):
        rows = fetch_samples(study)
        if not rows:
            empty_studies += 1
            continue
        for row in rows:
            acc = (row.get("ncbi_accession") or "").strip()
            if not acc:
                continue
            # a cMD row can list several runs for one sample
            for token in acc.replace(";", ",").split(","):
                token = token.strip()
                if token in ours:
                    matched[token] = row
                    per_study_hits[study] += 1
        if i % 25 == 0:
            print(f"  ...{i}/{len(studies)} studies, {len(matched):,} matches so far", flush=True)
        time.sleep(0.05)

    print(f"\nstudies with no usable file: {empty_studies}")
    print(f"our samples found in cMD: {len(matched):,}")

    fills_sex = fills_age = fills_bmi = 0
    conflicts_sex = conflicts_age = 0
    non_year_units: dict[str, int] = defaultdict(int)
    for acc, row in matched.items():
        cur_state = ours[acc]
        # cMD names the column "sex"; "gender" is checked only as a fallback
        sex_val = (row.get("sex") or row.get("gender") or "").strip()
        age_val = (row.get("age") or "").strip()
        unit = (row.get("age_unit") or "").strip()
        if unit and unit.lower() not in ("year", "years"):
            non_year_units[unit] += 1
        if sex_val:
            if cur_state["sex"] is None:
                fills_sex += 1
            elif sex_val.lower() != str(cur_state["sex"]).lower():
                conflicts_sex += 1
        if age_val:
            if cur_state["age"] is None:
                fills_age += 1
            else:
                try:
                    if abs(float(age_val) - float(cur_state["age"])) > 1.0:
                        conflicts_age += 1
                except ValueError:
                    pass
        if cur_state["bmi"] is None and (row.get("bmi") or "").strip():
            fills_bmi += 1

    print(f"\nWOULD FILL a currently-empty field:")
    print(f"  sex : {fills_sex:,}")
    print(f"  age : {fills_age:,}")
    print(f"  bmi : {fills_bmi:,}")
    print(f"\nDISAGREEMENTS with values we already hold:")
    print(f"  sex : {conflicts_sex:,}")
    print(f"  age : {conflicts_age:,} (differ by >1 year)")
    if non_year_units:
        print(f"\nnon-year age units present: {dict(non_year_units)}")

    print("\ntop contributing cMD studies:")
    for study, n in sorted(per_study_hits.items(), key=lambda kv: -kv[1])[:12]:
        print(f"  {n:5d}  {study}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
