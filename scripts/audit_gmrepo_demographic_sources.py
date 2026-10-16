#!/usr/bin/env python3
"""Trace filled demographic cells to surviving historical parsing audits.

This report distinguishes an audit match from independently verifying a primary
record. It does not change the sample CSV or database. URLs point to records
that can be checked; their presence alone is not proof of a matching value.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from pathlib import Path
from urllib.parse import urlencode


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CHANGES = ROOT / "data/archive/gmrepo_samples_cleanup_20261016T152817861919Z/field_changes.csv"


def read(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def same(field, a, b):
    if not a or not b:
        return False
    if field == "sex":
        return a.casefold() == b.casefold()
    return abs(float(a) - float(b)) < 1e-6


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--changes", type=Path, default=DEFAULT_CHANGES)
    parser.add_argument("--output", type=Path, default=ROOT / "data/gmrepo_demographic_provenance.csv")
    args = parser.parse_args()
    fills = [r for r in read(args.changes)
             if r["field"] in {"sex", "age_years", "bmi"} and not r["before"] and r["after"]]
    samples = {r["run_id"]: r for r in read(ROOT / "data/gmrepo_project_samples.csv")}
    index = defaultdict(list)
    sources = [
        ("NCBI_SRA_BioSample", "ncbi_biosample_demographics.csv",
         {"sex": "sex", "age_years": "age_years", "bmi": "bmi"}),
        ("curatedMetagenomicData", "cmd_demographics_audit.csv",
         {"sex": "cmd_sex", "age_years": "cmd_age_years", "bmi": "cmd_bmi"}),
        ("Ren_2020_CKD_supplement", "ckd_study_demographics_audit.csv",
         {"sex": "sex", "age_years": "age_years", "bmi": "bmi"}),
    ]
    for name, filename, mapping in sources:
        path = ROOT / "data" / filename
        if path.exists():
            for record in read(path):
                if name == "Ren_2020_CKD_supplement" and "age_years" not in record:
                    # Also accept original, archived audits before column cleanup.
                    record["age_years"] = record.get("age", "")
                index[record["run_accession"]].append((name, path, mapping, record))
    output = []
    flags = Counter()
    for change in fills:
        run, field, value = change["run_id"], change["field"], change["after"]
        # Historical fills that were withdrawn are not active demographic data.
        if run not in samples or not same(field, samples[run].get(field), value):
            continue
        names, audits, links, notes, warnings = [], [], [], [], set()
        for name, path, mapping, record in index[run]:
            if not same(field, record.get(mapping[field]), value):
                continue
            names.append(name)
            audits.append(str(path))
            if name == "NCBI_SRA_BioSample":
                links.append("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" +
                             urlencode(dict(db="sra", id=run, rettype="full", retmode="xml")))
                if field == "age_years" and record.get("age_note"):
                    notes.append(record["age_note"])
                    if "censored" in record["age_note"]:
                        warnings.add("censored_age_not_exact")
            elif name == "curatedMetagenomicData":
                study = record["cmd_study"]
                links.append("https://github.com/waldronlab/curatedMetagenomicDataCuration/"
                             "tree/master/inst/curated/" + study)
                notes.append("cMD study: " + study)
            else:
                links.append("https://pmc.ncbi.nlm.nih.gov/articles/PMC7578882/")
                notes.append("Supporting Excel 1; supplement_id=" + record["supplement_id"] +
                             "; library_name=" + record["library_name"] +
                             "; join_basis=" + record["join_basis"])
                if record["join_basis"] not in {"documented", "source_alias_matches_after_separator_normalization"}:
                    warnings.add("inferred_sample_link_requires_verification")
        if not names:
            warnings.add("no_matching_surviving_source_audit")
        for warning in warnings:
            flags[warning] += 1
        sample = samples[run]
        output.append(dict(run_id=run, sample_id=sample["sample_id"],
                           project_id=sample["project_id"], field=field, value=value,
                           immediate_source="mysql.samples",
                           historical_sources="; ".join(sorted(set(names))),
                           evidence_status="historical_audit_match" if names else "untraced_in_surviving_audits",
                           source_urls="; ".join(dict.fromkeys(links)),
                           audit_files="; ".join(dict.fromkeys(audits)),
                           notes="; ".join(notes), quality_flags="; ".join(sorted(warnings))))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    print(f"Traced {len(output):,} filled cells: {args.output}")
    print("Fields:", dict(Counter(r["field"] for r in output)))
    print("Historical audit matches:", sum(bool(r["historical_sources"]) for r in output))
    print("Quality flags (cells, not participants):", dict(flags))
    print("A historical audit match is not an independent primary-record revalidation.")


if __name__ == "__main__":
    main()
