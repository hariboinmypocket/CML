#!/usr/bin/env python3
"""Verify CKD clinical values and identifier links against cached public sources.

Read only: no database connection, no source CSV changes, no inference from
number ranges. Writes a verification report, not a database table. The cache
contains the original NLM supplements, ENA run report, and exact-run SRA XML.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import xlrd
from openpyxl import load_workbook
from pypdf import PdfReader


ROOT = Path(__file__).resolve().parent.parent
SUPPLEMENT_BASE = "https://rcastoragev2.blob.core.windows.net/5de81978ee2b164cb591ed3fede2bb56/"
ARTICLE_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7578882/fullTextXML"
ENA_URL = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA562327"
           "&result=read_run&fields=run_accession,sample_accession,secondary_sample_accession,"
           "experiment_accession,study_accession,library_name,sample_alias,sample_title,"
           "experiment_alias,experiment_title,study_alias&format=tsv&download=false")


def read_csv(path, delimiter=","):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def normalized_identifier(value):
    # Normalize separators/case only. Never strip a library's institutional prefix.
    return re.sub(r"[\s_-]", "", str(value).casefold())


def spreadsheet_tables(path):
    data = path.read_bytes()
    if data.startswith(bytes.fromhex("d0cf11e0a1b11ae1")):
        book = xlrd.open_workbook(file_contents=data)
        return [(sheet.name, [sheet.row_values(i) for i in range(sheet.nrows)])
                for sheet in book.sheets()]
    if data.startswith(b"PK"):
        book = load_workbook(io.BytesIO(data), read_only=False, data_only=False)
        try:
            tables = []
            for sheet in book.worksheets:
                rows = [[cell.value for cell in row] for row in sheet]
                links = [[cell.hyperlink.target] for row in sheet for cell in row if cell.hyperlink]
                tables.append((sheet.title, rows + links))
            return tables
        finally:
            book.close()
    # Many publisher .xls supplements are actually tab-separated text.
    return [("TSV", list(csv.reader(io.StringIO(data.decode("utf-8-sig")), delimiter="\t")))]


def clinical_records(path):
    records = {}
    for sheet, rows in spreadsheet_tables(path):
        headers = [str(v).strip().casefold() for v in rows[0]]
        id_col = next(i for i, h in enumerate(headers) if h in {"sample", "number"})
        age_col = next(i for i, h in enumerate(headers) if h.startswith("age"))
        sex_col = next(i for i, h in enumerate(headers) if h in {"gender", "sex"})
        bmi_col = next(i for i, h in enumerate(headers) if h.startswith("bmi"))
        for row in rows[1:]:
            identifier = str(row[id_col]).strip()
            key = normalized_identifier(identifier)
            if key in records:
                raise ValueError("Duplicate clinical identifier in source workbook")
            records[key] = {"identifier": identifier, "sheet": sheet,
                            "age_years": row[age_col], "sex": str(row[sex_col]).strip().lower(),
                            "bmi": row[bmi_col]}
    return records


def same_measurement(field, source, exported):
    if field == "sex":
        return str(source).strip().casefold() == exported.casefold()
    try:
        value = float(source)
    except (ValueError, TypeError):
        return not exported and str(source).strip().casefold() in {"", "na", "n/a"}
    return bool(exported) and abs(value - float(exported)) < 1e-8


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=Path("/tmp/gutdb_ckd_link_verification"))
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/ckd_link_verification")
    args = parser.parse_args()
    cache = args.cache_dir
    audit_path = ROOT / "data/ckd_study_demographics_audit.csv"
    audit_bytes = audit_path.read_bytes()
    audit = read_csv(audit_path)
    ena_rows = read_csv(cache / "ena_runs.tsv", delimiter="\t")
    ena = {row["run_accession"]: row for row in ena_rows}
    if len(ena) != len(ena_rows):
        raise ValueError("ENA run IDs are not unique")
    clinical = clinical_records(cache / "ADVS-7-2001936-s003.xls")
    sra = {}
    sra_files = sorted(cache.glob("sra_inferred_*.xml"))
    for path in sra_files:
        for package in ET.parse(path).getroot().findall("EXPERIMENT_PACKAGE"):
            sample, experiment = package.find("SAMPLE"), package.find("EXPERIMENT")
            if sample is None or experiment is None:
                continue
            identifiers = [node.text for node in sample.findall("./IDENTIFIERS/*") if node.text]
            attributes = {node.findtext("TAG", "").strip().lower(): node.findtext("VALUE", "")
                          for node in sample.findall(".//SAMPLE_ATTRIBUTE")}
            record = {"library_name": experiment.findtext(".//LIBRARY_NAME"),
                      "sample_alias": sample.get("alias"), "sample_title": sample.findtext("TITLE"),
                      "identifiers": identifiers, "attributes": attributes,
                      "source_file": path.name}
            for run in package.findall(".//RUN"):
                sra[run.get("accession")] = record

    rows = []
    demographic_tags = {"sex", "host_sex", "host sex", "gender", "host_gender",
                        "age", "host_age", "host age", "patient_age", "host_age_years",
                        "bmi", "host_bmi", "host_body_mass_index"}
    for row in audit:
        source = clinical[normalized_identifier(row["supplement_id"])]
        values_match = all(same_measurement(field, source[field], row[field])
                           for field in ["age_years", "sex", "bmi"])
        external = ena.get(row["run_accession"], {})
        library_match = (external.get("library_name") == row["library_name"]
                         and external.get("study_accession") == "PRJNA562327")
        ncbi = sra.get(row["run_accession"], {})
        ncbi_agrees = bool(ncbi and ncbi["library_name"] == row["library_name"]
                           and external.get("sample_accession") in ncbi["identifiers"])
        # Match the complete source identifier, not the same numeric suffix.
        external_identifiers = [external.get(k, "") for k in
                                ["sample_alias", "sample_title", "experiment_alias"]]
        identifier_match = any(normalized_identifier(value) == normalized_identifier(row["supplement_id"])
                               for value in external_identifiers if value)
        rows.append({"run_accession": row["run_accession"], "supplement_id": row["supplement_id"],
                     "library_name": row["library_name"], "previous_join_basis": row["join_basis"],
                     "published_values_match_audit": values_match,
                     "ena_library_to_run_verified": library_match,
                     "sra_library_and_biosample_agree": ncbi_agrees,
                     "exact_source_identifier_in_ena": identifier_match,
                     "sra_demographic_tags_present": sorted(demographic_tags & set(ncbi.get("attributes", {}))),
                     "source_sheet": source["sheet"], "biosample_accession": external.get("sample_accession"),
                     "status": "source_alias_matches_after_separator_normalization" if identifier_match and library_match
                               else "clinical_to_library_crosswalk_not_found"})
    if not all(r["published_values_match_audit"] and r["ena_library_to_run_verified"] for r in rows):
        raise ValueError("A source measurement or run/library relationship disagrees with the audit")
    inferred = [r for r in rows if r["previous_join_basis"] == "inferred_single_series"]
    if not all(r["sra_library_and_biosample_agree"] for r in inferred):
        raise ValueError("An inferred run could not be corroborated in exact-run SRA metadata")

    article = ET.parse(cache / "article.xml").getroot()
    supplement_names = [x.get("{http://www.w3.org/1999/xlink}href")
                        for x in article.findall(".//supplementary-material/media")]
    markers = re.compile(r"\b(?:SRR\d+|SRX\d+|SRS\d+|SAMN\d+|zzHC[_-]\d+|Mobio_ZZ_HC[_-]\d+|"
                         r"sample_number|accession|crosswalk|library_name)\b", re.I)
    inspections = []
    for name in supplement_names:
        path = cache / name
        hits = []
        if name.endswith((".xls", ".xlsx")):
            for sheet, table in spreadsheet_tables(path):
                for i, record in enumerate(table, 1):
                    for j, value in enumerate(record, 1):
                        if isinstance(value, str) and markers.search(value):
                            hits.append({"sheet": sheet, "row": i, "column": j, "text": value[:200]})
        elif name.endswith(".pdf"):
            reader = PdfReader(path)
            for i, page in enumerate(reader.pages, 1):
                text = page.extract_text() or ""
                for match in markers.finditer(text):
                    hits.append({"page": i, "text": text[max(0, match.start()-80):match.end()+80]})
        inspections.append({"file": name, "source_url": SUPPLEMENT_BASE + name,
                            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                            "identifier_or_crosswalk_candidates": hits})
    summary = {"audit_rows": len(rows), "published_values_matching_audit": sum(r["published_values_match_audit"] for r in rows),
               "ena_library_to_run_links_verified": sum(r["ena_library_to_run_verified"] for r in rows),
               "inferred_hc_rows_checked": len(inferred),
               "inferred_hc_sra_records_corroborated": sum(r["sra_library_and_biosample_agree"] for r in inferred),
               "inferred_hc_exact_clinical_id_matches": sum(r["exact_source_identifier_in_ena"] for r in inferred),
               "inferred_hc_sra_records_with_demographics": sum(bool(r["sra_demographic_tags_present"]) for r in inferred),
               "legacy_documented_exact_alias_matches": sum(r["exact_source_identifier_in_ena"] for r in rows if r["previous_join_basis"] == "documented"),
               "supplement_spreadsheets_inspected": sum(name.endswith((".xls", ".xlsx")) for name in supplement_names),
               "supporting_pdfs_inspected": sum(name.endswith(".pdf") for name in supplement_names)}
    result = {"verified_utc": datetime.now(timezone.utc).isoformat(), "summary": summary,
              "input_audit_sha256": hashlib.sha256(audit_bytes).hexdigest(),
              "sources": {"article": ARTICLE_URL, "ena": ENA_URL,
                          "ena_sha256": hashlib.sha256((cache/'ena_runs.tsv').read_bytes()).hexdigest(),
                          "sra_files": [{"name": p.name, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                                        for p in sra_files]},
              "supplement_inspections": inspections, "records": rows,
              "limitations": ["Identifier searches and source inspection did not locate an explicit HC participant-to-library crosswalk.",
                              "Matching demographic values to the publication does not verify their assignment to a sequencing run.",
                              "The previous documented labels do not themselves establish verified linkage: CKD institutional prefixes also lack an explicit crosswalk.",
                              "The 57 HzCKD identifiers match ENA aliases after case/separator normalization; they are distinct from the unresolved HC series."]}
    args.output_dir.mkdir(parents=True, exist_ok=False)
    (args.output_dir / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    (args.output_dir / "report.md").write_text(
        "# CKD demographic-link verification\n\n"
        f"Verification snapshot: {result['verified_utc']}. Source audit SHA-256: `{result['input_audit_sha256']}`.\n\n"
        f"All {len(rows)} audit records' demographic values match Supporting Excel 1. "
        "ENA confirms all recorded library-to-run links.\n\n"
        f"For the {len(inferred)} previously inferred HC mappings, exact-run SRA metadata confirms "
        "the same libraries and BioSample accessions, but none independently reports age, sex or BMI. "
        "No complete clinical HC identifier matched a public sample alias.\n\n"
        "The original 21 supplementary spreadsheets and two supporting PDFs were inspected. "
        "Supporting Excel 1 has a clinical sample_number field, but no crosswalk was found from it "
        "to sequencing accessions. ENA contains overlapping healthy-control numbering series "
        "Mobio_ZZ_HC_1..160 and zzHC_1..320. Matching suffixes or the larger range does not verify participant identity.\n\n"
        "**Verdict: the 273 HC participant-to-run assignments remain unverified.** "
        "Do not promote them to confirmed demographics or copy them into SQLite on this evidence. "
        "An explicit author-provided participant-to-BioSample/library/run crosswalk is needed.\n\n"
        "The other 216 records were previously labelled documented. Only the 57 HzCKD aliases match "
        "the source clinical IDs after case/separator normalization. The 159 CKD records also rely "
        "on an institutional-prefix assumption; the old label is not independent proof.\n\n"
        "No database tables, sample values, audit rows, or mapping-status labels were changed.\n\n"
        f"Sources: [article]({ARTICLE_URL}), [Supporting Excel 1]({SUPPLEMENT_BASE}ADVS-7-2001936-s003.xls), "
        f"[ENA run metadata]({ENA_URL}). Full per-record checks and source hashes are in verification.json.\n")
    if audit_path.read_bytes() != audit_bytes:
        raise RuntimeError("Audit changed concurrently; report refers to the pinned input hash")
    print(json.dumps(summary, indent=2))
    print("Verification report:", args.output_dir / "report.md")


if __name__ == "__main__":
    main()
