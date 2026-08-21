from __future__ import annotations

import csv
import io
from typing import Any, Iterator

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .transform import clean_str, nullable_float

CMD_API = "https://api.github.com/repos/waldronlab/curatedMetagenomicDataCuration/contents/inst/curated"
CMD_RAW = "https://raw.githubusercontent.com/waldronlab/curatedMetagenomicDataCuration/master/inst/curated"

# cMD records age in the unit named by its own age_unit column. Anything other
# than years must be converted, not assumed: an unconverted infant cohort is how
# this project previously ended up with 122-year-old neonates.
_AGE_UNIT_TO_YEARS = {
    "year": 1.0, "years": 1.0, "y": 1.0,
    "month": 1 / 12, "months": 1 / 12,
    "week": 7 / 365.25, "weeks": 7 / 365.25,
    "day": 1 / 365.25, "days": 1 / 365.25,
}
_MAX_PLAUSIBLE_AGE_YEARS = 122.0


class CMDCurationClient:
    """Reads curatedMetagenomicData's curation repo over plain HTTP.

    cMD demographics are transcribed from the papers by human curators, so they
    can carry age/sex for samples where BioSample and ENA hold nothing. The
    taxonomic profiles live in Bioconductor's ExperimentHub and are not reachable
    here; only the metadata is.
    """

    def __init__(self, timeout: int = 90) -> None:
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "gutdb-cmd-curation/1.0"})
        self.session.mount("https://", HTTPAdapter(max_retries=Retry(
            total=4, backoff_factor=1.5, status_forcelist=(429, 500, 502, 503, 504))))

    def studies(self) -> list[str]:
        response = self.session.get(CMD_API, timeout=self.timeout)
        response.raise_for_status()
        return [item["name"] for item in response.json() if item.get("type") == "dir"]

    def study_samples(self, study: str) -> list[dict[str, str]]:
        for suffix in ("_sample.tsv", "_metadata.tsv"):
            try:
                response = self.session.get(
                    f"{CMD_RAW}/{study}/{study}{suffix}", timeout=self.timeout
                )
            except requests.RequestException:
                continue
            if response.status_code == 200 and response.text.strip():
                return list(csv.DictReader(io.StringIO(response.text), delimiter="\t"))
        return []


def run_accessions(row: dict[str, str]) -> Iterator[str]:
    """A cMD row may list several runs for one biological sample."""
    raw = clean_str(row.get("ncbi_accession"))
    for token in raw.replace(";", ",").split(","):
        token = token.strip()
        if token:
            yield token


def parse_sex(row: dict[str, str]) -> str | None:
    """cMD names the column `sex`; `gender` is only a fallback."""
    value = clean_str(row.get("sex")) or clean_str(row.get("gender"))
    folded = value.casefold()
    if folded in {"male", "m"}:
        return "male"
    if folded in {"female", "f"}:
        return "female"
    return None


def parse_age_years(row: dict[str, str]) -> tuple[float | None, str]:
    """Convert cMD age to years using its own age_unit, or reject it."""
    value = nullable_float(row.get("age"))
    if value is None:
        return None, "empty"
    unit = clean_str(row.get("age_unit")).casefold() or "year"
    factor = _AGE_UNIT_TO_YEARS.get(unit)
    if factor is None:
        return None, f"unknown_unit:{unit}"
    years = value * factor
    if years < 0 or years > _MAX_PLAUSIBLE_AGE_YEARS:
        return None, f"out_of_range:{years:.1f}"
    return round(years, 2), f"unit={unit}"


def parse_bmi(row: dict[str, str]) -> float | None:
    value = nullable_float(row.get("bmi"))
    if value is None or value < 10 or value > 80:
        return None
    return round(value, 1)
