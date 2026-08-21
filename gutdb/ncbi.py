from __future__ import annotations

import re
import time
from typing import Iterator
from xml.etree import ElementTree as ET

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .transform import clean_str

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Ordered by how specific/trustworthy the tag is; first match wins.
SEX_TAGS = ("host_sex", "sex", "host sex", "gender", "host_gender")
AGE_TAGS = ("host_age", "age", "host age", "patient_age", "host_age_years")
BMI_TAGS = ("host_body_mass_index", "bmi", "host_bmi")

NA_TEXT = {
    "not collected", "not applicable", "missing", "n/a", "na", "unknown",
    "not provided", "restricted access", "-", "not determined", "not recorded",
}

# NCBI's harmonized host_age package defaults to years when no unit is given.
_UNIT_TO_YEARS = {
    "d": 1 / 365.25, "day": 1 / 365.25, "days": 1 / 365.25,
    "wk": 7 / 365.25, "wks": 7 / 365.25, "week": 7 / 365.25, "weeks": 7 / 365.25,
    "mo": 1 / 12, "mos": 1 / 12, "month": 1 / 12, "months": 1 / 12,
    "y": 1.0, "yr": 1.0, "yrs": 1.0, "year": 1.0, "years": 1.0,
}
_NUMBER_RE = re.compile(r"(\d+(?:\.\d+)?)")
_UNIT_RE = re.compile(r"\b(days?|d|wks?|weeks?|mos?|months?|y|yrs?|years?)\b", re.IGNORECASE)
_MAX_PLAUSIBLE_AGE_YEARS = 122.0


def parse_host_age_years(raw_value: str) -> tuple[float | None, str]:
    """Convert a free-text BioSample host_age value to years.

    Returns (years, note). `years` is None when the value can't be trusted;
    `note` records why, for the audit trail.
    """
    text = clean_str(raw_value)
    if not text:
        return None, "empty"
    if text.casefold() in NA_TEXT:
        return None, "not_collected"
    match = _NUMBER_RE.search(text)
    if not match:
        return None, "unparseable"
    number = float(match.group(1))
    unit_match = _UNIT_RE.search(text)
    unit = unit_match.group(1).casefold() if unit_match else None
    factor = _UNIT_TO_YEARS.get(unit, 1.0)
    years = number * factor
    if years < 0 or years > _MAX_PLAUSIBLE_AGE_YEARS:
        return None, f"out_of_range:{years:.1f}y"
    note = "censored_lower_bound" if "+" in text else (f"unit={unit}" if unit else "assumed_years")
    return round(years, 2), note


def parse_host_sex(raw_value: str) -> str | None:
    text = clean_str(raw_value).casefold()
    if text in {"male", "m"}:
        return "male"
    if text in {"female", "f"}:
        return "female"
    return None


def parse_bmi(raw_value: str) -> float | None:
    text = clean_str(raw_value)
    if not text or text.casefold() in NA_TEXT:
        return None
    match = _NUMBER_RE.search(text)
    if not match:
        return None
    value = float(match.group(1))
    if value < 10 or value > 80:
        return None
    return round(value, 1)


def batched(items: list[str], size: int) -> Iterator[list[str]]:
    for i in range(0, len(items), size):
        yield items[i : i + size]


class SRAClient:
    """Fetches BioSample attributes embedded in SRA experiment XML by run accession."""

    def __init__(
        self, timeout: int = 90, email: str | None = None, api_key: str | None = None
    ) -> None:
        self.timeout = timeout
        self.email = email
        self.api_key = api_key
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "gutdb-demographics-enrichment/1.0"})
        retry = Retry(
            total=5,
            connect=5,
            read=5,
            status=5,
            backoff_factor=1.5,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset({"GET"}),
        )
        self.session.mount("https://", HTTPAdapter(max_retries=retry))

    def fetch_run_attributes(
        self, run_accessions: list[str], attempts: int = 4
    ) -> dict[str, dict[str, str]]:
        """Return {run_accession: {tag: value}} for each run's BioSample attributes.

        Retries the whole request (not just the connect/status phase the
        urllib3 Retry adapter covers) because NCBI occasionally truncates a
        large chunked response mid-body, which only surfaces while reading
        `response.text` after the adapter has already handed back a 200.
        """
        params = {
            "db": "sra",
            "id": ",".join(run_accessions),
            "rettype": "full",
            "retmode": "xml",
        }
        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key
        last_error: Exception | None = None
        for attempt in range(1, attempts + 1):
            try:
                response = self.session.get(
                    f"{EUTILS_BASE}/efetch.fcgi", params=params, timeout=self.timeout
                )
                response.raise_for_status()
                return _parse_experiment_packages(response.text, wanted=set(run_accessions))
            except requests.exceptions.RequestException as error:
                last_error = error
                if attempt < attempts:
                    time.sleep(2**attempt)
        assert last_error is not None
        raise last_error


# NCBI renamed the top-level rank from "superkingdom" to "domain"; accept either
# so the value survives whichever vocabulary a given record was built with.
# NCBI renamed the top-level rank from "superkingdom" to "domain"; accept either
# so the value survives whichever vocabulary a given record was built with.
LINEAGE_RANKS = ("superkingdom", "domain", "phylum", "class", "order", "family")

# Root taxids used to constrain a name lookup to the right branch of the tree.
DOMAIN_SUBTREE = {
    "bacteria": 2,
    "archaea": 2157,
    "eukaryota": 2759,
    "viruses": 10239,
}

# Phyla that cannot belong to a prokaryote. A row carrying one of these while
# superkingdom says Bacteria/Archaea is internally contradictory and is the
# signature of a homonym mis-resolution.
NON_PROKARYOTIC_PHYLA = frozenset({
    "chordata", "arthropoda", "mollusca", "annelida", "nematoda", "echinodermata",
    "streptophyta", "chlorophyta", "cnidaria", "platyhelminthes", "porifera",
    "rotifera", "ascomycota", "basidiomycota", "apicomplexa", "ciliophora",
    "euglenozoa", "amoebozoa", "microsporidia", "bryozoa", "tardigrada",
    "nemertea", "brachiopoda", "mucoromycota", "zoopagomycota", "oomycota",
})


class TaxonomyClient(SRAClient):
    """NCBI Taxonomy lookups: taxid -> lineage, and verified name -> taxid."""

    def fetch_lineages(self, taxids: list[str]) -> dict[str, dict[str, str]]:
        """Return {taxid: {rank: name, '_scientific_name': ..., '_rank': ...}}."""
        params = {"db": "taxonomy", "id": ",".join(str(t) for t in taxids), "retmode": "xml"}
        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key
        response = self.session.get(
            f"{EUTILS_BASE}/efetch.fcgi", params=params, timeout=self.timeout
        )
        response.raise_for_status()
        try:
            root = ET.fromstring(response.text)
        except ET.ParseError:
            return {}
        results: dict[str, dict[str, str]] = {}
        for taxon in root.findall("Taxon"):
            taxid = clean_str(taxon.findtext("TaxId"))
            if not taxid:
                continue
            entry = {
                "_scientific_name": clean_str(taxon.findtext("ScientificName")),
                "_rank": clean_str(taxon.findtext("Rank")),
            }
            for node in taxon.findall("./LineageEx/Taxon"):
                rank = clean_str(node.findtext("Rank")).casefold()
                if rank in LINEAGE_RANKS:
                    entry[rank] = clean_str(node.findtext("ScientificName"))
            results[taxid] = entry
        return results

    def resolve_name(self, name: str, domain: str = "") -> str:
        """Return a taxid only when NCBI's record genuinely matches this organism.

        Two distinct failure modes are guarded here, both of which produce
        confident-looking wrong lineages:

        1. Fuzzy matching. The composite 16S label "Escherichia-Shigella"
           returns Shigella's taxid. Caught by comparing ScientificName.
        2. Cross-kingdom homonyms. "Proteus" is both a bacterial genus and a
           salamander genus; "Rothia" is both a bacterium and a plant. A name
           check passes happily and grafts an animal phylum onto a bacterium.
           Caught by constraining the search to the known domain subtree and
           re-checking the returned lineage's domain.
        """
        query = clean_str(name)
        if not query:
            return ""
        subtree = DOMAIN_SUBTREE.get(clean_str(domain).casefold())
        term = f"{query} AND txid{subtree}[Subtree]" if subtree else query
        params = {"db": "taxonomy", "term": term, "retmode": "json"}
        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key
        response = self.session.get(
            f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=self.timeout
        )
        response.raise_for_status()
        try:
            idlist = response.json().get("esearchresult", {}).get("idlist", [])
        except ValueError:
            return ""
        if not idlist:
            return ""
        candidate = idlist[0]
        record = self.fetch_lineages([candidate]).get(candidate)
        if not record:
            return ""
        if record.get("_scientific_name", "").casefold() != query.casefold():
            return ""
        if domain:
            found_domain = record.get("domain") or record.get("superkingdom") or ""
            if found_domain and found_domain.casefold() != clean_str(domain).casefold():
                return ""
        return candidate


def _parse_experiment_packages(xml_text: str, wanted: set[str]) -> dict[str, dict[str, str]]:
    results: dict[str, dict[str, str]] = {}
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return results
    for package in root.findall("EXPERIMENT_PACKAGE"):
        sample = package.find(".//SAMPLE")
        attrs: dict[str, str] = {}
        if sample is not None:
            for attribute in sample.findall(".//SAMPLE_ATTRIBUTE"):
                tag = clean_str(attribute.findtext("TAG")).casefold()
                value = clean_str(attribute.findtext("VALUE"))
                if tag and value:
                    attrs[tag] = value
        for run in package.findall(".//RUN"):
            run_accession = run.get("accession")
            if run_accession and run_accession in wanted:
                results[run_accession] = attrs
    return results
