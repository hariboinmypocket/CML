from __future__ import annotations

from typing import Any, Iterator

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .transform import clean_str

MGNIFY_BASE = "https://www.ebi.ac.uk/metagenomics/api/v1"


class MGnifyClient:
    """Client for EBI MGnify's REST API: study samples and precomputed taxonomy."""

    def __init__(self, timeout: int = 90) -> None:
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "gutdb-mgnify-sync/1.0"})
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

    def _paginated(self, url: str, page_size: int = 250) -> Iterator[dict[str, Any]]:
        params = {"format": "json", "page_size": page_size}
        while url:
            response = self.session.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            payload = response.json()
            yield from payload.get("data", [])
            url = (payload.get("links") or {}).get("next")
            params = {}

    def study_samples(self, study_accession: str) -> list[dict[str, Any]]:
        url = f"{MGNIFY_BASE}/studies/{study_accession}/samples"
        return list(self._paginated(url))

    def study_analyses(self, study_accession: str) -> list[dict[str, Any]]:
        url = f"{MGNIFY_BASE}/studies/{study_accession}/analyses"
        return list(self._paginated(url))

    def analysis_taxonomy_ssu(self, analysis_accession: str) -> list[dict[str, Any]]:
        url = f"{MGNIFY_BASE}/analyses/{analysis_accession}/taxonomy/ssu"
        return list(self._paginated(url))


def sample_metadata_dict(sample: dict[str, Any]) -> dict[str, str]:
    attrs = sample.get("attributes", {})
    metadata = {clean_str(m.get("key")).casefold(): clean_str(m.get("value")) for m in attrs.get("sample-metadata", [])}
    return metadata


def best_analysis_per_sample(analyses: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Pick the highest pipeline-version analysis for each sample accession."""
    best: dict[str, dict[str, Any]] = {}
    for analysis in analyses:
        sample_rel = (analysis.get("relationships", {}).get("sample", {}) or {}).get("data") or {}
        sample_id = sample_rel.get("id")
        if not sample_id:
            continue
        version = clean_str(analysis.get("attributes", {}).get("pipeline-version"))
        current = best.get(sample_id)
        if current is None or version > clean_str(current.get("attributes", {}).get("pipeline-version")):
            best[sample_id] = analysis
    return best


def genus_abundances(taxonomy: list[dict[str, Any]]) -> list[tuple[str, int]]:
    """Return (name, count) for genus-rank nodes.

    16S SSU amplicon data barely resolves species (a handful of calls per sample
    against dozens of genus calls), so genus is the one rank that's both well
    populated and non-overlapping — mixing genus- and species-level counts in the
    same normalization would double-count reads and break the sum-to-1 invariant
    a relative_abundance column implies.
    """
    rows = []
    for node in taxonomy:
        attrs = node.get("attributes", {})
        if attrs.get("rank") != "genus":
            continue
        # MGnify separates binomials with underscores ("Abiotrophia_defectiva").
        # Left as-is they parse as a single-token genus, which both blocks NCBI
        # lookups and creates a second row for an organism already in the table.
        name = clean_str(attrs.get("name")).replace("_", " ")
        count = attrs.get("count")
        if name and isinstance(count, (int, float)) and count > 0:
            rows.append((name, int(count)))
    return rows
