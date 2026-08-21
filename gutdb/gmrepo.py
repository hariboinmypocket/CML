from __future__ import annotations

import json
from typing import Any, Mapping

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


class GMrepoClient:
    """Small client around GMrepo's documented REST endpoints."""

    def __init__(self, base_url: str, timeout: int = 90) -> None:
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "gutdb-populator/1.0"})
        retry = Retry(
            total=4,
            connect=4,
            read=4,
            status=4,
            backoff_factor=1.0,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset({"POST"}),
        )
        self.session.mount("https://", HTTPAdapter(max_retries=retry))

    def _post(self, endpoint: str, payload: dict[str, Any] | None = None) -> Any:
        # GMrepo runs Django with APPEND_SLASH. A GET without the trailing slash is
        # redirected, but a POST cannot be redirected while preserving its body, so
        # the server answers 500 instead. Every endpoint is normalized here; without
        # it, exactly the calls whose paths lack a slash fail, which looks like an
        # intermittent server outage rather than a client bug.
        endpoint = endpoint if endpoint.endswith("/") else f"{endpoint}/"
        response = self.session.post(
            f"{self.base_url}/api/{endpoint}",
            data=json.dumps(payload or {}) if payload else {},
            headers={"Content-Type": "application/json"},
            timeout=self.timeout,
        )
        response.raise_for_status()
        return response.json()

    def all_phenotypes(self) -> list[dict[str, Any]]:
        data = self._post("get_all_phenotypes")
        return data.get("phenotypes", data) if isinstance(data, dict) else data

    def curated_projects(self) -> list[dict[str, Any]]:
        data = self._post("getCuratedProjectsList", {})
        return data if isinstance(data, list) else data.get("projects", [])

    def project_runs(
        self, project_id: str, limit: int = 1000, skip: int = 0
    ) -> list[dict[str, Any]]:
        data = self._post(
            "getAllRunsByProjectIDAsync/",
            {"project_id": project_id, "limit": limit, "skip": skip},
        )
        return data if isinstance(data, list) else data.get("data", [])

    def project_abundance(self, project_id: str, mesh_id: str = "") -> dict[str, Any]:
        return self._post(
            "getMicrobeAbundancesByPhenotypeMeshIDAndProjectID",
            {"project_id": project_id, "mesh_id": mesh_id},
        )

    def full_taxonomic_profile(self, run_id: str) -> dict[str, Any]:
        return self._post("getFullTaxonomicProfileByRunID", {"run_id": run_id})

    @staticmethod
    def parse_taxonomic_profile(payload: Any) -> list[dict[str, Any]]:
        """Pull per-taxon abundance rows out of a run's taxonomic profile.

        GMrepo nests the profile differently across endpoints/versions, so this
        scans the payload for the first list-of-records that actually carries a
        taxon name and an abundance. It raises on an unrecognized shape instead
        of returning an empty list: a silent schema change must surface as a
        failed ingestion run, not as a sample that quietly loads zero taxa.
        """
        name_keys = ("scientific_name", "taxon_name", "organism", "name", "species")
        abundance_keys = ("relative_abundance", "abundance", "mean_abundance", "value")
        rank_keys_attr = ("taxon_rank_level", "taxonomic_rank", "rank", "level")
        taxid_keys = ("ncbi_taxon_id", "ncbi_tax_id", "taxon_id", "taxid")

        def pick(record: Mapping[str, Any], candidates: tuple[str, ...]) -> Any:
            folded = {str(k).strip().casefold(): v for k, v in record.items()}
            for key in candidates:
                if folded.get(key) not in (None, ""):
                    return folded[key]
            return None

        # A GMrepo profile response is an envelope: {"run": [...], "phenotypes": [...],
        # "genus": [...]} with a list per rank. Runs that were never profiled come
        # back as a valid envelope with the rank keys simply absent, which is real
        # "no data" rather than a schema change, so it must return empty instead of
        # raising - otherwise every unprofiled run looks like a hard failure.
        rank_keys = ("superkingdom", "phylum", "class", "order", "family", "genus", "species")
        # "phenotypes_exist" alone marks a run GMrepo simply does not hold (for
        # example one loaded here from another source). That is a valid empty
        # answer, not a schema change, and must not be raised as a failure.
        if isinstance(payload, Mapping) and (
            "run" in payload
            or "phenotypes_exist" in payload
            or any(k in payload for k in rank_keys)
        ):
            envelope: list[dict[str, Any]] = []
            for key in rank_keys:
                value = payload.get(key)
                if isinstance(value, list):
                    for record in value:
                        if not isinstance(record, Mapping):
                            continue
                        name = pick(record, name_keys)
                        abundance = pick(record, abundance_keys)
                        if name is None or abundance is None:
                            continue
                        envelope.append({
                            "scientific_name": str(name).strip(),
                            "relative_abundance": abundance,
                            "taxonomic_rank": pick(record, rank_keys_attr) or key,
                            "ncbi_tax_id": pick(record, taxid_keys),
                        })
            return envelope

        candidate_lists: list[list[Any]] = []
        if isinstance(payload, list):
            candidate_lists.append(payload)
        elif isinstance(payload, dict):
            for value in payload.values():
                if isinstance(value, list) and value:
                    candidate_lists.append(value)

        rows: list[dict[str, Any]] = []
        for candidate in candidate_lists:
            for record in candidate:
                if not isinstance(record, Mapping):
                    continue
                name = pick(record, name_keys)
                abundance = pick(record, abundance_keys)
                if name is None or abundance is None:
                    continue
                rows.append({
                    "scientific_name": str(name).strip(),
                    "relative_abundance": abundance,
                    "taxonomic_rank": pick(record, rank_keys_attr),
                    "ncbi_tax_id": pick(record, taxid_keys),
                })

        if not rows:
            observed: set[str] = set()
            for candidate in candidate_lists:
                for record in candidate[:5]:
                    if isinstance(record, Mapping):
                        observed.update(str(k) for k in record)
            raise ValueError(
                "Unrecognized GMrepo taxonomic-profile shape; no record carried both a "
                f"taxon name and an abundance. Top-level type={type(payload).__name__}, "
                f"top-level keys={sorted(payload)[:12] if isinstance(payload, dict) else 'n/a'}, "
                f"record keys seen={sorted(observed)[:20]}"
            )
        return rows

    def phenotype_comparison(self, mesh_id1: str, mesh_id2: str) -> dict[str, Any]:
        """Return current GMrepo marker rows across curated projects."""
        return self._post(
            "getPhenotypeComparisonsDetails/",
            {"mesh_id1": mesh_id1, "mesh_id2": mesh_id2},
        )

    def all_phenotype_comparisons(self) -> dict[str, Any]:
        return self._post("get_all_phenotype_comparisons/", {})

    def curated_project_comparison(
        self, project_id: str, phenotype1: str, phenotype2: str
    ) -> dict[str, Any]:
        return self._post(
            "getCuratedProjectDetailsOnPhenotypeComparison/",
            {
                "ncbi_project_id": project_id,
                "phenotype1": phenotype1,
                "phenotype2": phenotype2,
            },
        )
