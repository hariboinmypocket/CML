"""Species-node lookup from a local NCBI taxonomy dump, read at most once.

A rank='species' row that carries a *strain* tax ID is a recurring defect here,
repaired by hand each time it is noticed: the `retaxid` action in
scripts/merge_duplicate_taxa.py's OVERRIDES does it for F. prausnitzii and
L. salivarius, scripts/fix_sphingobacterium_paucimobilis.py did it for
tax 1346330, and scripts/fix_concatenated_genus_taxa.py did it for 14 more.
MiMeDB is a strain-level source, so its records routinely carry strain IDs and
the collapse to species keeps whichever one `most_common` picked.

Resolving the strain to its species node at load time stops that at the source.
Only nodes.dmp is read -- rank and parent are all this needs, so names.dmp
(316MB) is never touched.

The dump is optional, as it is everywhere else in this project
(scripts/load_microbiomehd.py states the same contract): without it the loader
behaves exactly as before and says so once, rather than failing. That makes an
ingest's tax IDs depend on whether the dump is present, which is a real cost --
the alternative is a loader that cannot do what every repair script here does.
"""
from __future__ import annotations

import os
from pathlib import Path

from .config import require_local

_MISSING = object()
_cache: object = _MISSING


class SpeciesNodes:
    """rank and parent for every node, with a memoised species-ancestor walk."""

    def __init__(self, nodes_path: str | Path) -> None:
        self.rank_of: dict[int, str] = {}
        self.parent: dict[int, int] = {}
        with open(require_local(nodes_path), encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = line.split("\t|")
                taxid = int(parts[0])
                self.parent[taxid] = int(parts[1])
                self.rank_of[taxid] = parts[2].strip()
        self._memo: dict[int, int | None] = {}

    def species_node(self, taxid: int | None) -> int | None:
        """The species node at or above `taxid`, or None.

        None also when the tax ID is absent from the dump -- NCBI deletes nodes
        (1232666 and 1328388 are both gone), and a caller must leave those
        alone rather than guess.
        """
        if not taxid or taxid not in self.rank_of:
            return None
        if taxid in self._memo:
            return self._memo[taxid]
        chain: list[int] = []
        node, guard = taxid, 0
        found: int | None = None
        while node and node != 1 and guard < 60:
            if node in self._memo:
                found = self._memo[node]
                break
            chain.append(node)
            if self.rank_of.get(node) == "species":
                found = node
                break
            node = self.parent.get(node, 0)
            guard += 1
        for member in chain:
            self._memo[member] = found
        return found


def load(directory: str | None = None) -> SpeciesNodes | None:
    """The dump at GUTDB_TAXDUMP (or `directory`), or None when absent.

    Cached for the process, so a load that upserts thousands of taxa parses
    nodes.dmp once.
    """
    global _cache
    if _cache is not _MISSING and directory is None:
        return _cache  # type: ignore[return-value]
    base = directory or os.environ.get("GUTDB_TAXDUMP", "taxdump")
    nodes = Path(base) / "nodes.dmp"
    result = SpeciesNodes(nodes) if nodes.exists() else None
    if directory is None:
        _cache = result
    return result


def reset_cache() -> None:
    """Forget the cached dump; for tests that point GUTDB_TAXDUMP elsewhere."""
    global _cache
    _cache = _MISSING
