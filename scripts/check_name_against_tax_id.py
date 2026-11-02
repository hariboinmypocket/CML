#!/usr/bin/env python3
"""Find rows whose stored name and `ncbi_tax_id` describe different organisms.

Every repair in this project resolves identity through `ncbi_tax_id`, so a row
whose name and tax ID disagree is worse than a row with no tax ID: the ID wins
every automatic decision, silently. `Vibrio fischeri` carried 29494 -- which is
*Vibrio furnissii* -- in both stores, and nothing noticed until a duplicate
collision surfaced it from the side.

The checks that exist cannot see this class.
`scripts/repoint_strain_taxids.py` only examines rows whose tax ID *is* a strain
node, so a row holding the wrong species' ID outright is invisible to it.
`taxa_sharing_a_tax_id_within_one_rank` only fires on the subset where two rows
collide on one ID; a wrong ID that nothing else holds sits undetected.

WHY THE COMPARISON IS AGAINST SYNONYMS, NOT THE SCIENTIFIC NAME

Comparing a stored name with the dump's scientific name alone would flag the
correct rows in this database by the hundred. NCBI's scientific name for 668 is
*Aliivibrio fischeri* and this table says `Vibrio fischeri`, which the dump lists
as a synonym of that very node -- the row is right and the name is simply the
pre-2007 one. The same holds for every stale-but-real name the project
deliberately keeps (Firmicutes over Bacillota, Streptomycetales over
Kitasatosporales). So a row agrees with its tax ID when its name matches ANY name
the dump carries for that node, in any name class.

Two further tolerances, both already this project's rules rather than new ones:

  * A differing genus with an agreeing epithet is NOT flagged. That is the
    name-variant drift `merge_duplicate_taxa.py` canonicalizes -- `Prevotella
    ruminicola` against NCBI's *Xylanibacter ruminicola* is one organism under an
    older genus, not a wrong tax ID.
  * `same_epithet` supplies Latin gender agreement, so `glabrata`/`glabratus`
    and `spiroforme`/`spiroformis` agree.

What is left when both tolerances are applied is a row whose epithet matches
nothing the dump knows for its tax ID -- which is the Vibrio fischeri shape, and
is a wrong ID rather than an old name.

Reports only. A wrong tax ID can be wrong in either direction -- the name may be
right and the ID a typo, or the ID right and the name mislabelled -- and which
one to correct depends on the row's evidence and source. Deciding that per row is
not something to automate, and the Vibrio pair needed a human reading the dump to
settle.

Needs the taxdump; pass GUTDB_TAXDUMP or unzip taxdmp.zip into ./taxdump.
"""
from __future__ import annotations

import csv
import os
import re
import sqlite3
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.transform import clean_str, parse_scientific_name, same_epithet

TAXDUMP = os.environ.get("GUTDB_TAXDUMP", "taxdump")
SQLITE_PATH = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")
AUDIT = Path("data/name_against_tax_id_audit.csv")

# Name classes worth trusting as "this node has been called that". `authority`
# carries the publication string ("Vibrio fischeri (Beijerinck 1889) Lehmann and
# Neumann 1896"), which still begins with the name, and parse_scientific_name
# takes the first two tokens, so it contributes correctly.
NAME_CLASSES = ("scientific name", "synonym", "equivalent name", "includes",
                "genbank synonym", "authority", "type material")


def normalize(text: str) -> str:
    """Formatting variants of one name must not read as a disagreement.

    Three kinds of decoration, each of which produced false positives on the
    first run: underscores where the table uses spaces (`sp. 1_1_14`), brackets
    around a misapplied name (`[Clostridium] leptum`), and the quoting and
    parentheses the dump puts round names that are not validly published
    (`"candidatus fonsibacter ubiquis" henson et al. 2018`,
    `Anaerotruncus massiliensis (ex togo et al. 2026)`).
    """
    text = clean_str(text).replace("_", " ")
    text = re.sub(r"\(.*?\)", " ", text)
    text = text.replace("[", "").replace("]", "").replace('"', "").replace("'", "")
    return re.sub(r"\s+", " ", text).strip()


def load_names() -> dict[int, set[tuple[str, str]]]:
    """taxid -> every (genus, epithet) pair the dump carries for it."""
    path = f"{TAXDUMP}/names.dmp"
    if not os.path.exists(path):
        raise SystemExit(
            f"taxdump not found at {TAXDUMP}/. Unzip taxdmp.zip there or set GUTDB_TAXDUMP.")
    pairs: dict[int, set[tuple[str, str]]] = {}
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split("\t|")
            if len(parts) < 4 or parts[3].strip() not in NAME_CLASSES:
                continue
            genus, epithet = parse_scientific_name(normalize(parts[1]))
            if not genus:
                continue
            pairs.setdefault(int(parts[0]), set()).add(
                (genus.casefold(), epithet.casefold()))
    return pairs


def head(epithet: str) -> str:
    """The epithet proper, without a trailing strain designation or authority.

    The dump's name for a strain node puts the designation after the epithet
    ("Clostridium difficile CD196"), and a name that is not validly published
    carries its authority ("Candidatus Fonsibacter ubiquis" Henson et al. 2018),
    so only the leading word is the epithet.
    """
    return epithet.split()[0] if epithet else ""


def comparable(epithet: str) -> bool:
    """Whether a token can be read as an epithet at all.

    Catalogue numbers reach here as names: the dump lists "DSM 41672" and
    "As 4.1618" for taxid 285555, which parse to an epithet of "41672". Treating
    those as epithets makes every comparison against such a node a mismatch.
    """
    return len(re.findall(r"[A-Za-z]", head(epithet))) >= 3


def classify(name: str, genus_key: str, species: str, rank: str,
             known: set[tuple[str, str]] | None) -> str:
    """Whether a row's name and its tax ID describe the same organism.

    `known` is every (genus, epithet) pair the dump carries for the row's tax
    ID, in any name class -- which is the point: NCBI's scientific name for 668
    is *Aliivibrio fischeri* while this table says `Vibrio fischeri`, a synonym
    of that node, and the row is right.
    """
    if not known:
        return "tax_id_not_in_dump"
    stored_genus, stored_epithet = parse_scientific_name(normalize(name))
    stored_genus = stored_genus.casefold() or clean_str(genus_key).casefold()
    stored_epithet = stored_epithet.casefold()

    if (stored_genus, stored_epithet) in known:
        return "agrees_exactly"
    if rank != "species" or not stored_epithet:
        if any(stored_genus == genus for genus, _ in known):
            return "agrees_genus_only"
    # The dump's name for a STRAIN node carries the designation after the
    # epithet -- 645462 is "Clostridium difficile CD196" -- so a row named for
    # the species agrees with it. Compare the leading word of each epithet.
    if any(same_epithet(head(stored_epithet), head(epithet))
           for _, epithet in known if comparable(epithet)):
        return "agrees_on_the_epithet"
    if not any(comparable(epithet) for _, epithet in known):
        # Nodes whose every name is type material ("DSM 41672") offer no
        # binomial to compare, so this check has nothing to say about them.
        return "no_binomial_at_this_node"
    return "DISAGREES"


def run_store(store: str, cursor, names, audit: list[dict]) -> Counter:
    cursor.execute(
        "SELECT id, scientific_name, genus_key, species, ncbi_tax_id, taxonomic_rank "
        "FROM taxa WHERE ncbi_tax_id IS NOT NULL")
    stats: Counter = Counter()
    for taxon_id, name, genus_key, species, taxid, rank in cursor.fetchall():
        known = names.get(int(taxid))
        verdict = classify(name, genus_key, species, rank, known)
        stats[verdict] += 1
        if verdict != "DISAGREES":
            continue
        example = sorted(f"{g.capitalize()} {e}".strip() for g, e in known)[:3]
        audit.append({
            "store": store, "taxon_id": taxon_id, "scientific_name": name,
            "taxonomic_rank": rank, "ncbi_tax_id": taxid,
            "dump_calls_it": "; ".join(example),
        })
        print(f"    {taxon_id:>8} {name[:40]:<40} tax {taxid:<9} "
              f"dump: {', '.join(example)}")
    return stats


def main() -> int:
    print(f"loading taxdump from {TAXDUMP}/ ...")
    names = load_names()
    print(f"  nodes with at least one parseable name: {len(names):,}")

    audit: list[dict] = []
    totals = {}

    print("\nmysql")
    conn = connect(Settings.from_env(".env"))
    cursor = conn.cursor()
    totals["mysql"] = run_store("mysql", cursor, names, audit)
    cursor.close()
    conn.close()

    print("\nsqlite")
    if Path(SQLITE_PATH).exists():
        lite = sqlite3.connect(SQLITE_PATH)
        cursor = lite.cursor()
        totals["sqlite"] = run_store("sqlite", cursor, names, audit)
        cursor.close()
        lite.close()
    else:
        print(f"  (no SQLite mirror at {SQLITE_PATH}; skipped)")

    if audit:
        AUDIT.parent.mkdir(parents=True, exist_ok=True)
        with AUDIT.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(audit[0].keys()))
            writer.writeheader()
            writer.writerows(audit)
        print(f"\naudit written to {AUDIT}")

    for store, stats in totals.items():
        print(f"\n{store}")
        for key in ("agrees_exactly", "agrees_genus_only",
                    "agrees_on_the_epithet", "no_binomial_at_this_node",
                    "tax_id_not_in_dump", "DISAGREES"):
            print(f"  {key:30} {stats[key]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
