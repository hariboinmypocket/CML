#!/usr/bin/env python3
"""Give each rank='species' row the species tax ID, not a strain's.

283 species rows in MySQL and 103 in the mirror carry an `ncbi_tax_id` that is a
strain or subspecies node rather than the species it names: row 3
`Bifidobacterium longum` holds 1322347, which NCBI calls "B. longum subsp.
longum KACC 91563", under species 216816. The sources are strain-level --
MiMeDB records a deposited strain, GMrepo's markers carry whatever ID the study
used -- and nothing has ever reconciled the two.

It matters because the tax ID is this database's canonical key. Every repair
here resolves names, lineage and duplicates through it, so a strain ID makes a
row permanently unmatchable: `merge_duplicate_taxa.py` groups by shared tax ID
and cannot see that row 92 `Prevotella stercorea` and the row holding 363265 are
one organism; a reload inserts beside such a row rather than updating it, which
is the single case 665989e could not close; and NCBI's current name for a strain
ID keeps the strain designation, so the row reads as a rename candidate forever.
The fix is the `retaxid` action `merge_duplicate_taxa.py`'s OVERRIDES performs by
hand for F. prausnitzii and L. salivarius, and that two scripts here have since
performed 15 more times, applied as a rule.

It is only a rule because two guards make it one. Repointing to a species node
blindly is the F. prausnitzii mistake at scale:

  * `same_epithet`. NCBI moves strains between species, so the species node of a
    strain ID is not always the species the row names. Tax 1280705 sits under
    *Segatella baroniae* while the row says `Prevotella bryantii`; tax 1219013
    under a Rhodococcus subgenus while the row says `hoagii`; tax 216596 under
    *Rhizobium johnstonii* while the row says `leguminosarum`. Those five rows
    are reported, never written: each is a strain reidentified since the row was
    loaded, and following it would relabel the organism rather than correct its
    ID. Latin gender agreement still passes, as everywhere here.
  * Nothing else holds the target. 21 rows would land on a tax ID another row
    already carries, which makes them a merge -- two rows for one organism, with
    evidence on both -- and not a tax ID correction. They are listed in the audit
    as merge candidates for `merge_duplicate_taxa.py`, which can see them once
    somebody decides the merge is right.

Infraspecific rows are left alone outright: a row named `Kocuria marina subsp.
indica` or `... subsp. bulgaricus` legitimately carries a subspecies node, and
repointing it to the species would erase the distinction the name makes. The
epithet guard catches these too (`indica` is not `marina`), so the explicit check
is belt and braces for a row whose last token happens to match.

What this does NOT do is rename anything. Ten rows take the species ID while
their stored genus differs from NCBI's current one -- `Prevotella ruminicola`
onto 839, which NCBI calls *Xylanibacter ruminicola*. The ID is still right, and
a stale-but-real synonym name is a state `merge_duplicate_taxa.py` knows how to
canonicalize, whereas a strain ID is a state nothing can. Renaming is a separate
judgement and is left for one.

No evidence is relabelled: the row keeps its name, its rank and its dependents,
and a strain ID and its species node share every rank above species, so stored
lineage and the phylum reconciliation are untouched either way.

Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import csv
import os
import sqlite3
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb import taxdump
from gutdb.config import Settings
from gutdb.db import connect
from gutdb.transform import clean_str, parse_scientific_name, same_epithet

DRY_RUN = "--apply" not in sys.argv
TAXDUMP = os.environ.get("GUTDB_TAXDUMP", "taxdump")
SQLITE_PATH = os.environ.get("GUTDB_SQLITE", "gut_microbiome.sqlite")
AUDIT = Path("data/strain_taxid_audit.csv")
# A dry run writes beside the audit, never over it. Re-running in dry mode after
# --apply would otherwise replace the record of what was written with a list of
# the handful of rows still skipped -- which is exactly what happened the first
# time this ran, and the per-row record of 338 repointings was lost with it.
DRY_AUDIT = Path("data/strain_taxid_audit.dryrun.csv")

INFRASPECIFIC = (" subsp.", " subsp ", " var.", " var ", " biovar", " serovar", " pv.")


def species_names(dump) -> dict[int, str]:
    """Scientific names, for species-rank nodes only."""
    path = f"{TAXDUMP}/names.dmp"
    if not os.path.exists(path):
        raise SystemExit(
            f"taxdump not found at {TAXDUMP}/. Unzip taxdmp.zip there or set GUTDB_TAXDUMP.")
    names: dict[int, str] = {}
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split("\t|")
            if len(parts) < 4 or parts[3].strip() != "scientific name":
                continue
            taxid = int(parts[0])
            if dump.rank_of.get(taxid) == "species":
                names[taxid] = parts[1].strip()
    return names


def run_store(store: str, cursor, dump, names: dict[int, str],
              audit: list[dict]) -> Counter:
    placeholder = "%s" if store == "mysql" else "?"
    cursor.execute(
        "SELECT id, scientific_name, genus_key, species, ncbi_tax_id, taxonomic_rank "
        "FROM taxa WHERE ncbi_tax_id IS NOT NULL")
    rows = cursor.fetchall()
    held: Counter = Counter(int(r[4]) for r in rows if r[4])

    stats: Counter = Counter()
    for taxon_id, name, genus_key, species, taxid, rank in rows:
        if rank != "species":
            continue
        taxid = int(taxid)
        node = dump.species_node(taxid)
        if not node or node == taxid:
            continue
        stats["carries_a_strain_taxid"] += 1

        species_clean = clean_str(species)
        if any(marker in f" {species_clean.casefold()} " for marker in INFRASPECIFIC):
            stats["skipped_infraspecific"] += 1
            reason = "row names an infraspecific taxon; its node is correct"
        else:
            target_name = names.get(node, "")
            target_genus, target_epithet = parse_scientific_name(target_name)
            stored_epithet = species_clean.split()[-1] if species_clean else ""
            if not same_epithet(stored_epithet, target_epithet):
                stats["skipped_epithet_differs"] += 1
                reason = (f"species node {node} is {target_name!r}: the strain was "
                          f"reidentified, so this would relabel the organism")
            elif held.get(node, 0) >= 1:
                stats["skipped_target_held"] += 1
                reason = (f"tax ID {node} is already held by another row; this is a "
                          "merge for merge_duplicate_taxa.py, not a retaxid")
            else:
                note = ("stored genus differs from NCBI's "
                        f"({target_genus}); name left alone"
                        if target_genus.casefold() != genus_key else "")
                if not DRY_RUN:
                    cursor.execute(
                        f"UPDATE taxa SET ncbi_tax_id={placeholder} WHERE id={placeholder}",
                        (node, taxon_id))
                held[node] += 1
                stats["repointed"] += 1
                stats["repointed_genus_differs"] += bool(note)
                audit.append({
                    "store": store, "action": "repointed", "taxon_id": taxon_id,
                    "scientific_name": name, "from_tax_id": taxid, "to_tax_id": node,
                    "ncbi_name_for_target": target_name, "note": note,
                })
                continue
        audit.append({
            "store": store, "action": "skipped", "taxon_id": taxon_id,
            "scientific_name": name, "from_tax_id": taxid, "to_tax_id": node,
            "ncbi_name_for_target": names.get(node, ""), "note": reason,
        })
    return stats


def main() -> int:
    print(f"loading taxdump from {TAXDUMP}/ ...")
    dump = taxdump.load(TAXDUMP)
    if dump is None:
        raise SystemExit(f"taxdump not found at {TAXDUMP}/.")
    names = species_names(dump)
    print(f"  nodes {len(dump.rank_of):,}   species names {len(names):,}")

    audit: list[dict] = []
    totals = {}

    print("\nmysql")
    conn = connect(Settings.from_env(".env"))
    cursor = conn.cursor()
    totals["mysql"] = run_store("mysql", cursor, dump, names, audit)
    if not DRY_RUN:
        conn.commit()
    cursor.close()
    conn.close()

    print("\nsqlite")
    if Path(SQLITE_PATH).exists():
        lite = sqlite3.connect(SQLITE_PATH)
        cursor = lite.cursor()
        totals["sqlite"] = run_store("sqlite", cursor, dump, names, audit)
        if not DRY_RUN:
            lite.commit()
        cursor.close()
        lite.close()
    else:
        print(f"  (no SQLite mirror at {SQLITE_PATH}; skipped)")

    destination = DRY_AUDIT if DRY_RUN else AUDIT
    if audit:
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(audit[0].keys()))
            writer.writeheader()
            writer.writerows(audit)
        print(f"\naudit written to {destination}")
    else:
        print("\nno species row carries a strain tax ID - nothing to repoint")

    for store, stats in totals.items():
        print(f"\n{store}")
        for key in ("carries_a_strain_taxid", "repointed", "repointed_genus_differs",
                    "skipped_target_held", "skipped_epithet_differs",
                    "skipped_infraspecific"):
            print(f"  {key:30} {stats[key]}")

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
