#!/usr/bin/env python3
"""Merge taxa rows that share one NCBI tax ID into a single canonical row.

115 tax IDs were spread across 241 rows. The cause is name-variant drift: the
same organism arrives from different sources as a formatting variant
("Oribacterium sp. oral taxon 102" / "Oribacterium sp oral taxon 102"), a
parser artifact that kept both the old and new genus ("Lactobacillus
limosilactobacillus reuteri"), or a nomenclature synonym ("Ruminococcus
gnavus" / "Mediterraneibacter gnavus").

Two things go wrong because of it:

1. Evidence and measurements land on different rows. "Prevotella copri" holds
   4 associations and no abundances; "Segatella copri" holds 33 associations
   and 4,704 abundances. Joining the two tables finds no overlap for one of
   the most-studied organisms in the gut.
2. Eight pairs were loaded into the SAME sample twice, so one organism is
   counted twice in 558 samples. Every overlapping pair was verified to hold
   an identical value, which confirms it is one measurement stored under two
   spellings rather than two distinct source rows. The per-rank sums are
   nonetheless exactly 1.0, because renormalization happened after the
   duplication -- so `validate` cannot see this, and the duplicated organism
   holds double share while every other taxon in the sample is squeezed.

Canonicalization follows NCBI, which is the point of storing ncbi_tax_id:
the row whose name matches NCBI's current scientific name for that tax ID
wins, and the rest are merged into it. Where no row matches, the row carrying
the most data wins and is renamed to NCBI's current name.

Classification precedence puts a differing species epithet ahead of an NCBI
name match, because a renamed *strain* is indistinguishable from a renamed
*species* by name alone. NCBI tax ID 411483 is "Faecalibacterium
prausnitzii A2-165", renamed F. duncaniae in 2022, while the species
F. prausnitzii remains tax ID 853. Merging on the name match would have
relabelled 11,048 abundance rows of the most abundant commensal here as a
different organism. Groups like that are listed in OVERRIDES instead.

Attributes are backfilled, not dropped: where the canonical row is empty and
a merged row holds a value, the value moves across, so MiMeDB traits and
lineage sitting on a losing row survive the merge.

Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import csv
import os
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect
from gutdb.ncbi import TaxonomyClient
from gutdb.transform import parse_scientific_name, taxon_key

DRY_RUN = "--apply" not in sys.argv
AUDIT = Path("data/taxa_merge_audit.csv")
# NCBI asks E-utilities callers to identify themselves; requests still succeed
# without it, just under tighter rate limits. Supplied the same way the CLI
# takes it (--email) rather than hardcoded, so no personal address is tracked.
EMAIL = os.environ.get("NCBI_EMAIL") or next(
    (arg.split("=", 1)[1] for arg in sys.argv if arg.startswith("--email=")), None
)

# Columns moved from a losing row onto the canonical row when the canonical is
# empty. Lineage and traits only -- never keys, names, rank or tax ID.
TRAIT_COLUMNS = (
    "superkingdom", "phylum", "class_name", "order_name", "family",
    "gram_stain", "oxygen_requirement", "ph_preference", "sporulation",
    "shape", "cell_arrangement", "mobility", "flagella_presence",
    "number_of_membranes", "biotic_relationship", "habitat",
    "temperature_range", "optimal_temperature", "metabolism",
    "energy_source", "human_pathogen",
)
# These default to the string 'N/A' rather than NULL, so emptiness differs.
NA_COLUMNS = ("energy_mode", "primary_food_source")

# Groups the automatic rules flag for review, resolved by hand. Each entry is
# keyed by tax ID; see the module docstring for the reasoning.
OVERRIDES: dict[int, dict] = {
    # F. prausnitzii carries a strain tax ID that was renamed out from under
    # it. Not a duplicate: repoint the species row to the species tax ID and
    # leave F. duncaniae standing as its own organism.
    411483: {"action": "retaxid", "taxon_id": 17, "new_taxid": 853,
             "note": "species F. prausnitzii mis-assigned the A2-165 strain taxid"},
    # Both rows are dead old-name artifacts with no data. A correctly named
    # Ligilactobacillus salivarius already exists (row 1287) but carries the
    # strain tax ID 1423799 instead of the species tax ID.
    1624: {"action": "delete_and_retaxid", "delete": [326, 2842],
           "taxon_id": 1287, "new_taxid": 1624,
           "note": "drop old-name artifacts; give the kept row the species taxid"},
    # Conservative flags that are genuine synonyms after checking: merge as
    # the automatic rule proposed.
    187101: {"action": "merge", "note": "S. amnii renamed S. vaginalis (2021)"},
    2585118: {"action": "merge", "note": "loser is the 'bacteroidales bacterium ph8' artifact"},
    29361: {"action": "merge", "note": "C. nexile -> Tyzzerella -> Faecalimonas nexilis"},
    2093823: {"action": "merge", "note": "'Candidatus' is a status prefix, not a name"},
    33039: {"action": "merge", "note": "genuine 2020 Ruminococcus -> Mediterraneibacter move"},
}

INFLECTIONS = ("us", "um", "a", "is", "e", "ii", "i", "ae", "os", "on", "s", "ys")


def norm(name: str) -> str:
    s = name.casefold()
    s = re.sub(r"\(.*?\)", " ", s)
    s = s.replace("[", " ").replace("]", " ")
    s = re.sub(r"[._\-]", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def stem(word: str) -> str:
    for suffix in sorted(INFLECTIONS, key=len, reverse=True):
        if word.endswith(suffix) and len(word) - len(suffix) >= 4:
            return word[: -len(suffix)]
    return word


def same_epithet(a: str, b: str) -> bool:
    """True when two epithets differ only by Latin gender agreement.

    'ramosum'/'ramosa' and 'crossotus'/'crossota' are the same organism under a
    new genus; 'prausnitzii'/'duncaniae' are not.
    """
    if a == b:
        return True
    sa, sb = stem(a), stem(b)
    return len(sa) >= 4 and sa == sb


def load_groups(cur) -> dict[int, list[dict]]:
    cur.execute(
        """
        SELECT t.ncbi_tax_id, t.id, t.scientific_name, t.taxonomic_rank, t.source_database,
          (SELECT COUNT(*) FROM taxon_disease_associations a WHERE a.taxon_id = t.id),
          (SELECT COUNT(*) FROM sample_taxon_abundances b WHERE b.taxon_id = t.id)
        FROM taxa t
        WHERE t.ncbi_tax_id IN (
            SELECT ncbi_tax_id FROM taxa
            WHERE ncbi_tax_id IS NOT NULL GROUP BY 1 HAVING COUNT(*) > 1)
        ORDER BY t.ncbi_tax_id, t.id
        """
    )
    groups: dict[int, list[dict]] = {}
    for taxid, tid, sci, rank, src, n_assoc, n_abund in cur.fetchall():
        groups.setdefault(taxid, []).append(
            {"id": tid, "sci": sci, "rank": rank, "src": src,
             "assoc": n_assoc, "abund": n_abund, "norm": norm(sci)}
        )
    return groups


def fetch_ncbi_names(taxids: list[int]) -> dict[str, dict[str, str]]:
    client = TaxonomyClient(email=EMAIL)
    out: dict[str, dict[str, str]] = {}
    ids = [str(t) for t in taxids]
    for i in range(0, len(ids), 50):
        out.update(client.fetch_lineages(ids[i : i + 50]))
        time.sleep(0.4)
    return out


def plan(groups: dict[int, list[dict]], ncbi: dict) -> tuple[list[dict], list[dict]]:
    """Return (merge_plan, skipped). Each merge entry names a canonical and its losers."""
    merges: list[dict] = []
    skipped: list[dict] = []
    for taxid, rows in sorted(groups.items()):
        override = OVERRIDES.get(taxid, {})
        action = override.get("action", "merge")
        meta = ncbi.get(str(taxid), {})
        current_name = meta.get("_scientific_name", "")
        n_current = norm(current_name)

        if action != "merge":
            skipped.append({"taxid": taxid, **override})
            continue

        matches = [r for r in rows if r["norm"] == n_current]
        artifact = set()
        for a in rows:
            tail = " ".join(a["norm"].split()[1:])
            for b in rows:
                if b is not a and tail and tail == b["norm"]:
                    artifact.add(a["id"])

        epithets = [r["norm"].split()[-1] for r in rows if r["norm"]]
        conflict = any(not same_epithet(epithets[0], e) for e in epithets[1:])
        if conflict and len({r["norm"] for r in rows}) > 1 and taxid not in OVERRIDES:
            skipped.append({"taxid": taxid, "action": "unreviewed_epithet_conflict",
                            "note": " | ".join(r["sci"] for r in rows)})
            continue

        if matches:
            canon = max(matches, key=lambda r: r["assoc"] + r["abund"])
        else:
            pool = [r for r in rows if r["id"] not in artifact] or rows
            canon = max(pool, key=lambda r: (r["assoc"] + r["abund"], -r["id"]))

        losers = [r for r in rows if r["id"] != canon["id"]]
        if not losers:
            continue
        merges.append({
            "taxid": taxid,
            "canonical": canon,
            "losers": losers,
            "rename_to": current_name if (canon["norm"] != n_current and current_name) else "",
            "note": override.get("note", ""),
        })
    return merges, skipped


def backfill_traits(cur, canonical_id: int, loser_ids: list[int]) -> int:
    """Move lineage/trait values from losing rows onto an empty canonical field."""
    columns = TRAIT_COLUMNS + NA_COLUMNS
    cols = ", ".join(f"`{c}`" for c in columns)
    cur.execute(f"SELECT {cols} FROM taxa WHERE id = %s", (canonical_id,))
    canon = dict(zip(columns, cur.fetchone()))

    placeholders = ",".join(["%s"] * len(loser_ids))
    cur.execute(
        f"SELECT {cols} FROM taxa WHERE id IN ({placeholders}) ORDER BY id",
        tuple(loser_ids),
    )
    donors = [dict(zip(columns, row)) for row in cur.fetchall()]

    updates: dict[str, object] = {}
    for column in columns:
        empty = canon[column] is None or (
            column in NA_COLUMNS and str(canon[column]).strip() in ("", "N/A")
        )
        if not empty:
            continue
        for donor in donors:
            value = donor[column]
            if value is None:
                continue
            if column in NA_COLUMNS and str(value).strip() in ("", "N/A"):
                continue
            updates[column] = value
            break

    if updates:
        assignments = ", ".join(f"`{c}` = %s" for c in updates)
        cur.execute(
            f"UPDATE taxa SET {assignments} WHERE id = %s",
            (*updates.values(), canonical_id),
        )
    return len(updates)


def rename_canonical(cur, canonical_id: int, new_name: str) -> str:
    """Rename a canonical row to NCBI's current name, keys included.

    Returns a status string; a name already taken by another row is reported
    and skipped rather than raising on the unique key.
    """
    genus, species = parse_scientific_name(new_name)
    if not genus:
        return "skipped:unparseable"
    genus_key, species_key = taxon_key(genus, species)
    cur.execute(
        "SELECT id FROM taxa WHERE genus_key = %s AND species_key = %s AND id <> %s",
        (genus_key, species_key, canonical_id),
    )
    clash = cur.fetchone()
    if clash:
        return f"skipped:name_taken_by_{clash[0]}"
    cur.execute(
        """UPDATE taxa SET genus = %s, species = %s, genus_key = %s, species_key = %s,
                  scientific_name = %s WHERE id = %s""",
        (genus, species, genus_key, species_key, new_name, canonical_id),
    )
    return "renamed"


def main() -> int:
    settings = Settings.from_env(".env")
    conn = connect(settings)
    cur = conn.cursor()

    cur.execute("SELECT COUNT(*) FROM taxa")
    taxa_before = cur.fetchone()[0]
    cur.execute("SELECT COUNT(*) FROM sample_taxon_abundances")
    abund_before = cur.fetchone()[0]
    cur.execute("SELECT COUNT(*) FROM taxon_disease_associations")
    assoc_before = cur.fetchone()[0]

    groups = load_groups(cur)
    print(f"duplicate tax IDs: {len(groups)}  rows involved: {sum(len(v) for v in groups.values())}")
    ncbi = fetch_ncbi_names(list(groups))
    print(f"NCBI names resolved: {len(ncbi)}/{len(groups)}")

    merges, skipped = plan(groups, ncbi)
    print(f"merge groups: {len(merges)}   handled separately: {len(skipped)}")

    audit: list[dict] = []
    stats = {"taxa_deleted": 0, "assoc_repointed": 0, "assoc_collisions_dropped": 0,
             "abund_repointed": 0, "abund_duplicates_dropped": 0,
             "traits_backfilled": 0, "renamed": 0, "retaxid": 0}
    touched_samples: set[int] = set()

    for entry in merges:
        canon = entry["canonical"]
        loser_ids = [r["id"] for r in entry["losers"]]
        placeholders = ",".join(["%s"] * len(loser_ids))

        for loser in entry["losers"]:
            audit.append({
                "taxid": entry["taxid"], "action": "merged_away",
                "taxon_id": loser["id"], "scientific_name": loser["sci"],
                "rank": loser["rank"], "source": loser["src"],
                "assoc": loser["assoc"], "abund": loser["abund"],
                "into_taxon_id": canon["id"], "into_name": canon["sci"],
                "rename_to": entry["rename_to"], "note": entry["note"],
            })

        if DRY_RUN:
            stats["taxa_deleted"] += len(loser_ids)
            stats["assoc_repointed"] += sum(r["assoc"] for r in entry["losers"])
            stats["abund_repointed"] += sum(r["abund"] for r in entry["losers"])
            continue

        stats["traits_backfilled"] += backfill_traits(cur, canon["id"], loser_ids)

        # Record the samples that hold both rows before the duplicates go, so
        # only those rank groups are renormalized afterwards.
        cur.execute(
            f"""SELECT DISTINCT a.sample_id FROM sample_taxon_abundances a
                JOIN sample_taxon_abundances b
                  ON b.sample_id = a.sample_id AND b.taxon_id = %s
                WHERE a.taxon_id IN ({placeholders})""",
            (canon["id"], *loser_ids),
        )
        touched_samples.update(r[0] for r in cur.fetchall())

        # Drop a losing association only where the canonical already holds the
        # same (disease, comparison, effect_type); repoint the rest.
        cur.execute(
            f"""DELETE a FROM taxon_disease_associations a
                JOIN taxon_disease_associations b
                  ON b.taxon_id = %s AND b.disease_id = a.disease_id
                     AND b.comparison_id = a.comparison_id
                     AND b.effect_type = a.effect_type
                WHERE a.taxon_id IN ({placeholders})""",
            (canon["id"], *loser_ids),
        )
        stats["assoc_collisions_dropped"] += cur.rowcount
        cur.execute(
            f"UPDATE taxon_disease_associations SET taxon_id = %s WHERE taxon_id IN ({placeholders})",
            (canon["id"], *loser_ids),
        )
        stats["assoc_repointed"] += cur.rowcount

        # Same shape for abundances. Overlapping values were verified identical,
        # so dropping the losing copy loses no measurement.
        cur.execute(
            f"""DELETE a FROM sample_taxon_abundances a
                JOIN sample_taxon_abundances b
                  ON b.sample_id = a.sample_id AND b.taxon_id = %s
                WHERE a.taxon_id IN ({placeholders})""",
            (canon["id"], *loser_ids),
        )
        stats["abund_duplicates_dropped"] += cur.rowcount
        cur.execute(
            f"UPDATE sample_taxon_abundances SET taxon_id = %s WHERE taxon_id IN ({placeholders})",
            (canon["id"], *loser_ids),
        )
        stats["abund_repointed"] += cur.rowcount

        cur.execute(f"DELETE FROM taxa WHERE id IN ({placeholders})", tuple(loser_ids))
        stats["taxa_deleted"] += cur.rowcount

        if entry["rename_to"]:
            status = rename_canonical(cur, canon["id"], entry["rename_to"])
            audit.append({
                "taxid": entry["taxid"], "action": f"rename_canonical:{status}",
                "taxon_id": canon["id"], "scientific_name": canon["sci"],
                "rank": canon["rank"], "source": canon["src"],
                "assoc": canon["assoc"], "abund": canon["abund"],
                "into_taxon_id": canon["id"], "into_name": entry["rename_to"],
                "rename_to": entry["rename_to"], "note": entry["note"],
            })
            if status == "renamed":
                stats["renamed"] += 1

    # Hand-resolved groups: tax ID corrections and dead-artifact deletions.
    for entry in skipped:
        action = entry.get("action")
        if action == "retaxid":
            audit.append({"taxid": entry["taxid"], "action": "retaxid",
                          "taxon_id": entry["taxon_id"], "scientific_name": "",
                          "rank": "", "source": "", "assoc": "", "abund": "",
                          "into_taxon_id": entry["taxon_id"],
                          "into_name": f"ncbi_tax_id -> {entry['new_taxid']}",
                          "rename_to": "", "note": entry["note"]})
            if not DRY_RUN:
                cur.execute("UPDATE taxa SET ncbi_tax_id = %s WHERE id = %s",
                            (entry["new_taxid"], entry["taxon_id"]))
                stats["retaxid"] += cur.rowcount
        elif action == "delete_and_retaxid":
            placeholders = ",".join(["%s"] * len(entry["delete"]))
            for dead in entry["delete"]:
                audit.append({"taxid": entry["taxid"], "action": "deleted_dead_artifact",
                              "taxon_id": dead, "scientific_name": "", "rank": "",
                              "source": "", "assoc": 0, "abund": 0,
                              "into_taxon_id": entry["taxon_id"], "into_name": "",
                              "rename_to": "", "note": entry["note"]})
            audit.append({"taxid": entry["taxid"], "action": "retaxid",
                          "taxon_id": entry["taxon_id"], "scientific_name": "",
                          "rank": "", "source": "", "assoc": "", "abund": "",
                          "into_taxon_id": entry["taxon_id"],
                          "into_name": f"ncbi_tax_id -> {entry['new_taxid']}",
                          "rename_to": "", "note": entry["note"]})
            if not DRY_RUN:
                cur.execute(f"DELETE FROM taxa WHERE id IN ({placeholders})",
                            tuple(entry["delete"]))
                stats["taxa_deleted"] += cur.rowcount
                cur.execute("UPDATE taxa SET ncbi_tax_id = %s WHERE id = %s",
                            (entry["new_taxid"], entry["taxon_id"]))
                stats["retaxid"] += cur.rowcount
        else:
            print(f"  !! unresolved group {entry['taxid']}: {entry.get('note','')}")

    # Removing a duplicated copy drops that rank group below 1.0. The stored
    # values were normalized while the duplicate was present, so rescaling the
    # affected groups restores the true proportions. Only touched groups are
    # rescaled, to avoid perturbing clean samples with float noise.
    renormalized = 0
    if not DRY_RUN and touched_samples:
        for sample_id in sorted(touched_samples):
            for rank in ("genus", "species"):
                cur.execute(
                    """SELECT SUM(a.relative_abundance) FROM sample_taxon_abundances a
                       JOIN taxa t ON t.id = a.taxon_id
                       WHERE a.sample_id = %s AND t.taxonomic_rank = %s""",
                    (sample_id, rank),
                )
                total = cur.fetchone()[0]
                if not total or abs(float(total) - 1.0) < 1e-9:
                    continue
                cur.execute(
                    """UPDATE sample_taxon_abundances a
                       JOIN taxa t ON t.id = a.taxon_id
                       SET a.relative_abundance = a.relative_abundance / %s
                       WHERE a.sample_id = %s AND t.taxonomic_rank = %s""",
                    (float(total), sample_id, rank),
                )
                renormalized += 1
    stats["rank_groups_renormalized"] = renormalized

    AUDIT.parent.mkdir(parents=True, exist_ok=True)
    with AUDIT.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(audit[0].keys()))
        writer.writeheader()
        writer.writerows(audit)

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    else:
        conn.commit()

    print()
    for key, value in stats.items():
        print(f"  {key:30} {value}")

    cur.execute("SELECT COUNT(*) FROM taxa")
    taxa_after = cur.fetchone()[0]
    cur.execute("SELECT COUNT(*) FROM sample_taxon_abundances")
    abund_after = cur.fetchone()[0]
    cur.execute("SELECT COUNT(*) FROM taxon_disease_associations")
    assoc_after = cur.fetchone()[0]
    print(f"\n  taxa          {taxa_before} -> {taxa_after}")
    print(f"  abundances    {abund_before} -> {abund_after}")
    print(f"  associations  {assoc_before} -> {assoc_after}")
    print(f"  audit written to {AUDIT}")

    cur.close()
    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
