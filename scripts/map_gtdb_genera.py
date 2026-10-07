#!/usr/bin/env python3
"""Map the Muller 2022 collection's GTDB genus labels onto NCBI taxa.

This runs before any loader and writes nothing to the database. It produces a
reviewable decision per label in data/gtdb_ncbi_genus_map.csv, because the
alternative -- handing GTDB strings to resolve_name -- is how Ruminococcus2,
Escherichia/shigella and Clostridium_sensu_stricto got minted as taxa.

The collection's genus tables are keyed by full GTDB lineage strings:

    d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Monoglobales;f__Firm-18;g__UBA1775

12,263 distinct labels appear across the 14 cohorts, which is roughly thirty
times the number of genera a human gut actually contains. Three things are mixed
together in them, and only the first belongs in this database:

  Named genera that NCBI also recognizes (89.2% of abundance mass). These map.

  Names GTDB coined for uncultured lineages -- Fimenecus, Copromonas,
  Cryptobacteroides, Scatomorpha (2.6% of mass, 1,956 distinct). They look like
  ordinary Latin binomials, which is exactly what makes them dangerous: a
  fuzzy resolver will attach them to whatever is nearest. Each is verified
  against the local NCBI taxdump and rejected when NCBI has never heard of it.

  Genome-bin accessions -- g__UBA1775, g__JABGPL01, g__VGIO01, g__SURF-9
  (8.1% of mass, 8,023 distinct). No NCBI equivalent exists or could.

One subtler hazard is reported but not resolved here. GTDB splits polyphyletic
NCBI genera and marks the fragments with letter suffixes, so where the data
contains both g__Clostridium and g__Clostridium_A, GTDB's bare Clostridium is a
PROPER SUBSET of NCBI's Clostridium rather than the same thing. Mapping it onto
the NCBI genus is defensible -- it is the fragment NCBI's type species sits in
-- but the abundance is not interchangeable with an NCBI-profiled Clostridium
column, so those labels are flagged gtdb_split=yes for a human to accept or
drop. The suffix itself carries no NCBI meaning and is stripped for matching.

Usage:
    python scripts/map_gtdb_genera.py              # write the map, report coverage
    python scripts/map_gtdb_genera.py --unmapped   # also list rejected names by mass
"""
from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect

csv.field_size_limit(10 ** 9)

SOURCE = Path("data/muller2022")
OUT = Path("data/gtdb_ncbi_genus_map.csv")
TAXDUMP = Path("taxdump/names.dmp")

# A genus name in the Linnaean sense: one capitalized word, optionally carrying
# a GTDB polyphyly suffix (Fusibacter_A). Anything else -- UBA1775, JABGPL01,
# CAG-81, 1XD42-69 -- is a genome bin accession.
NAMED = re.compile(r"^[A-Z][a-z]+(_[A-Z]+)?$")


def gtdb_genus(label: str) -> str:
    """The g__ field of a GTDB lineage string, or '' if there is none."""
    return label.rsplit(";g__", 1)[-1] if ";g__" in label else ""


def read_labels() -> tuple[dict[str, float], dict[str, set[str]]]:
    """Total abundance mass and cohort set per distinct GTDB label."""
    mass: dict[str, float] = defaultdict(float)
    cohorts: dict[str, set[str]] = defaultdict(set)
    for study_dir in sorted(SOURCE.glob("*")):
        path = study_dir / "genera.tsv"
        if not path.exists():
            continue
        study = study_dir.name
        with path.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)[1:]
            for col in header:
                cohorts[col].add(study)
            for row in reader:
                for col, value in zip(header, row[1:]):
                    if value not in ("", "NA"):
                        mass[col] += float(value)
    return mass, cohorts


def ncbi_genus_names() -> set[str]:
    """Every name NCBI records, lowercased.

    Read from the local taxdump rather than queried, so the check works offline
    and covers synonyms and old names, not just the genera this database
    happens to hold already.
    """
    if not TAXDUMP.exists():
        print(f"   ! {TAXDUMP} absent -- cannot verify rejections against NCBI")
        return set()
    names: set[str] = set()
    with TAXDUMP.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.split("\t|\t")
            if len(parts) > 1:
                names.add(parts[1].strip().lower())
    return names


def main() -> int:
    mass, cohorts = read_labels()
    print(f"distinct GTDB labels across the collection: {len(mass):,}")

    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    cur.execute("""SELECT LOWER(scientific_name), id FROM taxa
                   WHERE taxonomic_rank = 'genus'""")
    in_taxa = {name: tid for name, tid in cur.fetchall()}
    conn.close()
    print(f"genus taxa already in this database: {len(in_taxa):,}")

    ncbi = ncbi_genus_names()
    print(f"names in the local NCBI taxdump: {len(ncbi):,}")

    # Which base names GTDB has split: present both bare and suffixed.
    bare: set[str] = set()
    suffixed: set[str] = set()
    for label in mass:
        g = gtdb_genus(label)
        if not NAMED.match(g):
            continue
        if "_" in g:
            suffixed.add(g.split("_")[0])
        else:
            bare.add(g)
    split_genera = bare & suffixed
    print(f"NCBI genera that GTDB splits (present bare and suffixed): {len(split_genera)}")

    total_mass = sum(mass.values())
    rows = []
    by_decision: dict[str, list[float]] = defaultdict(list)
    for label, m in mass.items():
        g = gtdb_genus(label)
        base = g.split("_")[0] if g else ""
        if label == "Unclassified" or not g:
            decision, taxon_id, ncbi_name = "reject_unclassified", None, ""
        elif not NAMED.match(g):
            decision, taxon_id, ncbi_name = "reject_mag_bin", None, ""
        elif base.lower() in in_taxa:
            decision = "map"
            taxon_id, ncbi_name = in_taxa[base.lower()], base
        elif ncbi and base.lower() in ncbi:
            # NCBI knows the name but this database has no genus row for it yet.
            decision, taxon_id, ncbi_name = "map_needs_taxon_row", None, base
        else:
            decision, taxon_id, ncbi_name = "reject_gtdb_coined", None, ""
        rows.append({
            "gtdb_label": label,
            "gtdb_genus": g,
            "base_genus": base,
            "decision": decision,
            "ncbi_taxon_id": taxon_id or "",
            "ncbi_name": ncbi_name,
            "gtdb_split": "yes" if base in split_genera else "",
            "mass_share": round(m / total_mass, 9),
            "n_cohorts": len(cohorts[label]),
        })
        by_decision[decision].append(m)

    rows.sort(key=lambda r: -r["mass_share"])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n{'decision':26} {'labels':>8} {'mass share':>12}")
    for decision in sorted(by_decision, key=lambda d: -sum(by_decision[d])):
        ms = by_decision[decision]
        print(f"{decision:26} {len(ms):>8} {100*sum(ms)/total_mass:>11.2f}%")
    mapped = sum(sum(by_decision[d]) for d in ("map", "map_needs_taxon_row"))
    print(f"\nmappable abundance mass: {100*mapped/total_mass:.2f}%")
    split_mass = sum(r["mass_share"] for r in rows
                     if r["gtdb_split"] and r["decision"].startswith("map"))
    print(f"   of which flagged gtdb_split (narrower than the NCBI genus): "
          f"{100*split_mass:.2f}%")
    print(f"\nmap written to {OUT} ({len(rows):,} rows)")

    if "--unmapped" in sys.argv:
        print("\nrejected GTDB-coined names by mass:")
        for r in rows:
            if r["decision"] == "reject_gtdb_coined" and r["mass_share"] > 0.0002:
                print(f"   {r['gtdb_genus']:28} {100*r['mass_share']:>6.3f}%  "
                      f"{r['n_cohorts']} cohorts")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
