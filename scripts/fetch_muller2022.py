#!/usr/bin/env python3
"""Fetch the Borenstein lab's curated microbiome-metabolome collection.

Muller, Algavi & Borenstein 2022, npj Biofilms and Microbiomes 8:79
"The gut microbiome-metabolome dataset collection: a curated resource for
integrative meta-analysis" (doi:10.1038/s41522-022-00345-5)

14 human faecal cohorts, 2,900 samples from 1,849 individuals, with paired
microbiome and metabolome profiles. The genus tables are already relative
abundances summing to 1 per sample, which is this database's rank invariant, and
they arrive as plain TSV rather than through Bioconductor the way
curatedMetagenomicData's profiles do.

Four tables per study are fetched and the rest are skipped:

  genera.tsv    genus relative abundances, samples as rows
  metadata.tsv  sample and subject characteristics, including study group
  mtb.tsv       metabolite levels, samples as rows
  mtb.map.tsv   original metabolite names -> KEGG / HMDB, with a
                High.Confidence.Annotation flag that must not be dropped

species.tsv is left behind deliberately: it exists for only 6 of the 14 cohorts,
and genus is the common denominator across 16S and shotgun anyway. The
*.counts.tsv tables are only needed to re-derive abundances that are already
provided. Together those would be another 535MB.

Two things make this collection harder to load than its format suggests, and
neither is handled here -- this script only fetches.

  Taxonomy is GTDB, not NCBI. Profiles came from kraken2/bracken against GTDB,
  which carries alphabetic suffixes for polyphyletic groups (Clostridium_A,
  Firmicutes_A) and circumscribes genera differently from NCBI. taxa is keyed on
  NCBI tax IDs, so a mapping with an explicit reject list has to come before any
  insert, or resolve_name will mint duplicates and inventions.

  Metabolite levels are not comparable across studies. Instruments, targeted
  versus untargeted designs and units all differ, and nothing sums to 1, so the
  levels cannot be pooled the way relative abundances can.

Downloads land in data/muller2022/<STUDY>/ which is gitignored, following
data/microbiomehd/: upstream data, fetched rather than redistributed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

REPO = "borenstein-lab/microbiome-metabolome-curated-data"
CONTENTS = f"https://api.github.com/repos/{REPO}/contents/data/processed_data"
DEST = Path("data/muller2022")
WANTED = ("genera.tsv", "metadata.tsv", "mtb.tsv", "mtb.map.tsv")


def session() -> requests.Session:
    s = requests.Session()
    s.headers.update({"User-Agent": "gutdb-muller-fetch/1.0"})
    s.mount("https://", HTTPAdapter(max_retries=Retry(
        total=4, backoff_factor=1.5, status_forcelist=(429, 500, 502, 503, 504))))
    return s


def main() -> int:
    s = session()
    listing = s.get(CONTENTS, timeout=90)
    listing.raise_for_status()
    studies = sorted(i["name"] for i in listing.json() if i["type"] == "dir")
    print(f"{len(studies)} studies in {REPO}")

    total = skipped = 0
    for study in studies:
        files = s.get(f"{CONTENTS}/{study}", timeout=90)
        files.raise_for_status()
        by_name = {i["name"]: i for i in files.json()}
        out = DEST / study
        out.mkdir(parents=True, exist_ok=True)
        got = []
        for name in WANTED:
            item = by_name.get(name)
            if item is None:
                got.append(f"{name}=absent")
                continue
            target = out / name
            if target.exists() and target.stat().st_size == item["size"]:
                skipped += 1
                got.append(f"{name}=cached")
                continue
            body = s.get(item["download_url"], timeout=300)
            body.raise_for_status()
            target.write_bytes(body.content)
            total += len(body.content)
            got.append(f"{name}={len(body.content):,}B")
        print(f"   {study:34} " + "  ".join(got))

    print(f"\nfetched {total/1e6:.1f}MB into {DEST}"
          + (f", {skipped} files already cached" if skipped else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
