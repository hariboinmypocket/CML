#!/usr/bin/env python3
"""Annotate taxa with MicrobiomeHD's non-specific response labels (file-S3).

Duvallet et al. 2017 call a genus part of the shared response to disease when it
is significant at q < 0.05 in the same direction in at least two DIFFERENT
diseases, labelling it 'health', 'disease', or 'mixed' (enriched in cases for two
diseases and in controls for two others).

This is not an association, so it does not belong in
taxon_disease_associations -- it is a statement about a genus across diseases,
not about a genus in one disease. It lands on taxa.nonspecific_response instead,
alongside human_pathogen, which is the same shape: a single source's trait call
where NULL means "not assessed" rather than a negative.

Its value is directional. v_taxon_specificity already reports that a genus is
pan-disease, but n_diseases carries no direction, so the pan_disease class mixes
health-associated generalists like Faecalibacterium with disease-associated ones
like Streptococcus. This column separates them.

The file lists 144 genera and labels 51; the unlabelled 93 are left NULL rather
than recorded as "specific", since the file does not say that. Five of the 51 are
RDP cluster labels rather than organisms (Clostridium_IV, Clostridium_XI,
Clostridium_XVIII, Clostridium_XlVb, Lachnospiracea_incertae_sedis) and one,
Ethanoligenens, is a real genus this database does not hold; none is created
here, because inventing a taxon to carry an annotation puts the annotation
first.

Dry run by default; pass --apply to write.
"""
from __future__ import annotations

import csv
import os
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gutdb.config import Settings
from gutdb.db import connect

DRY_RUN = "--apply" not in sys.argv
SOURCE = Path(os.environ.get("MICROBIOMEHD_DIR", "data/microbiomehd"))
LABELS = SOURCE / "file-S3.nonspecific_genera.txt"
AUDIT = Path("data/microbiomehd_s3_audit.csv")

RAW = "https://raw.githubusercontent.com/cduvallet/microbiomeHD/master"
VALID = {"health", "disease", "mixed"}


def main() -> int:
    if not LABELS.exists():
        raise SystemExit(
            f"{LABELS} not found. Fetch it with:\n\n"
            f"  mkdir -p {SOURCE}\n"
            f"  curl -sS -o {LABELS} \\\n"
            f"    {RAW}/final/supp-files/file-S3.nonspecific_genera.txt\n")

    rows = list(csv.reader(LABELS.open(), delimiter="\t"))[1:]
    labelled: dict[str, str] = {}
    listed = 0
    for row in rows:
        match = re.search(r"g__([A-Za-z0-9_\-\[\]]+)", row[0])
        if not match:
            continue
        listed += 1
        label = (row[1].strip().casefold() if len(row) > 1 else "")
        if label in VALID:
            labelled[match.group(1).strip("[]")] = label
    print(f"S3: {listed} genera listed, {len(labelled)} labelled")
    print(f"    {dict(Counter(labelled.values()))}")

    conn = connect(Settings.from_env(".env"))
    cur = conn.cursor()
    stats = Counter()
    audit: list[dict] = []

    for genus, label in sorted(labelled.items()):
        cur.execute(
            "SELECT id, scientific_name, nonspecific_response FROM taxa "
            "WHERE genus = %s AND species = '' AND taxonomic_rank = 'genus'",
            (genus,))
        hit = cur.fetchone()
        if not hit:
            stats["no_matching_taxon"] += 1
            audit.append({"genus": genus, "label": label, "taxon_id": "",
                          "action": "no matching genus row - not created"})
            continue
        taxon_id, name, current = hit
        if current == label:
            stats["already_set"] += 1
            continue
        stats["relabelled" if current else "annotated"] += 1
        audit.append({"genus": genus, "label": label, "taxon_id": taxon_id,
                      "action": f"{current or 'NULL'} -> {label}"})
        if not DRY_RUN:
            cur.execute("UPDATE taxa SET nonspecific_response = %s WHERE id = %s",
                        (label, taxon_id))

    if DRY_RUN:
        print("\n*** DRY RUN - nothing written. Re-run with --apply ***")
    else:
        conn.commit()

    for key in ("annotated", "relabelled", "already_set", "no_matching_taxon"):
        print(f"    {key:22} {stats[key]}")

    if audit:
        AUDIT.parent.mkdir(parents=True, exist_ok=True)
        with AUDIT.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=["genus", "label", "taxon_id", "action"])
            writer.writeheader()
            writer.writerows(audit)
        print(f"\naudit written to {AUDIT} ({len(audit)} rows)")

    cur.close()
    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
