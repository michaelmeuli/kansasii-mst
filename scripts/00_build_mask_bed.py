"""BED mask of regions that give false SNPs when reads are mapped: PE/PPE genes, insertion sequences, phage/integrase genes
(from the reference GFF3 product names) and self-repeats (the recipe of `snippy-core --mask auto`: minimap2 self-alignment).
Snippy itself has no --mask; snippy-core does, with the same BED semantics (0-based start, end exclusive). The mask is applied
in 02_core_snp_matrices.py. Output: output/mst/mask/{pe_ppe,is,phage,self_repeat,mask}.bed"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
from collections import defaultdict
from pathlib import Path

from common import MASK_BED, MASK_DIR, SNIPPY_BIN
from kansasii_paths import ATCC12478_FASTA, ATCC12478_GFF

CLASSES = {
    "pe_ppe": re.compile(r"PPE|PE family|PE domain|PE-PGRS|PGRS"),
    "is": re.compile(r"ransposase|nsertion sequence|IS[0-9]+"),
    "phage": re.compile(r"[Pp]hage|[Ii]ntegrase|[Pp]rophage|[Tt]erminase|[Pp]ortal|[Cc]apsid|[Ee]xcisionase"),
}
NOT_PHAGE = re.compile(r"phage shock|phage resistance")  # stress-response genes, not prophage

ap = argparse.ArgumentParser()
ap.add_argument("--ref-fasta", type=Path, default=ATCC12478_FASTA, help="snippy reference genome")
ap.add_argument("--ref-gff", type=Path, default=ATCC12478_GFF, help="its annotation")
ap.add_argument("--pad", type=int, default=0, help="extend every annotated feature by this many bp on both sides")
args = ap.parse_args()
os.environ["PATH"] = f"{SNIPPY_BIN}:{os.environ['PATH']}"
MASK_DIR.mkdir(parents=True, exist_ok=True)

lengths: dict[str, int] = {}
name = ""
for line in open(args.ref_fasta, encoding="utf-8"):
    if line.startswith(">"):
        name = line[1:].split()[0]
        lengths[name] = 0
    else:
        lengths[name] += len(line.strip())
genome = MASK_DIR / "ref.genome"
genome.write_text("".join(f"{k}\t{v}\n" for k, v in lengths.items()), encoding="utf-8")

# annotation classes: CDS (incl. pseudogene CDS) whose product name matches
rows: dict[str, list[tuple[str, int, int, str]]] = defaultdict(list)
for line in open(args.ref_gff, encoding="utf-8"):
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "CDS":
        continue
    attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    prod = attr.get("product", "")
    for cls, rx in CLASSES.items():
        if rx.search(prod) and not (cls == "phage" and NOT_PHAGE.search(prod)):
            start = max(int(f[3]) - 1 - args.pad, 0)
            rows[cls].append((f[0], start, min(int(f[4]) + args.pad, lengths[f[0]]), attr.get("locus_tag", "") + "|" + cls))
for cls in CLASSES:
    with open(MASK_DIR / f"{cls}.bed", "w", encoding="utf-8") as fh:
        for r in sorted(rows[cls]):
            fh.write("\t".join(map(str, r)) + "\n")

# self-repeats, as snippy-core --mask auto does
subprocess.run(
    f"minimap2 -D -P -w19 -m200 '{args.ref_fasta}' '{args.ref_fasta}' 2> /dev/null | cut -f 1,3,4 | sort -k1,1 -k2,2n"
    f" | bedtools slop -i stdin -g '{genome}' -r 1 -l 0 | bedtools merge -i stdin > '{MASK_DIR / 'self_repeat.bed'}'",
    shell=True, check=True, executable="/bin/bash")
beds = " ".join(f"'{MASK_DIR / (c + '.bed')}'" for c in [*CLASSES, "self_repeat"])
subprocess.run(f"cat {beds} | cut -f 1-3 | sort -k1,1 -k2,2n | bedtools merge -i stdin > '{MASK_BED}'",
               shell=True, check=True, executable="/bin/bash")


def bp(path: str) -> int:
    return sum(int(x[2]) - int(x[1]) for x in (ln.split("\t") for ln in open(path, encoding="utf-8")) if len(x) >= 3)


total = lengths["NC_022663.1"]
for cls in [*CLASSES, "self_repeat", "mask"]:
    p = str(MASK_BED if cls == "mask" else MASK_DIR / f"{cls}.bed")
    n = sum(1 for _ in open(p, encoding="utf-8"))
    print(f"{cls:12s} {n:5d} regions {bp(p):9d} bp ({100 * bp(p) / sum(lengths.values()):.1f}% of reference)")
print("chromosome", total, "bp; plasmid and others", sum(lengths.values()) - total, "bp")
