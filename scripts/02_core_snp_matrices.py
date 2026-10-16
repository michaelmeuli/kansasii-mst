"""Per-subset core-genome SNP distance matrices from snippy snps.aligned.fa (full-length alignment to the reference).
Core = sites with an unambiguous base (ACGT) in every genome of the subset; distance = Hamming count over variable core
sites. Positions in the mask BED (00_build_mask_bed.py: PE/PPE, IS, phage, self-repeats) are excluded from the core, like
snippy-core --mask; --no-mask skips this. No recombination filtering."""
from __future__ import annotations

import argparse

import numpy as np
import numpy.typing as npt
import pandas as pd
from pathlib import Path

from common import MASK_BED, OUT, REF_CHR, SUBSETS, alignment_path


def read_chr(path: Path) -> npt.NDArray[np.uint8]:
    seq: list[str] = []
    on = False
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            if on:
                break
            on = line[1:].split()[0] == REF_CHR
        elif on:
            seq.append(line.strip().upper())
    return np.frombuffer("".join(seq).encode(), dtype=np.uint8)


ap = argparse.ArgumentParser()
ap.add_argument("--no-mask", action="store_true", help="do not exclude the mask BED regions from the core")
args = ap.parse_args()


def load_mask(length: int) -> npt.NDArray[np.bool_]:
    """Boolean vector (index = 0-based reference position) of masked positions on REF_CHR."""
    m = np.zeros(length, dtype=bool)
    for line in open(MASK_BED, encoding="utf-8"):
        f = line.split("\t")
        if len(f) >= 3 and f[0] == REF_CHR:
            m[int(f[1]):int(f[2])] = True
    return m


summary: list[dict[str, object]] = []
s = pd.read_csv(OUT / "samples_used.tsv", sep="\t", dtype={"subsets": str}, keep_default_na=False)
acgt = np.frombuffer(b"ACGT", dtype=np.uint8)
for name in SUBSETS:
    names = s[s.subsets.str.split(",").apply(lambda x: name in x)]["sample"].tolist()
    A = np.stack([read_chr(alignment_path(n)) for n in names])
    core_all = np.isin(A, acgt).all(axis=0)
    var_all = int((A[:, core_all] != A[0, core_all]).any(axis=0).sum())
    masked = np.zeros(A.shape[1], dtype=bool) if args.no_mask else load_mask(A.shape[1])
    core = core_all & ~masked
    C = A[:, core]
    var = (C != C[0]).any(axis=0)
    V = C[:, var]
    print(f"{name}: {len(names)} genomes, core {core.sum()} sites ({100 * core.mean():.1f}%), variable {V.shape[1]} "
          f"(unmasked: core {core_all.sum()}, variable {var_all}; masked bp {masked.sum()})", flush=True)
    summary.append({"subset": name, "genomes": len(names), "masked_bp": int(masked.sum()), "core_sites_unmasked": int(core_all.sum()),
                    "variable_unmasked": var_all, "core_sites": int(core.sum()), "variable": V.shape[1]})
    D = np.array([(V != V[i]).sum(axis=1) for i in range(len(names))])
    pd.DataFrame(D, index=names, columns=names).to_csv(OUT / f"core_snp_distance_{name}.csv")
    pd.DataFrame({"pos": np.flatnonzero(core)[var] + 1}).to_csv(OUT / f"core_snp_sites_{name}.csv", index=False)
pd.DataFrame(summary).to_csv(OUT / "core_snp_summary.tsv", sep="\t", index=False)
