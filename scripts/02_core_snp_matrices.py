"""Per-subset core-genome SNP distance matrices from snippy snps.aligned.fa (full-length alignment to the reference).
Core = sites with an unambiguous base (ACGT) in every genome of the subset; distance = Hamming count over variable core
sites. No recombination filtering."""
from __future__ import annotations

import numpy as np
import numpy.typing as npt
import pandas as pd
from pathlib import Path

from common import OUT, REF_CHR, SUBSETS, alignment_path


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


s = pd.read_csv(OUT / "samples_used.tsv", sep="\t", dtype={"subsets": str}, keep_default_na=False)
acgt = np.frombuffer(b"ACGT", dtype=np.uint8)
for name in SUBSETS:
    names = s[s.subsets.str.split(",").apply(lambda x: name in x)]["sample"].tolist()
    A = np.stack([read_chr(alignment_path(n)) for n in names])
    core = np.isin(A, acgt).all(axis=0)
    C = A[:, core]
    var = (C != C[0]).any(axis=0)
    V = C[:, var]
    print(f"{name}: {len(names)} genomes, core {core.sum()} sites ({100 * core.mean():.1f}%), variable {V.shape[1]}", flush=True)
    D = np.array([(V != V[i]).sum(axis=1) for i in range(len(names))])
    pd.DataFrame(D, index=names, columns=names).to_csv(OUT / f"core_snp_distance_{name}.csv")
    pd.DataFrame({"pos": np.flatnonzero(core)[var] + 1}).to_csv(OUT / f"core_snp_sites_{name}.csv", index=False)
