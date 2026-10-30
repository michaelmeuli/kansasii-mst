"""Shared paths, QC exclusion and subset definitions."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from kansasii_paths import ASSEMBLIES, COMPLEX, OUT_MST, SCREENING_MAP_RESULTS

OUTPUT = OUT_MST / "results"
MASK_DIR = OUTPUT / "mask"
MASK_BED = MASK_DIR / "mask.bed"  # PE/PPE + IS + phage + self-repeats, merged (00_build_mask_bed.py)
SNIPPY_BIN = Path.home() / "data/conda/envs/snippy/bin"  # bedtools, minimap2
REF_CHR = "NC_022663.1"  # chromosome of the snippy reference (ATCC12478_ACCESSION)

# one MST per subset; None = every QC-passing complex member
SUBSETS: dict[str, list[str] | None] = {
    "complex": None,
    "kansasii": ["kansasii"],
    "pseudokansasii": ["pseudokansasii"],
    "persicum": ["persicum"],
}
# MST nodes joined by edges of <= N SNPs are merged into one node (03_mst_figures.py); N sits where the many short central
# edges end (see README)
MERGE_SNPS = {"complex": 1000, "kansasii": 50, "pseudokansasii": 20, "persicum": 100}


def alignment_path(sample: str, assemblies: Path = ASSEMBLIES) -> Path:
    return assemblies / sample / "5_typing/kansasii_snippy/snippy_out/snps.aligned.fa"


def load_samples(results: Path = SCREENING_MAP_RESULTS) -> pd.DataFrame:
    """All 186 genomes with columns sample, NR, species, keep (bool) and reason (why excluded, else '')."""
    m = pd.read_csv(results, dtype={"excluded": str},
                    usecols=["NR", "PROBENNUMMER", "LABEL", "species", "contamination_flag", "excluded"])
    reason = pd.Series("", index=m.index)
    reason[~m.species.isin(COMPLEX)] = "not in complex"
    reason[m.LABEL.notna() & (reason == "")] = "control strain"
    reason[(m.contamination_flag == "contaminated") & (reason == "")] = "contaminated"
    # `excluded` = contaminated or < 90 % of the cgMLST loci called (immensekansasii 06)
    reason[(m.excluded.str.lower() == "true") & (reason == "")] = "low cgMLST coverage"
    return pd.DataFrame({"sample": m.PROBENNUMMER, "NR": m.NR, "species": m.species,
                         "keep": reason == "", "reason": reason}).sort_values("NR").reset_index(drop=True)
