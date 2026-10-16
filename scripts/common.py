"""Shared paths, QC exclusion and subset definitions."""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

ROOT = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii"))
OUT = ROOT / "output/mst"
RES = ROOT / "runs/mkan329/assembly/results"
REF_DIR = ROOT / "data/gtdb_genomes/Mycobacteriaceae/ncbi_dataset/data/GCF_000157895.3"
REF_FASTA = REF_DIR / "GCF_000157895.3_ASM15789v2_genomic.fna"
REF_GFF = REF_DIR / "genomic.gff"
MASK_DIR = OUT / "mask"
MASK_BED = MASK_DIR / "mask.bed"  # PE/PPE + IS + phage + self-repeats, merged (00_build_mask_bed.py)
SNIPPY_BIN = Path.home() / "data/conda/envs/snippy/bin"  # bedtools, minimap2
REF_CHR = "NC_022663.1"  # chromosome of the snippy reference GCF_000157895.3

COMPLEX = ["kansasii", "persicum", "pseudokansasii", "attenuatum", "innocens", "ostraviense"]
# one MST per subset; None = every QC-passing complex member
SUBSETS: dict[str, list[str] | None] = {
    "complex": None,
    "kansasii": ["kansasii"],
    "pseudokansasii": ["pseudokansasii"],
    "persicum": ["persicum"],
}
# MST nodes joined by edges of <= N SNPs are merged into one node (03_mst_figures.py); N sits where the many short central
# edges end (see README)
MERGE_SNPS = {"complex": 50, "kansasii": 50, "pseudokansasii": 20, "persicum": 100}
# same colours as mlsa-kansasii/mlsa/__init__.py
SPECIES_COLORS = {
    "kansasii": "#1b9e77", "persicum": "#d95f02", "pseudokansasii": "#7570b3", "innocens": "#e7298a",
    "attenuatum": "#66a61e", "ostraviense": "#e6ab02", "gastri": "#a6761d",
}
# low-quality genomes missed by the contamination flag (also used in kansasii-gwas)
EXTRA_EXCL = {10, 21}


def alignment_path(sample: str) -> Path:
    return RES / sample / "5_typing/kansasii_snippy/snippy_out/snps.aligned.fa"


def load_samples() -> pd.DataFrame:
    """All 186 genomes with columns sample, NR, species, keep (bool) and reason (why excluded, else '')."""
    m = pd.read_csv(ROOT / "data/imm/screening_map_results.csv",
                    usecols=["NR", "PROBENNUMMER", "LABEL", "species", "contamination_flag"])
    reason = pd.Series("", index=m.index)
    reason[~m.species.isin(COMPLEX)] = "not in complex"
    reason[m.LABEL.notna() & (reason == "")] = "control strain"
    reason[(m.contamination_flag == "contaminated") & (reason == "")] = "contaminated"
    reason[m.NR.isin(EXTRA_EXCL) & (reason == "")] = "low quality (manual)"
    return pd.DataFrame({"sample": m.PROBENNUMMER, "NR": m.NR, "species": m.species,
                         "keep": reason == "", "reason": reason}).sort_values("NR").reset_index(drop=True)
