# kansasii-mst

Minimum spanning tree (MST) figures from core-genome SNPs, for four sample sets:
the whole *M. kansasii* complex, *M. kansasii* sensu stricto, *M. pseudokansasii* and *M. persicum*.

```bash
sbatch submit_all.sbatch   # scripts/00 .. 03, about 2 minutes (outputs: output/kansasii-mst/results/)
```

| Script | Does |
|---|---|
| `00_build_mask_bed.py` | mask BED: PE/PPE, IS, phage genes (reference GFF product names) + self-repeats -> `output/kansasii-mst/results/mask/mask.bed` |
| `01_select_samples.py` | QC exclusion, subset membership -> `output/kansasii-mst/results/samples_used.tsv` |
| `02_core_snp_matrices.py` | per-subset core SNPs and Hamming distance matrices |
| `03_mst_figures.py` | MST (scipy) drawn with networkx, labels repelled with adjustText (`--subset`, `--merge N`) -> `output/kansasii-mst/results/mst_<subset>.{png,pdf}` |

**Environment:** `environment.yml` (`env_mst`: pandas, numpy, scipy, matplotlib, networkx, adjusttext). Script 00 also needs
`bedtools` and `minimap2`, taken from the `snippy` conda env (`SNIPPY_BIN` in `scripts/common.py`). Type check:
`~/data/conda/envs/env_typecheck/bin/mypy`. Run heavy steps with `sbatch`, not on the login node.

**Excluded:** genomes with `excluded` = true in `screening_map_results.csv` (CheckM contamination above 5 %, < 90 % of the cgMLST loci called, or a non-complex clinical sample, NR 184; the same column in every repo), control
strains (`LABEL` set) and non-complex species.

**Core SNPs:** from the snippy alignments to NC_022663.1; core = ACGT in every genome of the subset, so the complex-wide core
is much smaller than the within-species cores and the complex MST resolves less. No recombination filtering.
Identical genomes are collapsed into one node (area = genome count); edge labels are SNP differences.
Nodes are labelled with the sample number only (`Mkan329-001` -> `1`; up to 8 per node, then `+k more`). Node area depends only on the number of genomes in the node (never on edges) and grows as its square root, so big nodes stay moderate; a size legend shows the scale. Nodes of 4+ genomes have a framed list box joined to the node by a leader line. `output/kansasii-mst/results/mst_nodes_<subset>.tsv` lists every sample per node (full ids and numbers).

**Node merging:** MST edges of at most N SNPs are contracted (single linkage), so crowded centres become fewer, larger
nodes (area = genome count, label shows `(n=...)`). Defaults in `scripts/common.py` (`MERGE_SNPS`): complex 1000 (coarse on purpose: the complex figure is for the overview; 172 genomes -> 75 nodes, persicum collapses to a few nodes), kansasii 50, kansasii 50,
pseudokansasii 20, persicum 100, chosen where the many short MST edges end (e.g. kansasii 117 haplotypes -> 74 nodes at 50).
Override with `python scripts/03_mst_figures.py --merge N` (`--merge 0` = identical genomes only). `mst_nodes_<subset>.tsv`
has `max_internal_snps`, the largest pairwise distance inside each node (single linkage can chain).

**Masking:** snippy has no `--mask` (only `snippy-core` does), and the per-sample `snps.aligned.fa` are unmasked, so the BED is
applied when the core is built (`02_core_snp_matrices.py`, same semantics as `snippy-core --mask`; `--no-mask` turns it off).
Masked on the reference GCF_000157895.3: PE/PPE genes 254 CDS (384 kb), IS/transposase 70 (56 kb), phage/integrase 40 (38 kb;
"phage shock" and "phage resistance" genes are not prophage and excluded) and self-repeats from a minimap2 self-alignment, the
`snippy-core --mask auto` recipe (111 regions, 351 kb); merged 375 regions, 695 kb = 10.6% of the reference (673 kb on the
chromosome). Effect on variable core sites: complex 260,020 -> 254,643, kansasii 222,680 -> 210,259, pseudokansasii 373 -> 338,
persicum 31,630 -> 29,826 (`output/kansasii-mst/results/core_snp_summary.tsv`). Node-merge thresholds were re-checked and kept.
Phage regions are annotation-based (no PHASTER/PhiSpy available), so unannotated prophage is not masked.

## Layout

```
scripts/          00_build_mask_bed.py ... 03_mst_figures.py, common.py; kansasii_paths.py = shared constants (vendored)
submit_all.sbatch runs 00-03 in order (log in output/kansasii-mst/logs/)
```

Paths come from `scripts/kansasii_paths.py` (same names in every repo: `ASSEMBLIES`, `ATCC12478_FASTA`, `SCREENING_MAP_RESULTS`, `OUTPUT_ROOT`) and are
flags of the scripts: `--assemblies`, `--screening-map-results` (01), `--ref-fasta`, `--ref-gff` (00). Edit
`kansasii-docs/templates/kansasii_paths.py`, not the copy (`bash repos/check_sync.sh --fix`).

## Provenance

The sbatch files write `MANIFEST.txt` (commit, clean/dirty flag, script, conda env hash, inputs) into the output
directory when they finish; see `kansasii-docs/docs/PROVENANCE.md`. Project index: `kansasii-docs/README.md`.
