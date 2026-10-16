# kansasii-mst

Minimum spanning tree (MST) figures from core-genome SNPs, for four sample sets:
the whole *M. kansasii* complex, *M. kansasii* sensu stricto, *M. pseudokansasii* and *M. persicum*.

```bash
sbatch run_all.sbatch      # scripts/01 .. 03, about minutes
```

| Script | Does |
|---|---|
| `01_select_samples.py` | QC exclusion, subset membership -> `output/mst/samples_used.tsv` |
| `02_core_snp_matrices.py` | per-subset core SNPs and Hamming distance matrices |
| `03_mst_figures.py` | MST (scipy) drawn with networkx -> `output/mst/mst_<subset>.{png,pdf}` |

**Excluded:** genomes flagged `contaminated` in `screening_map_results.csv`, NR 10 and 21 (as in kansasii-gwas), control
strains (`LABEL` set) and non-complex species.

**Core SNPs:** from the snippy alignments to NC_022663.1; core = ACGT in every genome of the subset, so the complex-wide core
is much smaller than the within-species cores and the complex MST resolves less. No recombination filtering.
Identical genomes are collapsed into one node (area = genome count); edge labels are SNP differences.
Nodes are labelled with `PROBENNUMMER` (up to 3 per node); `output/mst/mst_nodes_<subset>.tsv` lists every sample per node.

**Node merging:** MST edges of at most N SNPs are contracted (single linkage), so crowded centres become fewer, larger
nodes (area = genome count, label shows `(n=...)`). Defaults in `scripts/common.py` (`MERGE_SNPS`): complex 50, kansasii 50,
pseudokansasii 20, persicum 100, chosen where the many short MST edges end (e.g. kansasii 117 haplotypes -> 74 nodes at 50).
Override with `python scripts/03_mst_figures.py --merge N` (`--merge 0` = identical genomes only). `mst_nodes_<subset>.tsv`
has `max_internal_snps`, the largest pairwise distance inside each node (single linkage can chain).
