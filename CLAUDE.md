# kansasii-mst

> **OBSOLETE (2026-10-30).** Superseded by Ridom SeqSphere+ (cgMLST allele MST) as the source of truth; not used for the final report. Gaps (SNP distances, masking) and the export list: `kansasii-docs/docs/RIDOM.md`. Do not extend or rerun for the report; outputs are kept for comparison.

Core-SNP minimum spanning trees of the M. kansasii complex (all complex, kansasii, pseudokansasii, persicum).
Only what is easy to get wrong:

- QC exclusions live in `scripts/common.py` (`load_samples`): `excluded` (contaminated: CheckM > 5 %, or < 90 % of the cgMLST loci called; set in screening_map_results.csv by immensekansasii 06, the same column in every repo), controls, non-complex species. No manual per-isolate exclusions.
- Core is computed **per subset** from `snps.aligned.fa`, not from the complex-wide snippy-core alignment.
- Masking: `00_build_mask_bed.py` builds `output/kansasii-mst/results/mask/mask.bed` (PE/PPE, IS, phage by GFF product regex + minimap2 self-repeats); snippy has no `--mask`, so it is applied in `02` (`--no-mask` to disable). Rerun 00-03 after changing it.
- Node merging threshold `MERGE_SNPS` (common.py) was tuned on the MST edge-weight distribution; recheck if the mask or QC set changes.
- Slurm `--wrap` runs under sh (no `source`); use `submit_all.sbatch` or call `~/data/conda/envs/env_mst/bin/python` directly.
- Run scripts via `sbatch submit_all.sbatch` (env `env_mst`, see `environment.yml`); not on the login node.
- Type check: `~/data/conda/envs/env_typecheck/bin/mypy` from the repo root.
- Commit directly on `master`; push only when asked.
