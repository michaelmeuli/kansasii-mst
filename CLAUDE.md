# kansasii-mst

Core-SNP minimum spanning trees of the M. kansasii complex (all complex, kansasii, pseudokansasii, persicum).
Only what is easy to get wrong:

- QC exclusions live in `scripts/common.py` (`load_samples`): contaminated flag, manual NR 10/21, controls, non-complex species.
- Core is computed **per subset** from `snps.aligned.fa`, not from the complex-wide snippy-core alignment.
- Masking: `00_build_mask_bed.py` builds `output/mst/mask/mask.bed` (PE/PPE, IS, phage by GFF product regex + minimap2 self-repeats); snippy has no `--mask`, so it is applied in `02` (`--no-mask` to disable). Rerun 00-03 after changing it.
- Node merging threshold `MERGE_SNPS` (common.py) was tuned on the MST edge-weight distribution; recheck if the mask or QC set changes.
- Slurm `--wrap` runs under sh (no `source`); use `run_all.sbatch` or call `~/data/conda/envs/env_mst/bin/python` directly.
- Run scripts via `sbatch run_all.sbatch` (env `env_mst`, see `environment.yml`); not on the login node.
- Type check: `~/data/conda/envs/env_typecheck/bin/mypy` from the repo root.
- Commit directly on `master`; push only when asked.
