# kansasii-mst

Core-SNP minimum spanning trees of the M. kansasii complex (all complex, kansasii, pseudokansasii, persicum).
Only what is easy to get wrong:

- QC exclusions live in `scripts/common.py` (`load_samples`): contaminated flag, manual NR 10/21, controls, non-complex species.
- Core is computed **per subset** from `snps.aligned.fa`, not from the complex-wide snippy-core alignment.
- Run scripts via `sbatch run_all.sbatch` (env `env_mst`, see `environment.yml`); not on the login node.
- Type check: `~/data/conda/envs/env_typecheck/bin/mypy` from the repo root.
- Commit directly on `master`; push only when asked.
