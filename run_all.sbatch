#!/usr/bin/env bash
#SBATCH --job-name=kansasii-mst
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.log
# Core-SNP minimum spanning tree figures: scripts/00 .. 03 (00 builds the PE/PPE+IS+phage+repeat mask BED) in numeric order.
#
#   sbatch run_all.sbatch
#
# Inputs  : data/imm/screening_map_results.csv, runs/mkan329/assembly/results/*/5_typing/kansasii_snippy
# Outputs : $KANSASII_ROOT/output/mst/  (KANSASII_ROOT defaults to the cluster root)
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mst
cd /shares/sander.imm.uzh/MM/kansasii/repos/kansasii-mst
for s in scripts/[0-9][0-9]_*.py; do
  echo "=== $s ($(date +%T))"
  python "$s"
done
echo "done $(date +%T)"
