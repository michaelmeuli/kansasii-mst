#!/usr/bin/env bash
# Core-SNP minimum spanning tree figures: scripts/00 .. 03 in numeric order (00 builds the PE/PPE+IS+phage+repeat mask BED).
#   sbatch submit_all.sbatch
# Inputs : data/imm/screening_map_results.csv, runs/mkan329/assembly/results/*/5_typing/kansasii_snippy
# Outputs: output/mst/
#SBATCH --job-name=kansasii-mst
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/mst/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/kansasii-mst
OUTPUT=$OUTPUT_ROOT/mst
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mst
cd "$REPO"
for s in scripts/[0-9][0-9]_*.py; do
  echo "=== $s ($(date +%T))"
  python "$s"
done
echo "done $(date +%T)"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT" "$REPO" submit_all.sbatch "$SCREENING_MAP_RESULTS" "$ASSEMBLIES" "$ATCC12478_FASTA"
