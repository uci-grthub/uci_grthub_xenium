#!/bin/bash
#SBATCH --job-name=rds_metadata
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/rds_metadata.%A_%a.out
#SBATCH --error=logs/rds_metadata.%A_%a.err

# Cell metadata (cell type labels etc.) of Seurat .rds files, one array task
# per file (src/extract_rds_metadata.R):
#   sbatch --array=1-$(wc -l < <list>) src/extract_rds_metadata.sh <list> <outdir>
# <list> holds one .rds path per line. 32G covers .rds files up to ~4 GB;
# give larger ones more, e.g. --mem=96G --time=4:00:00 for a 12.6 GB file.
set -euo pipefail
list=$1
outdir=$2
rds=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$list")
rscript=$(grep -A3 '^azimuth:' config.yaml | sed -n 's/^ *rscript: *"\?\([^"]*\)"\?/\1/p')
echo "task ${SLURM_ARRAY_TASK_ID}: $rds"
"$rscript" src/extract_rds_metadata.R "$rds" "$outdir"
