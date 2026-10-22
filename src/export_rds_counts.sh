#!/bin/bash
#SBATCH --job-name=rds_counts
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=logs/rds_counts.%j.out
#SBATCH --error=logs/rds_counts.%j.err

# Tangram references from a labeled Seurat .rds (src/export_rds_counts.R):
# raw counts per tissue, subsampled per cell type, Xenium panel genes only.
# 96G fits the 12.6 GB GSE279468 whole-skin object (peak ~50 GB to load).
# Usage, from the repository root:
#   sbatch src/export_rds_counts.sh <file.rds> <outdir> <label_column> <group_column> \
#       <xenium features.tsv.gz> [max_per_label] [min_per_label] [overrides.csv]
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

rscript=$(python3 -c 'import yaml; print(yaml.safe_load(open("config.yaml"))["azimuth"]["rscript"])')
"$rscript" src/export_rds_counts.R "$@"
