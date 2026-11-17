#!/bin/bash
#SBATCH --job-name=seurat_obs
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/seurat_obs.%j.out
#SBATCH --error=logs/seurat_obs.%j.err

# Seurat (SCT) cell types into obs of every per-sample, merged/ and
# merged_nocorrection/ AnnData (rule xenium_seurat_obs_all,
# src/add_transfer_obs.py): seurat_celltype and seurat_celltype_score.
# Only obs is rewritten and mtimes are kept; --rerun-triggers mtime stops
# code/params provenance from rebuilding the AnnData objects first.
# Usage, from the repository root:
#   sbatch src/submit_seurat_obs.sh
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

snakemake xenium_seurat_obs_all --cores 1 --rerun-triggers mtime
