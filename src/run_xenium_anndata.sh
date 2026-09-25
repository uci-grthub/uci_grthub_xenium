#!/bin/bash
#SBATCH --job-name=xenium_anndata
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/xenium_anndata.%j.out
#SBATCH --error=logs/xenium_anndata.%j.err

# Build one processed AnnData per Xenium sample plus its plots
# (rule xenium_anndata_all).
# Snakemake runs the per-sample jobs inside this allocation.
set -euo pipefail

source /opt/apps/mamba/24.3.0/etc/profile.d/conda.sh
conda activate snakemake_dfs

snakemake xenium_anndata_all --cores "${SLURM_CPUS_PER_TASK}" \
    --resources mem_mb=120000 --rerun-incomplete
