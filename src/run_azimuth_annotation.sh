#!/bin/bash
#SBATCH --job-name=xenium_azimuth
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --time=24:00:00
#SBATCH --output=logs/xenium_azimuth.%j.out
#SBATCH --error=logs/xenium_azimuth.%j.err

# Azimuth cell type labels for every Xenium sample (rule xenium_azimuth_all):
# downloads and builds the reference once, then maps each sample.
# Snakemake runs the per-sample jobs inside this allocation. Pass a target to
# build only part of it, e.g. one sample:
#   sbatch src/run_azimuth_annotation.sh output/xenium/anndata/annotated/<sample>.h5ad
# Needs the R env from src/install_azimuth_env.sh.
set -euo pipefail

source /opt/apps/mamba/24.3.0/etc/profile.d/conda.sh
conda activate snakemake_dfs

# mtime only: existing per-sample AnnData objects are inputs, not rebuilt
# because the Snakefile changed
snakemake "${@:-xenium_azimuth_all}" --cores "${SLURM_CPUS_PER_TASK}" \
    --resources mem_mb=170000 --rerun-incomplete --rerun-triggers mtime
