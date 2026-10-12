#!/bin/bash
#SBATCH --job-name=xenium_report
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --error=logs/xenium_report.%j.err
#SBATCH --output=logs/xenium_report.%j.out

# Build output/Xenium_Project_Report.pdf (rule xenium_report) and anything it
# still needs (merged objects and their plots for every merge group).
# Extra snakemake options go in SNAKEMAKE_ARGS. Submit from the repository root.
cd "${SLURM_SUBMIT_DIR:-.}"
set -euo pipefail

source /opt/apps/mamba/24.3.0/etc/profile.d/conda.sh
conda activate snakemake_dfs

snakemake --cores "${SLURM_CPUS_PER_TASK}" --resources mem_mb=60000 \
    --rerun-triggers mtime --rerun-incomplete ${SNAKEMAKE_ARGS:-} -- xenium_report
