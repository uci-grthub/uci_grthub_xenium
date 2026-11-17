#!/bin/bash
#SBATCH --job-name=gene_filter_spatial
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/gene_filter_spatial.%j.out
#SBATCH --error=logs/gene_filter_spatial.%j.err

# Tissue before/after filtering cells at 10 and 50 genes for every Xenium
# section (src/plot_gene_filter_spatial.py), from the raw Xenium output.
# Plots: output/xenium/anndata/plots/<sample>/gene_filter_spatial.png;
# counts: output/xenium/qc/gene_filter_counts.csv.
# Usage, from the repository root:
#   sbatch src/submit_gene_filter_spatial.sh [xenium output dir ...]
# (default: every experiment.xenium directory under data/)
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

PYTHON=/dfs9/ucightf-lab/kstachel/.conda/envs/squidpy_dfs/bin/python

if [ "$#" -gt 0 ]; then
    DIRS=("$@")
else
    mapfile -t DIRS < <(find data -name experiment.xenium -printf '%h\n' | sort)
fi
echo "${#DIRS[@]} sections"

PYTHONNOUSERSITE=1 "$PYTHON" src/plot_gene_filter_spatial.py \
    --xenium-dir "${DIRS[@]}" --min-genes 10 50 \
    --plots-root output/xenium/anndata/plots \
    --summary output/xenium/qc/gene_filter_counts.csv
