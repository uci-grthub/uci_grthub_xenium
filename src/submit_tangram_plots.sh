#!/bin/bash
#SBATCH --job-name=tangram_plots
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/tangram_plots.%j.out
#SBATCH --error=logs/tangram_plots.%j.err

# Tangram cell type UMAPs for merged objects from existing predictions
# (src/plot_label_transfer.py --method tangram), e.g. for merged_nocorrection/ after
# src/submit_tangram.sh ran on merged/ (same cells, different UMAP).
# Reads output/xenium/tangram/<group>_predictions.csv; plots go to
# <dir>/plots/<group>/.
# Usage, from the repository root:
#   sbatch src/submit_tangram_plots.sh <merged.h5ad> [merged.h5ad ...]
# e.g. sbatch src/submit_tangram_plots.sh output/xenium/anndata/merged_nocorrection/*.h5ad
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

if [ $# -eq 0 ]; then
    echo "usage: sbatch $0 <merged.h5ad> [merged.h5ad ...]" >&2
    exit 1
fi
# built by src/install_tangram_env.sh
python=${TANGRAM_PYTHON:-envs/tangram/.pixi/envs/default/bin/python}

for h5ad in "$@"; do
    group=$(basename "$h5ad" .h5ad)
    PYTHONNOUSERSITE=1 "$python" src/plot_label_transfer.py --method tangram --h5ad "$h5ad" \
        --predictions "output/xenium/tangram/${group}_predictions.csv" \
        --outdir "$(dirname "$h5ad")/plots/$group"
done
