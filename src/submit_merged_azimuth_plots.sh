#!/bin/bash
#SBATCH --job-name=merged_azimuth_plots
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/merged_azimuth_plots.%j.out
#SBATCH --error=logs/merged_azimuth_plots.%j.err

# Azimuth cell type UMAPs for merged objects (src/plot_merged_azimuth.py),
# from the per-sample predictions in output/xenium/azimuth/ (run
# src/run_azimuth_annotation.sh first). Plots go to <dir>/plots/<group>/.
# Usage, from the repository root:
#   sbatch src/submit_merged_azimuth_plots.sh <merged.h5ad> [merged.h5ad ...]
# e.g. sbatch src/submit_merged_azimuth_plots.sh output/xenium/anndata/merged_nocorrection/*.h5ad
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

if [ $# -eq 0 ]; then
    echo "usage: sbatch $0 <merged.h5ad> [merged.h5ad ...]" >&2
    exit 1
fi

# interpreter and threshold from config.yaml, as the rules use them
read -r PYTHON MIN_SCORE < <(python3 -c '
import yaml; c = yaml.safe_load(open("config.yaml"))
print(c["xenium"]["python"], (c.get("azimuth") or {}).get("min_prediction_score", 0.5))')

for h5ad in "$@"; do
    outdir="$(dirname "$h5ad")/plots/$(basename "$h5ad" .h5ad)"
    PYTHONNOUSERSITE=1 "$PYTHON" src/plot_merged_azimuth.py --h5ad "$h5ad" \
        --predictions output/xenium/azimuth/*_predictions.csv \
        --min-score "$MIN_SCORE" --outdir "$outdir"
done
