#!/bin/bash
#SBATCH --job-name=seurat_transfer
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=logs/seurat_transfer.%j.out
#SBATCH --error=logs/seurat_transfer.%j.err

# Seurat label transfer of one merged Xenium slide object
# (src/run_seurat_transfer.R) with a reference from src/export_rds_counts.sh,
# then cell type plots (src/plot_label_transfer.py, score < 0.5 unassigned)
# for merged/ and merged_nocorrection/ (same cells, different UMAP).
# Normalization SCT (default; method seurat_sct) or LogNormalize (method
# seurat). Predictions go to output/xenium/seurat_transfer/
# <group>_predictions_sct.csv (SCT) or <group>_predictions.csv (LogNormalize),
# plots (<method> in the file names) to
# output/xenium/anndata/{merged,merged_nocorrection}/plots/<group>/.
# The skin slide peaked at 3.9 GB with LogNormalize; SCT query residuals for
# the 563k-cell yum17 sample add ~9 GB, so 96G for the tumor slide.
# Usage, from the repository root:
#   sbatch src/submit_seurat_transfer.sh <group> <reference dir> [SCT|LogNormalize]
# e.g.
#   sbatch src/submit_seurat_transfer.sh slide_0120636 data/reference/tangram/Skin
#   sbatch src/submit_seurat_transfer.sh slide_0120633 data/reference/tangram/Tumor
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    echo "usage: sbatch $0 <group> <reference dir> [SCT|LogNormalize]" >&2
    exit 1
fi
group=$1
ref=$2
normalization=${3:-SCT}
case $normalization in
    SCT) method=seurat_sct; suffix=_sct ;;
    LogNormalize) method=seurat; suffix= ;;
    *) echo "normalization must be SCT or LogNormalize" >&2; exit 1 ;;
esac
label=cell_type_015
min_score=0.5

# Rscript (Seurat 5) and the squidpy Python from config.yaml, as the rules use
# them; the plotting Python is built by src/install_tangram_env.sh
read -r RSCRIPT XENIUM_PYTHON < <(python3 -c '
import yaml; c = yaml.safe_load(open("config.yaml"))
print(c["azimuth"]["rscript"], c["xenium"]["python"])')
PLOT_PYTHON=${PLOT_PYTHON:-envs/tangram/.pixi/envs/default/bin/python}

anndata=output/xenium/anndata
outdir=output/xenium/seurat_transfer
predictions=$outdir/${group}_predictions${suffix}.csv
query=$outdir/query_mtx/${group}_${method}
trap 'rm -rf "$query"' EXIT

PYTHONNOUSERSITE=1 "$XENIUM_PYTHON" src/h5ad_to_mtx.py --h5ad "$anndata/merged/$group.h5ad" \
    --layer counts --obs-columns sample_id --outdir "$query"
"$RSCRIPT" src/run_seurat_transfer.R --query-mtx "$query" --reference "$ref" \
    --out "$predictions" --label "$label" --drop-labels RBC --min-counts 50 \
    --normalization "$normalization"

for dir in merged merged_nocorrection; do
    PYTHONNOUSERSITE=1 "$PLOT_PYTHON" src/plot_label_transfer.py --method "$method" \
        --h5ad "$anndata/$dir/$group.h5ad" --predictions "$predictions" \
        --label "$label" --min-score "$min_score" --outdir "$anndata/$dir/plots/$group"
done
