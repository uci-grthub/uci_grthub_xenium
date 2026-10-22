#!/bin/bash
#SBATCH --job-name=tangram
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/tangram.%j.out
#SBATCH --error=logs/tangram.%j.err

# Tangram cell type annotation of one merged Xenium object
# (src/run_tangram.py) with a reference from src/export_rds_counts.sh.
# RBC is dropped from the reference (no hemoglobin genes on the 5K panel) and
# cells with < 50 transcripts are labeled low_counts, not mapped.
# Predictions go to output/xenium/tangram/<group>_predictions.csv, plots to
# <dir>/plots/<group>/. 32G is enough for the skin slide (18k cells); give
# the tumor slide (862k cells) more, e.g. --mem=128G.
# Usage, from the repository root:
#   sbatch src/submit_tangram.sh <merged.h5ad> <reference dir> [label column] [tag]
# A tag (e.g. subtype, for data/reference/tangram_subtype) keeps the coarse
# results: predictions go to <group>_<tag>_predictions.csv, plots to
# plots/<group>/<tag>/.
# e.g.
#   sbatch src/submit_tangram.sh output/xenium/anndata/merged/slide_0120636.h5ad data/reference/tangram/Skin
#   sbatch --mem=128G src/submit_tangram.sh output/xenium/anndata/merged/slide_0120633.h5ad data/reference/tangram/Tumor
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

if [ $# -lt 2 ]; then
    echo "usage: sbatch $0 <merged.h5ad> <reference dir> [label column] [tag]" >&2
    exit 1
fi
h5ad=$1
ref=$2
label=${3:-cell_type_015}
tag=${4:-}
# built by src/install_tangram_env.sh
python=${TANGRAM_PYTHON:-envs/tangram/.pixi/envs/default/bin/python}
group=$(basename "$h5ad" .h5ad)

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
PYTHONNOUSERSITE=1 "$python" src/run_tangram.py --h5ad "$h5ad" --reference "$ref" \
    --label "$label" --drop-labels RBC --min-counts 50 --predictions "output/xenium/tangram/${group}${tag:+_$tag}_predictions.csv" \
    --outdir "$(dirname "$h5ad")/plots/$group${tag:+/$tag}"
