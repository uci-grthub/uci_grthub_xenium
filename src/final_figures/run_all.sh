#!/usr/bin/env bash
# Run all final-figure generation scripts in order.
#
# Usage:
#   bash run_all.sh                    # Python/plotnine scripts (default)
#   bash run_all.sh --r                # R/ggplot2 scripts
#   bash run_all.sh --skip-extract     # skip h5ad extraction step
#   bash run_all.sh --r --skip-extract
#
# Python scripts: p01–p05_*.py  (require: conda activate squidpy)
# R scripts:    01–05_*.R  (require: module load R/4.4.2)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTDIR="${PROJECT_DIR}/src/final_figures/output/figures"

USE_R=false
SKIP_EXTRACT=false
for arg in "$@"; do
  [[ "$arg" == "--r"            ]] && USE_R=true
  [[ "$arg" == "--skip-extract" ]] && SKIP_EXTRACT=true
done

echo "=================================================="
echo " Final Figure Generation"
echo " Mode:    $([ "$USE_R" = true ] && echo 'R / ggplot2' || echo 'Python / plotnine')"
echo " Project: ${PROJECT_DIR}"
echo " Output:  ${OUTDIR}"
echo "=================================================="

# ── Environment ────────────────────────────────────────────────────────────────
source /opt/apps/mamba/24.3.0/etc/profile.d/conda.sh

conda activate squidpy 2>/dev/null || true
PYTHON=$(which python3 2>/dev/null || echo "")
if [[ -z "$PYTHON" ]]; then
  echo "ERROR: python3 not found in squidpy env."; exit 1
fi
echo "Python: $PYTHON"

if [[ "$USE_R" == true ]]; then
  module load R/4.4.2 2>/dev/null || true
  RSCRIPT=$(which Rscript 2>/dev/null || echo "")
  if [[ -z "$RSCRIPT" ]]; then
    echo "ERROR: Rscript not found. Run: module load R/4.4.2"; exit 1
  fi
  echo "Rscript: $RSCRIPT"
fi
echo ""

# ── Step 0: Extract data from h5ad (always Python / squidpy) ──────────────────
if [[ "$SKIP_EXTRACT" == false ]]; then
  echo "[1/6] Extracting h5ad data (UMAP + expression)..."
  $PYTHON "${SCRIPT_DIR}/00_extract_h5ad_data.py"
  echo "      Done."
else
  echo "[1/6] Skipping h5ad extraction (--skip-extract)."
fi

# ── Helper to run one script ───────────────────────────────────────────────────
run_step() {
  local step="$1" label="$2" r_script="$3" py_script="$4"
  echo ""
  echo "[${step}/6] ${label}..."
  if [[ "$USE_R" == true ]]; then
    $RSCRIPT "${SCRIPT_DIR}/${r_script}"
  else
    $PYTHON "${SCRIPT_DIR}/${py_script}"
  fi
  echo "      Done."
}

run_step 2 "Generating volcano plots"      01_volcano_plots.R     p01_volcano_plots.py
run_step 3 "Generating heatmaps"           02_heatmaps.R          p02_heatmaps.py
run_step 4 "Generating UMAP plots"         03_umap_plots.R        p03_umap_plots.py
run_step 5 "Running pathway enrichment"    04_pathway_enrichment.R p04_pathway_enrichment.py
run_step 6 "Generating network figures"    05_network_viz.R       p05_network_viz.py

# ── Summary ────────────────────────────────────────────────────────────────────
echo ""
echo "=================================================="
echo " All figures saved to:"
echo "   ${OUTDIR}"
echo ""
ls -lh "${OUTDIR}"/*.pdf 2>/dev/null || echo " (no PDFs yet)"
echo "=================================================="
