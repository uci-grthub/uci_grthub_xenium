#!/bin/bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

ADULT_SRC="${REPO_DIR}/data/output-XETG00221__0069979__Adult_cohort__20251003__181017"
AGED_SRC="${REPO_DIR}/data/output-XETG00221__0069982__Aged_cohort__20251003__181018"
PROCESSED_SRC="${REPO_DIR}/output/scanpy"

STAGING="${REPO_DIR}/geo_staging"

echo "=== GEO staging ==="
echo "Destination: ${STAGING}"
echo "Date: $(date)"

# ---------------------------------------------------------------------------
# 1. Create directory layout
# ---------------------------------------------------------------------------
mkdir -p \
    "${STAGING}/Adult_cohort" \
    "${STAGING}/Aged_cohort" \
    "${STAGING}/processed" \
    "${STAGING}/metadata"

# ---------------------------------------------------------------------------
# 2. Copy raw Xenium output files (both cohorts)
# ---------------------------------------------------------------------------
copy_cohort() {
    local src="$1"
    local dst="$2"
    local label="$3"

    echo "--- Copying ${label} raw files ---"

    # Morphology images
    cp -v "${src}/morphology.ome.tif" "${dst}/"
    cp -rv "${src}/morphology_focus" "${dst}/"

    # Transcript coordinates (raw equivalent)
    cp -v "${src}/transcripts.parquet" "${dst}/"

    # Cell-feature matrix (MTX + HDF5)
    cp -rv "${src}/cell_feature_matrix" "${dst}/"
    cp -v "${src}/cell_feature_matrix.h5" "${dst}/"

    # Segmentation files
    cp -v "${src}/cell_boundaries.csv.gz" "${dst}/"
    cp -v "${src}/nucleus_boundaries.csv.gz" "${dst}/"
    cp -v "${src}/cells.csv.gz" "${dst}/"
    cp -v "${src}/coordinates.csv" "${dst}/"

    # Panel and run metadata
    cp -v "${src}/gene_panel.json" "${dst}/"
    cp -v "${src}/experiment.xenium" "${dst}/"
    cp -v "${src}/metadata.csv" "${dst}/"
}

copy_cohort "${ADULT_SRC}" "${STAGING}/Adult_cohort" "Adult"
copy_cohort "${AGED_SRC}"  "${STAGING}/Aged_cohort"  "Aged"

# ---------------------------------------------------------------------------
# 3. Copy processed analysis files
# ---------------------------------------------------------------------------
echo "--- Copying processed h5ad files ---"

for area in area_1 area_2 area_3 area_4; do
    cp -v "${PROCESSED_SRC}/${area}_harmony_integrated_processed.h5ad" \
          "${STAGING}/processed/"
done
cp -v "${PROCESSED_SRC}/all_harmony_integrated_processed.h5ad" \
      "${STAGING}/processed/"

# ---------------------------------------------------------------------------
# 4. Consolidated sample mapping
# ---------------------------------------------------------------------------
echo "--- Building combined sample map ---"

# Write header once, then append both cohorts (skipping their headers)
head -1 "${ADULT_SRC}/metadata.csv" > "${STAGING}/metadata/sample_map.csv"
tail -n +2 "${ADULT_SRC}/metadata.csv" >> "${STAGING}/metadata/sample_map.csv"
tail -n +2 "${AGED_SRC}/metadata.csv"  >> "${STAGING}/metadata/sample_map.csv"

# Copies of experiment manifests for provenance
cp -v "${ADULT_SRC}/experiment.xenium" "${STAGING}/metadata/Adult_cohort_experiment.xenium"
cp -v "${AGED_SRC}/experiment.xenium"  "${STAGING}/metadata/Aged_cohort_experiment.xenium"

# ---------------------------------------------------------------------------
# 5. Checksum manifest
# ---------------------------------------------------------------------------
echo "--- Generating MD5 checksums ---"
find "${STAGING}" -type f ! -name "MANIFEST.md5" \
    | sort \
    | xargs md5sum > "${STAGING}/MANIFEST.md5"

echo ""
echo "=== Done ==="
du -sh "${STAGING}"
echo "File count: $(wc -l < "${STAGING}/MANIFEST.md5")"
