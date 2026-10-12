# Preparing Xenium data for NCBI GEO

Xenium is imaging-based: the raw data GEO asks for are the morphology images and
decoded transcripts, not FASTQ files. Stage the original Xenium output folders
(plus any processed analysis files) in one shared location for submission, and
include a checksum manifest so integrity can be verified after transfer.

## What GEO requires (10x Genomics Xenium In Situ)

Raw data files, per region (`output-XETG*` directory):

- High-resolution morphology images (`morphology.ome.tif`, `morphology_focus/`)
- Decoded transcripts (`transcripts.parquet` or `transcripts.csv.gz`) — the
  per-transcript gene calls with x/y/z coordinates and quality scores

Processed data files:

- Transcript counts, MEX (`cell_feature_matrix/`) or HDF5 (`cell_feature_matrix.h5`)
- Cell summary (`cells.parquet` or `cells.csv.gz`)
- Cell segmentation (`cell_boundaries.parquet` or `.csv.gz`)
- Nucleus segmentation (`nucleus_boundaries.parquet` or `.csv.gz`)
- Optional supplementary objects (e.g. `.h5ad` or Seurat `.rds`)

Metadata:

- `docs/Xenium_template.xlsx` — GEO's Xenium metadata template (study,
  samples, protocols). The `METADATA EXAMPLE` sheet shows a filled-in record.
- Panel and run information: each region's `experiment.xenium` and
  `gene_panel.json`.
- A sample-to-group mapping (sample id, slide, region, condition), and for
  regions cut from a larger section, the ROI coordinates.

Upload individual files over FTP; do not bundle them into `.zip`/`.tar`
archives.

## Staging checklist

1. Copy each region's raw and processed files under `geo_staging/<sample>/`
   (`geo_staging/` is gitignored).
2. Put shared metadata (`experiment.xenium` per region, sample mapping table)
   under `geo_staging/metadata/`.
3. Write a checksum manifest:
   `(cd geo_staging && find . -type f ! -name md5sums.txt -exec md5sum {} + > md5sums.txt)`.
4. Fill in `docs/Xenium_template.xlsx`, one row per sample, listing each data
   file in the data file columns.

Example Xenium GEO records: GSE243168, GSE294965, GSE276273.

Project-specific requests and notes go in `docs/local/` (gitignored).
