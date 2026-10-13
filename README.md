# Xenium workflow

Snakemake workflow for 10x Xenium data on HPC3 (SLURM). The Xenium
instrument runs its onboard analysis, so the delivery already holds cells,
transcripts and the cell-feature matrix per region; this workflow discovers
those region outputs and gathers their QC into one table. The repository also
carries the lab's Cell Ranger and downstream single-cell rules, which Xenium
projects do not use.

Project-specific notes (run ids, slides, panel, results) live in `docs/local/`,
which is gitignored.

## Layout
- `Snakefile` — the workflow. All rules live here except the includes below.
- `config.example.yaml` — template for `config.yaml` (paths, per-rule
  parameters, merge groups, mail address, SLURM accounts). `config.yaml` itself
  is per-project and gitignored.
- `rules/` — `scenic.smk` and `azimuth.smk` are included; `jbrowse.smk`, `numbat.smk` and
  `pileup_and_phase.smk` are not wired into the Snakefile.
- `src/` — R, Python and shell scripts called by the rules (not `scripts/`).
- `data/` — the Xenium delivery, mirrored by `src/rsync_xenium.sh`:
  `data/<order>/<run>/output-XETG<serial>__<slide>__<region>__<datetime>/`.
- `metadata/` — sample metadata.
- `envs/` — conda environment YAMLs.
- `output/` — results; `output/logs/` for per-rule logs.

## Getting the data

`src/rsync_xenium.sh` copies a delivery into `data/` and verifies every file by
checksum (Xenium deliveries ship no `md5sums.txt`):

```bash
mkdir -p logs/rsync_xenium
sbatch src/rsync_xenium.sh /path/to/delivery/<order> data/<order>
```

## Running

Create the project config first; the Snakefile stops if `config.yaml` is
missing:

```bash
cp config.example.yaml config.yaml   # then set paths, email, slurm_account
```

Snakemake lives in a conda environment, not a module:

```bash
source /opt/apps/mamba/24.3.0/etc/profile.d/conda.sh
conda activate snakemake_dfs
snakemake -n          # prints the discovered Xenium samples
snakemake --cores 2   # summary tables; local, take seconds
sbatch src/run_xenium_anndata.sh   # per-sample AnnData (rule xenium_anndata_all)
```

## Xenium

### How samples are discovered

There is no sample list to maintain. Every `experiment.xenium` found under
`xenium.data_root` (recursively) marks one region output directory. The sample
id is `<slide_id>_<region_name>`, both read from that manifest, e.g.
`0012345_tumor_1`; the slide prefix keeps ids unique if a region name
is reused on another slide. A duplicate id stops the workflow at parse time.

```yaml
xenium:
  data_root: "data"
```

### Outputs

Under `output/xenium/`:

- `analysis_summary_table.csv` / `.html` — every metric table shown in the
  per-region `analysis_summary.html` (run info, region info, panel, software,
  decoding yield, segmentation, stain segmentation, negative controls, debug
  info) plus the raised alerts; one row per metric, one column per sample.
  Values keep the report's formatting (`94.2%`, `44,187`). Built by
  `src/aggregate_xenium_summaries.py`, which reads the JSON the report embeds.
- `metrics_summary.csv` — every region's `metrics_summary.csv` stacked, one row
  per sample, unformatted numbers; use this one for plotting or thresholds.
- `analysis_summaries/<sample>_analysis_summary.html` — the per-region reports
  copied side by side.

The script also runs on its own:

```bash
python src/aggregate_xenium_summaries.py \
    --summaries data/*/*/output-*/analysis_summary.html \
    --csv summary.csv --html summary.html
```

### Per-sample AnnData

`snakemake xenium_anndata_all` (submitted by `src/run_xenium_anndata.sh`, which
runs every sample inside one 16-core allocation) writes
`output/xenium/anndata/<sample>.h5ad` with `src/xenium_to_anndata.py`:

1. `cell_feature_matrix.h5`, gene expression features only (control probes
   and codewords dropped); cell metadata from `cells.csv.gz`, which reads
   reliably where some deliveries ship truncated `cells.parquet` files.
2. Keep cells with at least 10 genes (`xenium.min_genes`) and genes seen in at
   least 5 cells (`xenium.min_cells`).
3. Normalize by cell size: counts / `cell_area` x 100 µm²
   (`xenium.target_area`, a fixed constant like Seurat's `scale.factor`), then
   `log1p`.
4. PCA (30 PCs), neighbors, UMAP, Leiden (`xenium.leiden_resolution`, 0.3).

Raw counts are in `layers["counts"]`, centroids in `obsm["spatial"]`, the
per-cell scale in `obs["size_factor"]`, and Xenium's own total (genes plus
controls) in `obs["xenium_total_counts"]`. Because genes are filtered per
sample, `var` differs between samples; concatenate with `join="outer"` or
re-filter from `layers["counts"]`.

### Merged AnnData

The same target also writes merged objects under `output/xenium/anndata/merged/` with
`src/merge_xenium_anndata.py`. Every sample is rebuilt from
`layers["counts"]` (per-sample normalization is discarded), genes are
outer-joined, and all cells are normalized with the same fixed
`xenium.target_area` (100 µm²) as the per-sample objects. Any other merging
code should call `cell_size_normalize()` from `src/xenium_to_anndata.py`
with that constant, never a per-sample median.

PCA/UMAP/Leiden run on the merged object; `obs["sample_id"]` is the batch and
per-sample clusters are kept as `obs["leiden_sample"]`. Batch correction is
set by `xenium.batch_correction`: `none` (default) or `harmony`, which runs
Harmony on the PCs over `sample_id` (`obsm["X_pca_harmony"]`) and builds
neighbors/UMAP/Leiden from those; `X` and `X_pca` stay uncorrected. The choice
is recorded in `uns["preprocessing"]["batch_correction"]` and in the
`umap.png` title. Changing it applies to every merge group.

Which samples are merged is set by `xenium.merge_groups` (group name -> glob
over sample ids); each group writes `merged/<group>.h5ad` and plots in
`merged/plots/<group>/` (UMAP by sample and cluster, one UMAP panel per
sample to judge batch effects, cluster composition per sample, marker
dotplot):

```yaml
merge_groups:
  merged: "*"                      # all samples
  slide_0012345: "0012345_*"       # every region on slide 0012345
```

Build one group with
`snakemake output/xenium/anndata/merged/plots/<group>/markers.csv`, or on
SLURM `sbatch src/submit_merge_plots.sh <group> [group ...]`.

### Cell type annotation (Azimuth)

`snakemake xenium_azimuth_all` (submitted by `src/run_azimuth_annotation.sh`)
labels every per-sample object with Azimuth (`rules/azimuth.smk`). Azimuth
has no mouse skin reference, so one is built from a labeled CELLxGENE atlas,
by default Tabula Muris Senis 10x (`azimuth.reference_url`):

1. `azimuth_reference_download`: the atlas `.h5ad` to
   `data/reference/azimuth/reference.h5ad`.
2. `azimuth_reference_export` (`src/h5ad_to_mtx.py`): raw counts of the cell
   types in `azimuth.keep_labels` (types that occur in skin and skin tumors;
   organ-specific types dropped), at most `azimuth.max_cells_per_label` each,
   panel genes only (matched by Ensembl id from `gene_panel.json`, renamed to
   the panel's symbols).
3. `azimuth_reference_build` (`src/build_azimuth_reference.R`): SCTransform
   v2, PCA, UMAP and `AzimuthReference` -> `output/xenium/azimuth/reference/`
   (`ref.Rds`, `idx.annoy`, `reference_umap.png`).
4. `xenium_azimuth` (`src/run_azimuth.R`): `RunAzimuth` on each sample's
   `layers["counts"]` -> `output/xenium/azimuth/<sample>_predictions.csv`.
5. `xenium_azimuth_annotate` (`src/add_azimuth_labels.py`):
   `output/xenium/anndata/annotated/<sample>.h5ad` with
   `obs["azimuth_cell_type"]` (thresholded), `obs["azimuth_cell_type_raw"]`,
   `obs["azimuth_cell_type_score"]`,
   `obs["azimuth_mapping_score"]` and `obsm["X_azimuth_refumap"]`, and plots
   in `plots/<sample>/` (`azimuth_spatial.png`, `azimuth_umap.png`,
   `azimuth_scores.png`, `azimuth_vs_leiden.png`).

Cells with a prediction score below `azimuth.min_prediction_score` (0.5)
are labeled `unassigned` in `obs["azimuth_cell_type"]`; Azimuth's label for
every cell stays in `obs["azimuth_cell_type_raw"]`. Without the threshold,
low-count cells land on immune types they do not express.

The reference is normal tissue: it has no melanocytes or tumor cells, so
tumor cells get the nearest reference label with a low mapping score. Check
`azimuth_mapping_score` as well before trusting a label in tumor samples.

The R steps use the pixi env in `envs/azimuth/` (Seurat 5, Azimuth 0.5.0;
locked in `pixi.lock`, GitHub packages pinned to commits), built once with
`sbatch src/install_azimuth_env.sh`.

### Project report

`snakemake xenium_report` (or `sbatch src/submit_report.sh`) writes
`output/Xenium_Project_Report.pdf`: project and run information, summary
statistics, panel and software versions, the pipeline with its parameters,
per-sample Xenium QC and filtering (cells before/after QC, median transcripts
and genes, Q20 decoded, transcripts assigned, clusters), every merged object
with its UMAP, and the output files. Set the author in `config.yaml`
(`report.author`).

Two rules build it: `xenium_report_stats` reads every per-sample and merged
`.h5ad` (backed, in the `xenium.python` env) into
`output/xenium/report_stats.json` with `src/xenium_report_stats.py`;
`xenium_report` renders the PDF with `src/generate_xenium_report.py`
(reportlab; interpreter set by `report.python`, by default the pixi env in
`envs/report/`, built once with `pixi install --manifest-path envs/report/pixi.toml`). The report describes each merged
object's batch correction as it was built and flags objects that no longer
match `xenium.batch_correction`.

Other downstream analysis of the region outputs is done in notebooks: see `src/read_in_xenium.ipynb`, `src/samui_xenium_process.py` and
`src/subset_xenium_by_coordinates_and_combine.py`.

## Requirements

- Snakemake 9 with the `cluster-generic` executor plugin, and pandas (conda env
  `snakemake_dfs`).
- scanpy, leidenalg and igraph for the AnnData rule (`xenium.python`, default
  the `squidpy_dfs` env).
- R 4.4 with Seurat 5 and Azimuth for the cell type annotation rules
  (`azimuth.rscript`, the pixi env in `envs/azimuth/`).
- CellRanger 8.0.1 (module) for the single-cell rules only.

## Other rules in this repository

Single-cell rules, working but idle in this project (no Cell Ranger samples
configured; target `cellranger_all`): `fastqc`, `multiqc`, `collect_web_summaries`, `cellranger_count`, `multi_sample_summary`,
`tenx_harmony_integration`, `tenx_scvi_integration`, `tenx_harmony_notebook`,
`plot_filtering_timeline`, `cellbender`, `scrublet`.

Legacy, and broken as written — do not expect these to run without repair:

- `scvelo` and `pileup_and_phasing` call `src/compute_velocity.py` and
  `src/pileup_and_phase.R`, neither of which is in the repository.
- `velocyto` passes `config["..."]` as literal shell text rather than an
  evaluated value.
- `rules/scenic.smk` consumes `output/loom/<sample>.loom`, which was produced by
  a `loompy` rule that has been removed.
- `rules/jbrowse.smk`, `rules/numbat.smk` and `rules/pileup_and_phase.smk` are
  not included by the Snakefile.

## Documentation
- `rules/README.md` — downstream analysis rules.
- `docs/prep_geo_data.md` — GEO submission staging.
- `docs/Xenium_template.xlsx` — GEO's blank Xenium metadata template.
- `docs/local/` — project-specific notes and results (gitignored).
