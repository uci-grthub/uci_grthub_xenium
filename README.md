# GaneA Xenium (order 0926X-02)

Snakemake workflow for this project's Xenium data on HPC3 (SLURM). The Xenium
instrument runs its onboard analysis, so the delivery already holds cells,
transcripts and the cell-feature matrix per region; this workflow discovers
those region outputs and gathers their QC into one table. The repository also
carries the lab's Cell Ranger and downstream single-cell rules, which this
project does not use.

Run: `20260916_0926X-02_GaneA`, two slides (`0120633`, `0120636`), 15 regions,
Xenium Mouse 5K Pan Tissue & Pathways panel, FFPE.

## Layout
- `Snakefile` — the workflow. All rules live here except the includes below.
- `config.yaml` — paths and per-rule parameters.
- `rules/` — `scenic.smk` is included; `jbrowse.smk`, `numbat.smk` and
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
sbatch src/rsync_xenium.sh /path/to/delivery/0926X-02_GaneA data/0926X-02_GaneA
```

## Running

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
`0120633_b16_10724_tumor`; the slide prefix keeps ids unique if a region name
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
   and codewords dropped); cell metadata from `cells.csv.gz`. The
   `cells.parquet` files in this delivery fail to read ("Unexpected end of
   stream") for 12 of 15 regions, so the CSV is used.
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

PCA/UMAP/Leiden run on the merged object with no batch correction
(`obs["sample_id"]` is the batch); per-sample clusters are kept as
`obs["leiden_sample"]`.

Which samples are merged is set by `xenium.merge_groups` (group name -> glob
over sample ids); each group writes `merged/<group>.h5ad` and plots in
`merged/plots/<group>/` (UMAP by sample and cluster, one UMAP panel per
sample to judge batch effects, cluster composition per sample, marker
dotplot):

```yaml
merge_groups:
  merged: "*"                      # all 15 samples
  slide_0120636: "0120636_*"       # the 8 regions on slide 0120636
```

Build one group with
`snakemake output/xenium/anndata/merged/plots/slide_0120636/markers.csv`.

Other downstream analysis of the region outputs is done in notebooks: see `src/read_in_xenium.ipynb`, `src/samui_xenium_process.py` and
`src/subset_xenium_by_coordinates_and_combine.py`.

## Requirements

- Snakemake 9 with the `cluster-generic` executor plugin, and pandas (conda env
  `snakemake_dfs`).
- scanpy, leidenalg and igraph for the AnnData rule (`xenium.python`, default
  the `squidpy_dfs` env).
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
