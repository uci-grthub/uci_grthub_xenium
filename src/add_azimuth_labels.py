#!/usr/bin/env python3
"""
Add Azimuth predictions (from run_azimuth.R) to a per-sample Xenium AnnData
and plot them (rule xenium_azimuth_annotate).

For each reference label column <label> (e.g. cell_type) obs gains
azimuth_<label>, azimuth_<label>_raw and azimuth_<label>_score; every label
shares azimuth_mapping_score, and obsm["X_azimuth_refumap"] holds the cells
projected onto the reference UMAP. azimuth_<label> is "unassigned" where the
prediction score is below --min-score (or Azimuth returned no prediction);
azimuth_<label>_raw keeps Azimuth's label for every cell. Low-count cells
with no clear match tend to land on an arbitrary immune type with a low
score, which the threshold removes.

Writes --out (the annotated .h5ad) and into --outdir, for the first label:

  azimuth_spatial.png    predicted label on the tissue
  azimuth_umap.png       predicted label on the sample UMAP and on the
                         reference UMAP
  azimuth_scores.png     prediction and mapping score distributions, and
                         mapping score on the tissue
  azimuth_vs_leiden.png  fraction of each Leiden cluster per predicted label

Low mapping scores flag cells unlike anything in the reference; with a
normal-skin reference that includes tumor cells.

Usage:
  python src/add_azimuth_labels.py --h5ad <sample>.h5ad \
      --predictions <sample>_predictions.csv --labels cell_type \
      --min-score 0.5 --out annotated/<sample>.h5ad --outdir plots/<sample>
"""
import argparse
import os

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from plot_xenium_anndata import save, spatial, spot_size

SEQ_CMAP = "viridis"
UNASSIGNED = "unassigned"


def add_predictions(adata, predictions, labels, min_score):
    pred = predictions.reindex(adata.obs_names)
    n_missing = int(pred["mapping.score"].isna().sum())
    if n_missing:
        print(f"WARNING: {n_missing} cells have no Azimuth prediction")
    for label in labels:
        raw = pred[f"predicted.{label}"]
        score = pred[f"predicted.{label}.score"]
        assigned = raw.where(score >= min_score, UNASSIGNED)
        adata.obs[f"azimuth_{label}_raw"] = pd.Categorical(raw)
        adata.obs[f"azimuth_{label}"] = pd.Categorical(
            assigned, categories=sorted(set(assigned) - {UNASSIGNED}) + [UNASSIGNED])
        adata.obs[f"azimuth_{label}_score"] = score.to_numpy()
        # gray for unassigned, scanpy's default palette for the labels
        n = len(adata.obs[f"azimuth_{label}"].cat.categories) - 1
        palette = sc.pl.palettes.default_20 if n <= 20 else sc.pl.palettes.default_102
        adata.uns[f"azimuth_{label}_colors"] = list(palette[:n]) + ["#d3d3d3"]
    adata.obs["azimuth_mapping_score"] = pred["mapping.score"].to_numpy()
    adata.obsm["X_azimuth_refumap"] = pred[["refUMAP_1", "refUMAP_2"]].to_numpy()


def plot_spatial(adata, key, sample, outdir):
    fig, ax = plt.subplots(figsize=(11, 9))
    spatial(adata, key, ax, title=f"{sample}: Azimuth {key.removeprefix('azimuth_')}")
    save(fig, outdir, "azimuth_spatial.png")


def plot_umaps(adata, key, sample, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    sc.pl.umap(adata, color=key, ax=axes[0], show=False, legend_loc="on data",
               legend_fontsize=7, title=f"{sample}: sample UMAP")
    sc.pl.embedding(adata, basis="X_azimuth_refumap", color=key, ax=axes[1], show=False,
                    size=spot_size(adata), title=f"{sample}: projected on reference UMAP")
    save(fig, outdir, "azimuth_umap.png")


def plot_scores(adata, key, sample, outdir):
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.5))
    for ax, col in zip(axes[:2], [f"{key}_score", "azimuth_mapping_score"]):
        ax.hist(adata.obs[col].dropna(), bins=50, color="#4c72b0")
        ax.set_xlabel(col)
        ax.set_ylabel("cells")
        ax.set_xlim(0, 1)
    spatial(adata, "azimuth_mapping_score", axes[2], cmap=SEQ_CMAP, vmin=0, vmax=1,
            title="Mapping score")
    fig.suptitle(sample)
    save(fig, outdir, "azimuth_scores.png")


def plot_vs_leiden(adata, key, sample, outdir):
    table = pd.crosstab(adata.obs["leiden"], adata.obs[key], normalize="index")
    fig, ax = plt.subplots(figsize=(0.45 * table.shape[1] + 3, 0.35 * table.shape[0] + 2))
    im = ax.imshow(table.to_numpy(), cmap=SEQ_CMAP, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(table.shape[1]), table.columns, rotation=90)
    ax.set_yticks(np.arange(table.shape[0]), table.index)
    ax.set_xlabel(f"Azimuth {key.removeprefix('azimuth_')}")
    ax.set_ylabel("Leiden cluster")
    ax.set_title(f"{sample}: fraction of each cluster")
    fig.colorbar(im, ax=ax, shrink=0.6)
    save(fig, outdir, "azimuth_vs_leiden.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--predictions", required=True, help="CSV from run_azimuth.R")
    parser.add_argument("--labels", nargs="+", default=["cell_type"],
                        help="reference label columns; plots use the first")
    parser.add_argument("--min-score", type=float, default=0.5,
                        help="prediction score below which a cell is unassigned")
    parser.add_argument("--out", required=True, help="annotated .h5ad")
    parser.add_argument("--outdir", required=True, help="plot directory")
    args = parser.parse_args()

    adata = ad.read_h5ad(args.h5ad)
    sample = adata.uns.get("sample_id", os.path.basename(args.h5ad).removesuffix(".h5ad"))
    predictions = pd.read_csv(args.predictions, index_col="cell_id")
    add_predictions(adata, predictions, args.labels, args.min_score)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    adata.write_h5ad(args.out, compression="gzip")

    key = f"azimuth_{args.labels[0]}"
    os.makedirs(args.outdir, exist_ok=True)
    plot_spatial(adata, key, sample, args.outdir)
    plot_umaps(adata, key, sample, args.outdir)
    plot_scores(adata, key, sample, args.outdir)
    plot_vs_leiden(adata, key, sample, args.outdir)

    print(f"{sample}: {adata.n_obs} cells, median mapping score "
          f"{np.nanmedian(adata.obs['azimuth_mapping_score']):.2f}, "
          f"{(adata.obs[key] == UNASSIGNED).mean():.0%} unassigned "
          f"(score < {args.min_score}) -> {args.out}")
    print(adata.obs[key].value_counts().to_string())


if __name__ == "__main__":
    main()
