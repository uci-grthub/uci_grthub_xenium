#!/usr/bin/env python3
"""
Azimuth cell type plots for a merged Xenium AnnData (from
merge_xenium_anndata.py), using the per-sample predictions of run_azimuth.R.

Merged cell ids are the per-sample ids (<sample>_<cell>), so predictions are
joined on them directly; cells below --min-score are "unassigned", as in
add_azimuth_labels.py. Only obs and the UMAP are read, not the expression
matrix.

Writes into --outdir:

  umap_azimuth.png         UMAP colored by Azimuth cell type (top) and by Leiden
  umap_azimuth_scores.png  UMAP colored by prediction and mapping score
                           (low mapping score: unlike any reference type,
                           e.g. tumor cells)
  azimuth_composition.png  fraction of each sample's cells per cell type
  azimuth_vs_leiden.png    fraction of each Leiden cluster per cell type

Usage:
  python src/plot_merged_azimuth.py --h5ad merged/<group>.h5ad \
      --predictions output/xenium/azimuth/*_predictions.csv \
      --labels cell_type --min-score 0.5 --outdir merged/plots/<group>
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

from add_azimuth_labels import SEQ_CMAP, UNASSIGNED, add_predictions, plot_vs_leiden
from plot_xenium_anndata import save


def read_obs_umap(path):
    """obs and X_umap of an .h5ad, without its expression matrix."""
    backed = ad.read_h5ad(path, backed="r")
    adata = ad.AnnData(obs=backed.obs.copy(),
                       obsm={"X_umap": np.asarray(backed.obsm["X_umap"])})
    backed.file.close()
    return adata


def read_predictions(paths, obs_names):
    """Predictions for the cells in obs_names, from per-sample CSVs."""
    frames = [pd.read_csv(p, index_col="cell_id") for p in paths]
    pred = pd.concat(frames)
    pred = pred[pred.index.isin(obs_names)]
    print(f"{len(pred)} of {len(obs_names)} merged cells have a prediction")
    return pred


def point_size(adata):
    return max(0.5, min(20.0, 120000 / adata.n_obs))


def plot_umap(adata, key, name, outdir):
    # unassigned cells first, so the labeled cells are drawn on top of them
    order = np.argsort((adata.obs[key] != UNASSIGNED).to_numpy(), kind="stable")
    view = adata[order]
    # stacked, so a long multi-column legend never runs into the other panel
    fig, axes = plt.subplots(2, 1, figsize=(14, 20), gridspec_kw={"hspace": 0.15})
    sc.pl.umap(view, color=key, ax=axes[0], show=False, size=point_size(adata),
               legend_fontsize=8, title=f"{name}: Azimuth {key.removeprefix('azimuth_')}")
    sc.pl.umap(view, color="leiden", ax=axes[1], show=False, size=point_size(adata),
               legend_loc="on data", legend_fontsize=8, title=f"{name}: Leiden")
    save(fig, outdir, "umap_azimuth.png")


def plot_umap_scores(adata, key, name, outdir):
    size = point_size(adata)
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for ax, col, title in zip(axes, [f"{key}_score", "azimuth_mapping_score"],
                              ["Prediction score", "Mapping score"]):
        sc.pl.umap(adata, color=col, ax=ax, show=False, size=size, cmap=SEQ_CMAP,
                   vmin=0, vmax=1, title=f"{name}: {title}")
    save(fig, outdir, "umap_azimuth_scores.png")


def plot_composition(adata, key, name, outdir):
    table = pd.crosstab(adata.obs["sample_id"], adata.obs[key], normalize="index")
    colors = dict(zip(adata.obs[key].cat.categories, adata.uns[f"{key}_colors"]))
    fig, ax = plt.subplots(figsize=(0.6 * len(table) + 8, 6))
    table.plot.bar(stacked=True, ax=ax, width=0.85,
                   color=[colors[c] for c in table.columns])
    ax.set_ylabel("fraction of cells")
    ax.set_xlabel("")
    ax.set_ylim(0, 1)
    ax.set_title(f"{name}: Azimuth {key.removeprefix('azimuth_')} per sample")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, frameon=False)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    save(fig, outdir, "azimuth_composition.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ad", required=True, help="merged .h5ad")
    parser.add_argument("--predictions", nargs="+", required=True,
                        help="per-sample CSVs from run_azimuth.R")
    parser.add_argument("--labels", nargs="+", default=["cell_type"],
                        help="reference label columns; plots use the first")
    parser.add_argument("--min-score", type=float, default=0.5,
                        help="prediction score below which a cell is unassigned")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    name = os.path.basename(args.h5ad).removesuffix(".h5ad")
    adata = read_obs_umap(args.h5ad)
    add_predictions(adata, read_predictions(args.predictions, adata.obs_names),
                    args.labels, args.min_score)

    key = f"azimuth_{args.labels[0]}"
    os.makedirs(args.outdir, exist_ok=True)
    plot_umap(adata, key, name, args.outdir)
    plot_umap_scores(adata, key, name, args.outdir)
    plot_composition(adata, key, name, args.outdir)
    plot_vs_leiden(adata, key, name, args.outdir)

    print(f"{name}: {adata.n_obs} cells, "
          f"{(adata.obs[key] == UNASSIGNED).mean():.0%} unassigned -> {args.outdir}")
    print(adata.obs[key].value_counts().to_string())


if __name__ == "__main__":
    main()
