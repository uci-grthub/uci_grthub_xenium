#!/usr/bin/env python3
"""
Cell type UMAPs of a merged Xenium AnnData from label transfer predictions
(src/run_tangram.py, src/run_seurat_transfer.R).

The predictions CSV has one row per cell: cell_id, <method>_<label>,
<method>_<label>_score, total_counts and prob_<type>; cells that were not
transferred are labeled "low_counts". With --min-score, transferred cells
whose score is below it are "unassigned". Only obs and the UMAP of the .h5ad
are read, so the same predictions can be plotted on any merged object of the
same cells (e.g. merged/ and merged_nocorrection/).

Writes into --outdir:

  umap_<method>.png         UMAP colored by cell type (top) and by Leiden
  umap_<method>_score.png   UMAP colored by the label's score
  <method>_composition.png  fraction of each sample's cells per cell type
  <method>_vs_leiden.png    fraction of each Leiden cluster per cell type

Usage:
  python src/plot_label_transfer.py --h5ad merged/<group>.h5ad --method seurat \
      --predictions <csv> --min-score 0.5 --outdir merged/plots/<group>
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from plot_merged_azimuth import point_size, read_obs_umap
from plot_xenium_anndata import save

SEQ_CMAP = "viridis"
LOW_COUNTS = "low_counts"
UNASSIGNED = "unassigned"
GREYS = {LOW_COUNTS: "#d3d3d3", UNASSIGNED: "#8c8c8c"}
TITLES = {"tangram": "Tangram", "seurat": "Seurat transfer",
          "seurat_sct": "Seurat transfer (SCT)"}


def add_predictions(adata, pred, key, min_score):
    """Labels and scores of pred (indexed by cell id) as obs columns of adata."""
    missing = ~adata.obs_names.isin(pred.index)
    if missing.any():
        raise ValueError(f"{missing.sum()} cells have no prediction")
    pred = pred.reindex(adata.obs_names)
    labels = pred[key].astype(str)
    score = pred[f"{key}_score"]
    if min_score:
        labels = labels.where((labels == LOW_COUNTS) | (score >= min_score), UNASSIGNED)
    types = sorted(c.removeprefix("prob_") for c in pred.columns if c.startswith("prob_"))
    extra = [g for g in (UNASSIGNED, LOW_COUNTS) if (labels == g).any()]
    adata.obs[key] = pd.Categorical(labels, categories=types + extra)
    adata.obs[f"{key}_score"] = score.to_numpy()
    palette = sc.pl.palettes.default_20 if len(types) <= 20 else sc.pl.palettes.default_102
    colors = dict(zip(types, palette)) | GREYS
    adata.uns[f"{key}_colors"] = [colors[c] for c in adata.obs[key].cat.categories]


def plot_umap(adata, key, title, name, method, outdir):
    size = point_size(adata)
    # low_counts and unassigned cells first, so labeled cells are drawn on top
    order = np.argsort((~adata.obs[key].isin(list(GREYS))).to_numpy(), kind="stable")
    view = adata[order]
    fig, axes = plt.subplots(2, 1, figsize=(14, 20), gridspec_kw={"hspace": 0.15})
    sc.pl.umap(view, color=key, ax=axes[0], show=False, size=size,
               legend_fontsize=8, title=f"{name}: {title}")
    sc.pl.umap(view, color="leiden", ax=axes[1], show=False, size=size,
               legend_loc="on data", legend_fontsize=8, title=f"{name}: Leiden")
    save(fig, outdir, f"umap_{method}.png")


def plot_umap_score(adata, key, title, name, method, outdir):
    fig, ax = plt.subplots(figsize=(8, 7))
    sc.pl.umap(adata, color=f"{key}_score", ax=ax, show=False, size=point_size(adata),
               cmap=SEQ_CMAP, vmin=0, vmax=1, title=f"{name}: {title} label score")
    save(fig, outdir, f"umap_{method}_score.png")


def plot_composition(adata, key, title, name, method, outdir):
    table = pd.crosstab(adata.obs["sample_id"], adata.obs[key], normalize="index")
    colors = dict(zip(adata.obs[key].cat.categories, adata.uns[f"{key}_colors"]))
    fig, ax = plt.subplots(figsize=(0.6 * len(table) + 8, 6))
    table.plot.bar(stacked=True, ax=ax, width=0.85,
                   color=[colors[c] for c in table.columns])
    ax.set_ylabel("fraction of cells")
    ax.set_xlabel("")
    ax.set_ylim(0, 1)
    ax.set_title(f"{name}: {title} per sample")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, frameon=False)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    save(fig, outdir, f"{method}_composition.png")


def plot_vs_leiden(adata, key, title, name, method, outdir):
    table = pd.crosstab(adata.obs["leiden"], adata.obs[key], normalize="index")
    fig, ax = plt.subplots(figsize=(0.45 * table.shape[1] + 3, 0.35 * table.shape[0] + 2))
    im = ax.imshow(table.to_numpy(), cmap=SEQ_CMAP, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(table.shape[1]), table.columns, rotation=90)
    ax.set_yticks(np.arange(table.shape[0]), table.index)
    ax.set_xlabel(title)
    ax.set_ylabel("Leiden cluster")
    ax.set_title(f"{name}: fraction of each cluster")
    fig.colorbar(im, ax=ax, shrink=0.6)
    save(fig, outdir, f"{method}_vs_leiden.png")


def plot_all(adata, pred, method, label, min_score, name, outdir):
    key = f"{method}_{label}"
    add_predictions(adata, pred, key, min_score)
    title = f"{TITLES.get(method, method)} {label}"
    os.makedirs(outdir, exist_ok=True)
    for plot in (plot_umap, plot_umap_score, plot_composition, plot_vs_leiden):
        plot(adata, key, title, name, method, outdir)
    print(f"{name}: {adata.n_obs} cells -> {outdir}")
    print(adata.obs[key].value_counts().to_string())


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ad", required=True, help="merged .h5ad")
    parser.add_argument("--predictions", required=True, help="predictions CSV")
    parser.add_argument("--method", required=True, help="column prefix: tangram, seurat, seurat_sct")
    parser.add_argument("--label", default="cell_type_015", help="reference label column")
    parser.add_argument("--min-score", type=float, default=0.0,
                        help="score below which a transferred cell is unassigned (0: none)")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    name = os.path.basename(args.h5ad).removesuffix(".h5ad")
    pred = pd.read_csv(args.predictions, index_col="cell_id", low_memory=False)
    plot_all(read_obs_umap(args.h5ad), pred, args.method, args.label, args.min_score,
             name, args.outdir)


if __name__ == "__main__":
    main()
