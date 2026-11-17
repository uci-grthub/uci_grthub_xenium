#!/usr/bin/env python3
"""
Cell types of one Xenium section on the tissue, one panel per annotation
method (rule xenium_celltype_methods_spatial).

Reads the per-sample AnnData (spatial coordinates only) and:

  --azimuth   per-sample CSV from src/run_azimuth.R (predicted.<label>,
              predicted.<label>.score)
  --transfer  <method>=<csv> pairs from src/run_tangram.py or
              src/run_seurat_transfer.R (<method>_<label>,
              <method>_<label>_score, prob_<type>); these are per merged
              group, so only the section's cells are used

Cells whose score is below --min-score are "unassigned"; cells the transfer
methods skipped (< 50 transcripts) stay "low_counts". Transfer methods share
one palette (same reference types); Azimuth has its own.

Long thin sections (the skin strips) get one wide row per method, others one
column per method. Writes --out.

Usage:
  python src/plot_celltype_methods_spatial.py --h5ad anndata/<sample>.h5ad \
      --azimuth azimuth/<sample>_predictions.csv \
      --transfer tangram=tangram/<group>_predictions.csv \
                 seurat_sct=seurat_transfer/<group>_predictions_sct.csv \
      --label cell_type_015 --min-score 0.5 --out plots/<sample>/celltype_methods_spatial.png
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

from plot_xenium_anndata import save, style_spatial

LOW_COUNTS = "low_counts"
UNASSIGNED = "unassigned"
GREYS = {UNASSIGNED: "#8c8c8c", LOW_COUNTS: "#d3d3d3"}
TITLES = {"azimuth": "Azimuth", "tangram": "Tangram", "seurat": "Seurat transfer",
          "seurat_sct": "Seurat transfer (SCT)"}
# sections whose height / width is below this are drawn as stacked rows
STRIP_ASPECT = 0.4


def read_spatial(path):
    """obs names and spatial coordinates of an .h5ad, without its matrix."""
    backed = ad.read_h5ad(path, backed="r")
    adata = ad.AnnData(obs=pd.DataFrame(index=backed.obs_names.copy()),
                       obsm={"spatial": np.asarray(backed.obsm["spatial"])})
    sample = backed.uns.get("sample_id", os.path.basename(path).removesuffix(".h5ad"))
    backed.file.close()
    return adata, sample


def palette(types):
    colors = sc.pl.palettes.default_20 if len(types) <= 20 else sc.pl.palettes.default_102
    return dict(zip(types, colors)) | GREYS


def assign(labels, score, min_score):
    """Labels with transferred cells below min_score set to unassigned."""
    labels = labels.fillna(UNASSIGNED).astype(str)
    low = (labels != LOW_COUNTS) & ~(score >= min_score)
    return labels.where(~low, UNASSIGNED)


def read_azimuth(path, label, obs_names, min_score):
    key = f"predicted.{label}"
    pred = pd.read_csv(path, index_col="cell_id", usecols=["cell_id", key, f"{key}.score"])
    pred = pred.reindex(obs_names)
    return assign(pred[key], pred[f"{key}.score"], min_score)


def read_transfer(path, method, label, obs_names, min_score):
    """Labels for obs_names and the reference types (from the prob_ columns)."""
    key = f"{method}_{label}"
    header = pd.read_csv(path, nrows=0).columns
    types = sorted(c.removeprefix("prob_") for c in header if c.startswith("prob_"))
    pred = pd.read_csv(path, index_col="cell_id", usecols=["cell_id", key, f"{key}_score"])
    missing = ~obs_names.isin(pred.index)
    if missing.any():
        raise ValueError(f"{path}: {missing.sum()} of {len(obs_names)} cells have no prediction")
    pred = pred.reindex(obs_names)
    return assign(pred[key], pred[f"{key}_score"], min_score), types


def add_labels(adata, key, labels, colors):
    types = [t for t in colors if t not in GREYS and (labels == t).any()]
    extra = sorted(set(labels) - set(types) - set(GREYS))  # labels missing from the palette
    cats = types + extra + [g for g in GREYS if (labels == g).any()]
    colors = colors | palette(extra)
    adata.obs[key] = pd.Categorical(labels.to_numpy(), categories=cats)
    adata.uns[f"{key}_colors"] = [colors[c] for c in cats]


def point_size(adata, width, height):
    """Marker area (pt^2) filling about the panel's area per cell."""
    span = np.ptp(adata.obsm["spatial"], axis=0)
    # data extent inside a width x height inch panel, in points
    scale = min(width / span[0], height / span[1]) * 72
    area = span[0] * span[1] * scale ** 2
    return float(np.clip(0.6 * area / adata.n_obs, 0.2, 10.0))


def draw(ax, adata, key, title, size):
    # unassigned / low_counts first so labeled cells are drawn on top
    order = np.argsort((~adata.obs[key].isin(list(GREYS))).to_numpy(), kind="stable")
    view = adata[order]
    labels = adata.obs[key]
    n_lab = int((~labels.isin(list(GREYS))).sum())
    sc.pl.embedding(view, basis="spatial", color=key, ax=ax, show=False, size=size,
                    frameon=True, legend_fontsize=7,
                    title=f"{title}: {n_lab:,} of {len(labels):,} cells labeled "
                          f"({100 * n_lab / len(labels):.0f}%)")
    style_spatial(ax)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ad", required=True, help="per-sample AnnData")
    parser.add_argument("--azimuth", help="Azimuth predictions CSV of the sample")
    parser.add_argument("--azimuth-label", default="cell_type")
    parser.add_argument("--transfer", nargs="*", default=[],
                        help="<method>=<predictions csv> of the sample's merge group")
    parser.add_argument("--label", default="cell_type_015", help="transfer label column")
    parser.add_argument("--min-score", type=float, default=0.5)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    adata, sample = read_spatial(args.h5ad)
    panels = []
    if args.azimuth:
        labels = read_azimuth(args.azimuth, args.azimuth_label, adata.obs_names, args.min_score)
        types = sorted(set(labels) - set(GREYS))
        add_labels(adata, "azimuth", labels, palette(types))
        panels.append(("azimuth", f"Azimuth {args.azimuth_label}"))
    transfer = dict(item.split("=", 1) for item in args.transfer)
    shared = None
    for method, path in transfer.items():
        labels, types = read_transfer(path, method, args.label, adata.obs_names, args.min_score)
        shared = shared or palette(types)
        add_labels(adata, method, labels, shared)
        panels.append((method, f"{TITLES.get(method, method)} {args.label}"))
    if not panels:
        parser.error("give --azimuth and/or --transfer")

    span = np.ptp(adata.obsm["spatial"], axis=0)
    aspect = span[1] / span[0]
    n = len(panels)
    if aspect < STRIP_ASPECT:
        width = 18
        row = max(2.5, width * 0.75 * aspect + 1)
        size = point_size(adata, width * 0.75, row - 1)
        fig, axes = plt.subplots(n, 1, figsize=(width, row * n), layout="constrained",
                                 squeeze=False)
    else:
        panel = 8
        height = min(max(panel * aspect, 3.0), 12.0)
        size = point_size(adata, panel, height)
        fig, axes = plt.subplots(1, n, figsize=(panel * 1.6 * n, height + 1),
                                 layout="constrained", squeeze=False)
    for ax, (key, title) in zip(axes.ravel(), panels):
        draw(ax, adata, key, title, size)
    fig.suptitle(f"{sample}: cell types by method (score < {args.min_score:g} unassigned)",
                 fontsize=14)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    save(fig, os.path.dirname(os.path.abspath(args.out)), os.path.basename(args.out))

    for key, _ in panels:
        counts = adata.obs[key].value_counts()
        print(f"{sample} {key}: " + ", ".join(f"{k} {v}" for k, v in counts.items() if v))
    print(f"-> {args.out}", flush=True)


if __name__ == "__main__":
    main()
