#!/usr/bin/env python3
"""
Xenium sections on the tissue before and after filtering cells by the number
of genes detected (sc.pp.filter_cells(min_genes=...)).

Reads the raw Xenium output (cell_feature_matrix.h5 + cells.csv.gz), not the
per-sample h5ad, which is already filtered. The sample id is
<slide_id>_<region_name> from the directory's experiment.xenium, as in the
Snakefile.

Writes <plots-root>/<sample>/gene_filter_spatial.png per section:

  per stage   kept cells colored by genes per cell (raw, then each cutoff)
  per cutoff  all cells, the ones the cutoff removes in red
  histogram   genes per cell of the raw cells with the cutoffs

Long thin sections get one wide row per panel, others a grid with one column
per stage. With --summary, also writes a CSV of cells kept per section and
cutoff.

Usage:
  python src/plot_gene_filter_spatial.py --xenium-dir <xenium output dir> [...] \
      --min-genes 10 50 --plots-root output/xenium/anndata/plots \
      --summary output/xenium/qc/gene_filter_counts.csv
"""
import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from plot_xenium_anndata import SEQ_CMAP, save, spot_size, style_spatial
from xenium_to_anndata import load_xenium

KEY = "n_genes_by_counts"
# sections whose height / width is below this are drawn as stacked rows
STRIP_ASPECT = 0.4


def sample_id(xenium_dir):
    with open(os.path.join(xenium_dir, "experiment.xenium")) as fh:
        meta = json.load(fh)
    region = meta.get("region_name") or os.path.basename(xenium_dir)
    slide = meta.get("slide_id") or ""
    return f"{slide}_{region}" if slide else region


def filter_stages(raw, cutoffs):
    stages = [("raw", raw)]
    for cutoff in cutoffs:
        kept = raw.copy()
        sc.pp.filter_cells(kept, min_genes=cutoff)
        stages.append((f">= {cutoff} genes", kept))
    return stages


def draw_kept(ax, name, adata, n_raw, size, vmax):
    sc.pl.embedding(adata, basis="spatial", color=KEY, ax=ax, show=False, size=size,
                    frameon=True, cmap=SEQ_CMAP, vmin=0, vmax=vmax,
                    title=f"{name}: {adata.n_obs:,} cells "
                          f"({100 * adata.n_obs / n_raw:.1f}%), genes per cell")
    style_spatial(ax)


def draw_removed(ax, name, adata, raw, size):
    xy = raw.obsm["spatial"]
    removed = ~raw.obs_names.isin(adata.obs_names)
    ax.scatter(*xy[~removed].T, s=size, c="#bbbbbb", rasterized=True, linewidths=0)
    ax.scatter(*xy[removed].T, s=size, c="#d62728", rasterized=True, linewidths=0)
    ax.set_title(f"{name}: removed in red ({removed.sum():,}, "
                 f"{100 * removed.mean():.1f}%)")
    style_spatial(ax)


def draw_histogram(ax, raw, cutoffs):
    genes = raw.obs[KEY].to_numpy()
    top = np.percentile(genes, 99.5)
    ax.hist(genes, bins=np.arange(0, top + 5, max(1, round(top / 120))), color="#8da0cb")
    for cutoff in cutoffs:
        ax.axvline(cutoff, color="#d62728", linestyle="--")
        ax.text(cutoff, ax.get_ylim()[1] * 0.95, f" {cutoff}", color="#d62728", va="top")
    ax.set_xlabel("genes per cell")
    ax.set_ylabel("cells")
    ax.set_title(f"raw: genes per cell, median {np.median(genes):.0f}")


def plot_section(raw, stages, cutoffs, sample, outdir):
    n_raw = raw.n_obs
    vmax = float(np.percentile(raw.obs[KEY], 99))
    span = np.ptp(raw.obsm["spatial"], axis=0)
    aspect = span[1] / span[0]
    size = min(spot_size(raw), 8.0)
    n_stages = len(stages)

    if aspect < STRIP_ASPECT:
        # one wide row per panel
        width = 16
        row = max(2.5, width * 0.8 * aspect + 1)
        n_rows = 2 * n_stages - 1
        fig, axes = plt.subplots(n_rows + 1, 1, figsize=(width, row * n_rows + 3.5),
                                 layout="constrained",
                                 gridspec_kw={"height_ratios": [row] * n_rows + [3.5]})
        kept_axes, removed_axes, hist_ax = axes[:n_stages], axes[n_stages:-1], axes[-1]
    else:
        # one column per stage; histogram under the raw column
        panel = 6.5
        height = min(max(panel * aspect, 3.0), 12.0) + 1
        fig, axes = plt.subplots(2, n_stages, figsize=(panel * 1.2 * n_stages, 2 * height),
                                 layout="constrained", squeeze=False)
        kept_axes, removed_axes, hist_ax = axes[0], axes[1, 1:], axes[1, 0]

    for ax, (name, adata) in zip(kept_axes, stages):
        draw_kept(ax, name, adata, n_raw, size, vmax)
    for ax, (name, adata) in zip(removed_axes, stages[1:]):
        draw_removed(ax, name, adata, raw, size)
    draw_histogram(hist_ax, raw, cutoffs)

    fig.suptitle(f"{sample}: filtering cells by genes detected")
    os.makedirs(outdir, exist_ok=True)
    save(fig, outdir, "gene_filter_spatial.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--xenium-dir", required=True, nargs="+",
                        help="Xenium output directories (output-XETG*)")
    parser.add_argument("--min-genes", type=int, nargs="+", default=[10, 50])
    parser.add_argument("--plots-root", default="output/xenium/anndata/plots",
                        help="plots go to <plots-root>/<sample>/")
    parser.add_argument("--summary", help="CSV of cells kept per section and cutoff")
    args = parser.parse_args()
    cutoffs = sorted(args.min_genes)

    rows = []
    for xenium_dir in args.xenium_dir:
        sample = sample_id(xenium_dir)
        raw = load_xenium(xenium_dir, sample)
        sc.pp.calculate_qc_metrics(raw, percent_top=None, inplace=True)
        stages = filter_stages(raw, cutoffs)
        outdir = os.path.join(args.plots_root, sample)
        plot_section(raw, stages, cutoffs, sample, outdir)

        row = {"sample_id": sample, "cells_raw": raw.n_obs,
               "median_genes": float(np.median(raw.obs[KEY]))}
        for cutoff, (_, adata) in zip(cutoffs, stages[1:]):
            row[f"cells_min{cutoff}"] = adata.n_obs
            row[f"pct_min{cutoff}"] = round(100 * adata.n_obs / raw.n_obs, 1)
        rows.append(row)
        print(" ".join(f"{k}={v}" for k, v in row.items()), flush=True)
        print(f"  -> {os.path.join(outdir, 'gene_filter_spatial.png')}", flush=True)

    if args.summary:
        os.makedirs(os.path.dirname(args.summary) or ".", exist_ok=True)
        table = pd.DataFrame(rows)
        table.to_csv(args.summary, index=False)
        print(table.to_string(index=False))


if __name__ == "__main__":
    main()
