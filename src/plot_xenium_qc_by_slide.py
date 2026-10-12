#!/usr/bin/env python3
"""
Per-slide QC violin plots across all Xenium samples.

One page with one facet per per-cell statistic (genes per cell, transcripts
per cell, cell area); each facet shows one violin per slide, colored by
slide. Samples of the same slide are pooled into that slide's violin.

Usage:
  python src/plot_xenium_qc_by_slide.py --h5ads a.h5ad b.h5ad ... \
      --slide-map '{"sample_id": "slide_id", ...}' \
      --out output/xenium/anndata/plots/qc_violin_by_slide.png
"""
import argparse
import json
import os

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

QC_KEYS = ["n_genes_by_counts", "total_counts", "cell_area"]
QC_TITLES = {
    "n_genes_by_counts": "Genes per cell",
    "total_counts": "Transcripts per cell",
    "cell_area": "Cell area (µm²)",
}
# Okabe-Ito palette: distinct and colorblind-safe from the first two colors
# on (tab20's first two are both blue)
SLIDE_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#F0E442",
                "#56B4E9", "#E69F00", "#000000"]
# y-axis clipped at this percentile so outlier tails don't squash the violins
YMAX_PERCENTILE = 99.5


def load_qc_table(h5ads, slide_map):
    """Pool the QC columns of every sample into one table: slide, sample_id,
    one row per cell."""
    frames = []
    for path in h5ads:
        adata = ad.read_h5ad(path, backed="r")
        sample = str(adata.obs["sample_id"].iloc[0])
        if sample not in slide_map:
            raise ValueError(f"No slide in --slide-map for sample '{sample}' "
                             f"({path})")
        obs = adata.obs[["sample_id"] + QC_KEYS].copy()
        obs["slide"] = slide_map[sample]
        frames.append(obs)
        adata.file.close()
    return pd.concat(frames, ignore_index=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ads", required=True, nargs="+")
    parser.add_argument("--slide-map", required=True,
                        help="JSON object mapping sample_id -> slide_id")
    parser.add_argument("--out", required=True, help="output PNG path")
    args = parser.parse_args()

    slide_map = json.loads(args.slide_map)
    qc = load_qc_table(args.h5ads, slide_map)

    slides = sorted(qc["slide"].unique())
    n_cells = qc.groupby("slide").size()
    colors = dict(zip(slides,
                      (SLIDE_COLORS * (len(slides) // len(SLIDE_COLORS) + 1))
                      [:len(slides)]))

    fig, axes = plt.subplots(1, len(QC_KEYS),
                             figsize=(4.0 * len(QC_KEYS), 4.2))
    for ax, key in zip(axes, QC_KEYS):
        bodies = ax.violinplot([qc.loc[qc["slide"] == slide, key].to_numpy()
                                for slide in slides],
                               showmedians=True, showextrema=True)
        for body, slide in zip(bodies["bodies"], slides):
            body.set_facecolor(colors[slide])
            body.set_edgecolor("black")
            body.set_linewidth(0.5)
            body.set_alpha(0.9)
        for part in ("cbars", "cmins", "cmaxes", "cmedians"):
            bodies[part].set_color("black")
            bodies[part].set_linewidth(1.0)
        ax.set_ylim(0, qc[key].quantile(YMAX_PERCENTILE / 100) * 1.05)
        ax.set_xticks(range(1, len(slides) + 1))
        ax.set_xticklabels([f"{slide}\n({n_cells[slide]:,} cells)"
                            for slide in slides], fontsize=9)
        ax.set_title(QC_TITLES[key])
        ax.set_ylabel("")
        ax.grid(axis="y", color="#dddddd", linewidth=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"QC by slide ({len(slides)} slides, "
                 f"{len(qc):,} cells after filtering; "
                 f"y clipped at {YMAX_PERCENTILE}th percentile)")
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"QC by slide -> {args.out} ({len(slides)} slides)")


if __name__ == "__main__":
    main()
