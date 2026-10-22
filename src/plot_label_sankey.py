#!/usr/bin/env python3
"""
Sankey (alluvial) diagram of cell type labels across label transfer runs of
the same cells (src/run_tangram.py, src/run_seurat_transfer.R), e.g. Tangram
and Seurat with the coarse and the tumor subtype reference.

Each --predictions NAME=CSV is one column, in the order given. Labels are
thresholded as in src/plot_label_transfer.py: transferred cells whose score is
below --min-score are "unassigned". Cells labeled "low_counts" were not
transferred, are the same in every run and are left out. Labels below
--min-frac of the cells in a column are lumped into "other".

Writes --out (PNG) and, beside it, one crosstab CSV of cell counts per pair
of adjacent columns: <out stem>_<a>_vs_<b>.csv.

Usage:
  python src/plot_label_sankey.py --out output/xenium/label_transfer_sankey/<group>_sankey.png \
      --predictions "Seurat SCT prior=output/xenium/seurat_transfer/<group>_predictions_sct.csv" \
      --predictions "Seurat SCT subtype=output/xenium/seurat_transfer/<group>_subtype_predictions_sct.csv" \
      --predictions "Tangram subtype=output/xenium/tangram/<group>_subtype_predictions.csv" \
      --predictions "Tangram prior=output/xenium/tangram/<group>_predictions.csv" \
      --min-score 0.5 --title <group>
"""
import argparse
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LOW_COUNTS = "low_counts"
UNASSIGNED = "unassigned"
OTHER = "other"
GREYS = {OTHER: "#bdbdbd", UNASSIGNED: "#8c8c8c"}
# melanocytic lineage first and in warm colors, so tumor flows stand out
TUMOR = ["Tumor", "common_tumor", "proliferating_tumor", "BcaPfl_enriched_tumor",
         "BcaPwt_enriched_tumor", "NCSC", "Melanocyte", "Melanocytes",
         "Keratinocyte/Melanocyte", "Schwann_cell"]
TUMOR_COLORS = ["#d62728", "#e6550d", "#fd8d3c", "#a50f15", "#fb6a4a", "#9e9ac8",
                "#8c564b", "#8c564b", "#c49c94", "#e377c2"]
NODE_WIDTH = 0.08
PAD = 0.006
SHOW_FRAC = 0.015


def read_labels(path, label, min_score):
    """Thresholded labels of one predictions CSV, indexed by cell id."""
    pred = pd.read_csv(path, index_col="cell_id")
    key = next(c for c in pred.columns if c.endswith(f"_{label}"))
    labels = pred[key].astype(str)
    if min_score:
        labels = labels.where((labels == LOW_COUNTS) | (pred[f"{key}_score"] >= min_score),
                              UNASSIGNED)
    return labels


def lump(labels, min_frac):
    frac = labels.value_counts(normalize=True)
    small = frac.index[(frac < min_frac) & ~frac.index.isin([UNASSIGNED])]
    return labels.where(~labels.isin(small), OTHER)


def label_order(df):
    totals = pd.concat([df[c] for c in df]).value_counts()
    rest = [l for l in totals.index if l not in TUMOR and l not in GREYS]
    return [l for l in TUMOR if l in totals.index] + rest + [l for l in GREYS if l in totals.index]


def label_colors(order):
    rest = [l for l in order if l not in TUMOR and l not in GREYS]
    # tab20b without its reds, tab20c without its oranges and greys: no clash with tumor/greys
    b, c = plt.get_cmap("tab20b").colors, plt.get_cmap("tab20c").colors
    groups = [b[0:4], b[4:8], b[8:12], b[16:20], c[0:4], c[8:12], c[12:16]]
    cool = [g[i] for i in (0, 2) for g in groups]
    return dict(zip(TUMOR, TUMOR_COLORS)) | dict(zip(rest, cool)) | GREYS


def node_positions(counts, order, n):
    """Top y of each label's node in a column, stacked top to bottom."""
    pos, y = {}, 1.0
    for l in order:
        if counts.get(l, 0):
            pos[l] = y
            y -= counts[l] / n + PAD
    return pos


def band(ax, x0, x1, y0, y1, h, color):
    """Sigmoid band of height h from top y0 at x0 to top y1 at x1."""
    x = np.linspace(x0, x1, 50)
    s = 1 / (1 + np.exp(-12 * ((x - x0) / (x1 - x0) - 0.5)))
    top = y0 + (y1 - y0) * s
    ax.fill_between(x, top - h, top, color=color, alpha=0.45, lw=0)


def plot(df, order, colors, title, out):
    n = len(df)
    names = list(df.columns)
    counts = {c: df[c].value_counts() for c in names}
    pos = {c: node_positions(counts[c], order, n) for c in names}
    xs = np.arange(len(names)) * 1.0
    fig, ax = plt.subplots(figsize=(4 * len(names), 10))

    for i, (a, b) in enumerate(zip(names, names[1:])):
        tab = pd.crosstab(df[a], df[b])
        out_y = dict(pos[a])
        in_y = dict(pos[b])
        for la in order:
            for lb in order:
                if la not in tab.index or lb not in tab.columns or not tab.at[la, lb]:
                    continue
                h = tab.at[la, lb] / n
                band(ax, xs[i] + NODE_WIDTH / 2, xs[i + 1] - NODE_WIDTH / 2,
                     out_y[la], in_y[lb], h, colors[la])
                out_y[la] -= h
                in_y[lb] -= h

    for x, c in zip(xs, names):
        for l, top in pos[c].items():
            h = counts[c][l] / n
            ax.add_patch(plt.Rectangle((x - NODE_WIDTH / 2, top - h), NODE_WIDTH, h,
                                       color=colors[l], lw=0))
            if h >= SHOW_FRAC:
                last = x == xs[-1]
                ax.text(x + (-1 if last else 1) * (NODE_WIDTH / 2 + 0.02), top - h / 2,
                        f"{l} ({100 * h:.1f}%)", va="center", fontsize=8,
                        ha="right" if last else "left")
        ax.text(x, 1.02, c, ha="center", va="bottom", fontsize=11, weight="bold")

    ax.set_xlim(xs[0] - 0.3, xs[-1] + 0.3)
    bottom = min(top - counts[c][l] / n for c in names for l, top in pos[c].items())
    ax.set_ylim(bottom - 0.01, 1.06)
    ax.axis("off")
    ax.set_title(title, fontsize=13)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)


def slug(name):
    return re.sub(r"\W+", "_", name).strip("_").lower()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--predictions", action="append", required=True,
                        help="NAME=CSV, one per column, in order")
    parser.add_argument("--label", default="cell_type_015", help="reference label column")
    parser.add_argument("--min-score", type=float, default=0.5,
                        help="score below which a transferred cell is unassigned (0: none)")
    parser.add_argument("--min-frac", type=float, default=0.005,
                        help="labels below this fraction of a column are lumped into 'other'")
    parser.add_argument("--title", default="")
    parser.add_argument("--out", required=True, help="output PNG")
    args = parser.parse_args()

    runs = dict(p.split("=", 1) for p in args.predictions)
    df = pd.DataFrame({name: read_labels(path, args.label, args.min_score)
                       for name, path in runs.items()})
    if df.isna().any().any():
        raise ValueError("predictions CSVs are not of the same cells")
    low = (df == LOW_COUNTS).any(axis=1)
    if not (df[low] == LOW_COUNTS).all().all():
        raise ValueError("low_counts cells differ between runs")
    df = df[~low].apply(lump, min_frac=args.min_frac)

    order = label_order(df)
    title = (f"{args.title}\n{len(df):,} cells (score < {args.min_score} = unassigned); "
             f"{low.sum():,} low_counts cells not shown")
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    plot(df, order, label_colors(order), title, args.out)

    stem = os.path.splitext(args.out)[0]
    names = list(df.columns)
    for a, b in zip(names, names[1:]):
        pd.crosstab(df[a], df[b]).to_csv(f"{stem}_{slug(a)}_vs_{slug(b)}.csv")
    print("wrote", args.out)


if __name__ == "__main__":
    main()
