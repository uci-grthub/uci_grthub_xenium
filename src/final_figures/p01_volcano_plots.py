#!/usr/bin/env python3
# conda activate squidpy
"""
Volcano plots for pseudobulk and cluster-wise (Wilcoxon) DEGs using plotnine.

Outputs (src/final_figures/output/figures/):
  volcano_pb_<comparison>.pdf   – pseudobulk, faceted by Leiden cluster
  volcano_std_<file>.pdf        – cluster-wise Wilcoxon, faceted by cluster
"""

import os
import re
import warnings

import numpy as np
import pandas as pd
from plotnine import (
    ggplot, aes, geom_point, geom_hline, geom_vline, geom_text,
    scale_color_manual, scale_size_manual, facet_wrap, labs,
    theme_bw, theme, element_text, element_blank, element_rect,
    element_line, guides, guide_legend, coord_cartesian, save_as_pdf_pages,
)
from plotnine.geoms import geom_label
import matplotlib
matplotlib.use("Agg")

warnings.filterwarnings("ignore")

BASE   = "/dfs9/ucightf-lab/projects/ChatS/251006_0925Bio-06_ChatS_xenium"
DEG_PB = os.path.join(BASE, "results/to_client/all/degs/pseudobulk")
DEG_ST = os.path.join(BASE, "results/to_client/all/degs/standard")
OUTDIR = os.path.join(BASE, "src/final_figures/output/figures")
os.makedirs(OUTDIR, exist_ok=True)

N_CLUSTERS = 17
FC_THRESH  = 0.5
P_THRESH   = 0.05
TOP_N      = 6       # genes labeled per panel

COL_UP = "#C0392B"
COL_DN = "#2980B9"
COL_NS = "#BDBDBD"

THEME = (
    theme_bw(base_size=7) +
    theme(
        panel_grid_minor   = element_blank(),
        panel_grid_major   = element_line(size=0.25, color="#EEEEEE"),
        strip_background   = element_rect(fill="#2C3E50", color="#2C3E50"),
        strip_text         = element_text(color="white", size=6.5),
        legend_position    = "none",
        plot_title         = element_text(size=9, face="bold"),
        axis_title         = element_text(size=7),
        axis_text          = element_text(size=6),
    )
)


def _sig_col(fc: pd.Series, pval: pd.Series) -> pd.Series:
    return pd.Series(
        np.where(
            (pval < P_THRESH) & (fc > FC_THRESH), "Up",
            np.where(
                (pval < P_THRESH) & (fc < -FC_THRESH), "Down", "NS"
            ),
        ),
        index=fc.index,
    )


def _top_labels(df: pd.DataFrame, n: int = TOP_N) -> pd.DataFrame:
    return (
        df[df["sig"] != "NS"]
        .nlargest(n, "nlp")
        .copy()
    )


def volcano_df(df: pd.DataFrame, fc_col="logFC", pval_col="adj.P.Val",
               gene_col="gene") -> pd.DataFrame:
    df = df.copy().rename(columns={fc_col: "logFC", pval_col: "adj.P.Val",
                                    gene_col: "gene"})
    df["logFC"]     = pd.to_numeric(df["logFC"],     errors="coerce")
    df["adj.P.Val"] = pd.to_numeric(df["adj.P.Val"], errors="coerce")
    df = df.dropna(subset=["logFC", "adj.P.Val"])
    df["nlp"] = -np.log10(df["adj.P.Val"].clip(lower=1e-300))
    df["sig"] = _sig_col(df["logFC"], df["adj.P.Val"])
    return df


# ── Multi-panel volcano (faceted by cluster) ──────────────────────────────────
def faceted_volcano(combined: pd.DataFrame, title: str, outfile: str,
                    cluster_col: str = "cluster",
                    n_col: int = 6) -> None:
    """combined must have cols: gene, logFC, adj.P.Val, nlp, sig, <cluster_col>"""
    if combined.empty:
        print(f"  [skip] no data for {title}")
        return

    # Top labels per cluster
    label_df = (
        combined[combined["sig"] != "NS"]
        .groupby(cluster_col, group_keys=False)
        .apply(lambda g: g.nlargest(TOP_N, "nlp"))
        .reset_index(drop=True)
    )

    p = (
        ggplot(combined, aes("logFC", "nlp", color="sig"))
        + geom_point(size=0.4, alpha=0.5)
        + geom_hline(yintercept=-np.log10(P_THRESH),
                     linetype="dashed", color="#7F8C8D", size=0.35)
        + geom_vline(xintercept=[-FC_THRESH, FC_THRESH],
                     linetype="dashed", color="#7F8C8D", size=0.35)
        + geom_text(
            data=label_df,
            mapping=aes(label="gene"),
            size=4, color="black", ha="left", adjust_text={"avoid_points": True,
                                                             "expand_points": (1.2, 1.2)}
        )
        + scale_color_manual(values={"Up": COL_UP, "Down": COL_DN, "NS": COL_NS})
        + facet_wrap(cluster_col, ncol=n_col)
        + labs(title=title, x="log2 FC", y="-log₁₀(adj. p-value)")
        + THEME
    )

    n_clusters = combined[cluster_col].nunique()
    nrows = int(np.ceil(n_clusters / n_col))
    try:
        p.save(outfile, width=3.2 * n_col, height=3.2 * nrows, dpi=150,
               verbose=False)
        print(f"  Saved: {os.path.basename(outfile)}")
    except Exception as e:
        print(f"  Error saving {outfile}: {e}")


# ── Pseudobulk volcanos ────────────────────────────────────────────────────────
def pb_volcano_pdf(comparison: str) -> None:
    pattern = re.compile(rf"DEG_pb_{re.escape(comparison)}_cluster(\d+)\.csv$")
    files = sorted(
        [(int(m.group(1)), f)
         for f in os.listdir(DEG_PB)
         if (m := pattern.match(f))],
        key=lambda x: x[0],
    )
    if not files:
        print(f"  No files for: {comparison}")
        return

    frames = []
    for cl_n, fname in files:
        df = pd.read_csv(os.path.join(DEG_PB, fname), index_col=0)
        df["gene"]    = df.index
        df["cluster"] = f"Cluster {cl_n:02d}"
        frames.append(df)

    combined = volcano_df(pd.concat(frames, ignore_index=True))
    combined["cluster"] = pd.Categorical(
        combined["cluster"],
        categories=sorted(combined["cluster"].unique()),
        ordered=True,
    )

    label = (comparison
             .replace("comp1_", "").replace("comp2_", "")
             .replace("comp3_", "").replace("_", " ").title())
    outfile = os.path.join(OUTDIR, f"volcano_pb_{comparison}.pdf")
    faceted_volcano(combined, f"Pseudobulk DEGs – {label}", outfile)


PB_COMPS = [
    "comp1_heatvsctrl",
    "comp1_agedvsadult",
    "comp2_adultonly_heatvsctrl",
    "comp2_agedonly_heatvsctrl",
]

print("=== Pseudobulk volcano plots ===")
for comp in PB_COMPS:
    print(f"  {comp}")
    pb_volcano_pdf(comp)


# ── Standard (Wilcoxon) cluster-wise volcanos ──────────────────────────────────
def std_volcano_pdf(xlsx_path: str) -> None:
    if not os.path.exists(xlsx_path):
        print(f"  Not found: {xlsx_path}")
        return

    xl = pd.ExcelFile(xlsx_path)
    frames = []
    for sh in xl.sheet_names:
        df = xl.parse(sh)
        if "names" not in df.columns:
            continue
        df = df.rename(columns={
            "names":          "gene",
            "logfoldchanges": "logFC",
            "pvals_adj":      "adj.P.Val",
        })
        df["cluster"] = f"Cluster {int(sh):02d}"
        frames.append(df[["gene", "logFC", "adj.P.Val", "cluster"]])

    if not frames:
        return

    combined = volcano_df(pd.concat(frames, ignore_index=True),
                          fc_col="logFC", pval_col="adj.P.Val")
    combined["cluster"] = pd.Categorical(
        combined["cluster"],
        categories=sorted(combined["cluster"].unique()),
        ordered=True,
    )

    stem = os.path.splitext(os.path.basename(xlsx_path))[0]
    label = re.sub(r"_ranked.*", "", stem).replace("_", " ").title()
    outfile = os.path.join(
        OUTDIR, f"volcano_std_{stem}.pdf"
    )
    faceted_volcano(combined, f"Wilcoxon DEGs – {label}", outfile,
                    n_col=6)


STD_FILES = [
    "adult_heatstress_v_control_ranked_genes_per_cluster.xlsx",
    "aged_heatstress_v_control_ranked_genes_per_cluster.xlsx",
    "all_adult_v_aged_ranked_genes_per_cluster.xlsx",
    "all_control_v_heatstress_ranked_genes_per_cluster.xlsx",
]

print("\n=== Standard (Wilcoxon) volcano plots ===")
for f in STD_FILES:
    print(f"  {f}")
    std_volcano_pdf(os.path.join(DEG_ST, f))

print("\nDone.")
