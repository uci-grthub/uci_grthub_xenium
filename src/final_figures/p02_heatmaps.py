#!/usr/bin/env python3
# conda activate squidpy
"""
Heatmaps of top DEGs and curated pathway gene sets.

Uses seaborn clustermap for annotated heatmaps (plotnine's geom_tile lacks
hierarchical clustering + annotation sidebars).

Outputs (src/final_figures/output/figures/):
  heatmap_logFC_heatvsctrl.pdf      – logFC per gene × cluster, heat vs ctrl
  heatmap_logFC_agedvsadult.pdf     – logFC per gene × cluster, aged vs adult
  heatmap_expression_curated.pdf    – mean log-norm expression of curated
                                       pathway genes (requires 00_extract first)
"""

import os
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram

warnings.filterwarnings("ignore")

BASE    = "/dfs9/ucightf-lab/projects/ChatS/251006_0925Bio-06_ChatS_xenium"
DEG_PB  = os.path.join(BASE, "results/to_client/all/degs/pseudobulk")
DATADIR = os.path.join(BASE, "src/final_figures/output")
OUTDIR  = os.path.join(BASE, "src/final_figures/output/figures")
os.makedirs(OUTDIR, exist_ok=True)

N_CLUSTERS = 17
P_THRESH   = 0.05
FC_THRESH  = 0.25
TOP_N      = 6       # top genes per cluster to keep

# ── Color maps ────────────────────────────────────────────────────────────────
FC_CMAP    = "RdBu_r"
EXPR_CMAP  = "RdYlBu_r"
AGE_PAL    = {"adult": "#2166AC", "aged": "#D73027"}
COND_PAL   = {"control": "#66C2A5", "heatstress": "#FC8D62"}
AREA_PAL   = {"area_1": "#7B2D8B", "area_2": "#E6932A",
              "area_3": "#2E8B57", "area_4": "#C0392B"}
PATHWAY_PAL = sns.color_palette("Set2", 10)


# ── Helper: load all cluster files for one pseudobulk comparison ──────────────
def load_pb_comp(comparison: str) -> pd.DataFrame:
    frames = []
    for cl in range(N_CLUSTERS):
        f = os.path.join(DEG_PB, f"DEG_pb_{comparison}_cluster{cl}.csv")
        if not os.path.exists(f):
            continue
        df = pd.read_csv(f, index_col=0)
        df["gene"]    = df.index
        df["cluster"] = cl
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


# ── Helper: build gene × cluster logFC matrix ─────────────────────────────────
def build_logfc_matrix(df: pd.DataFrame) -> pd.DataFrame | None:
    sig = df[(df["adj.P.Val"] < P_THRESH) & (df["logFC"].abs() > FC_THRESH)]
    if sig.empty:
        return None

    # select TOP_N genes per cluster by absolute logFC (avoid pandas
    # DataFrame.nlargest(..., key=...) which may not be available on older
    # pandas versions). Use sort_values on the absolute values instead.
    top_genes = (
        sig.groupby("cluster", group_keys=False)
        .apply(lambda g: g.loc[
            g["logFC"].abs().sort_values(ascending=False).head(TOP_N).index
        ])
        .reset_index(drop=True)
        .drop_duplicates("gene")["gene"]
        .tolist()
    )

    mat = (
        df[df["gene"].isin(top_genes)]
        .pivot_table(index="gene", columns="cluster",
                     values="logFC", aggfunc="mean")
        .fillna(0)
    )
    mat.columns = [f"C{c}" for c in sorted(mat.columns)]
    return mat


# ── logFC heatmap with seaborn clustermap ─────────────────────────────────────
def logfc_heatmap_pdf(comparison: str, outfile: str, title: str) -> None:
    print(f"  {comparison}")
    df  = load_pb_comp(comparison)
    if df.empty:
        print("    No data found."); return

    mat = build_logfc_matrix(df)
    if mat is None or mat.empty:
        print("    No significant genes."); return

    vmax = min(np.percentile(mat.abs().values, 95), 3.0)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    nrow, ncol = mat.shape
    fig_h = max(6, nrow * 0.22 + 2)
    fig_w = max(6, ncol * 0.5  + 2.5)

    g = sns.clustermap(
        mat,
        cmap        = FC_CMAP,
        norm        = norm,
        figsize     = (fig_w, fig_h),
        row_cluster = nrow > 1,
        col_cluster = False,
        linewidths  = 0,
        yticklabels = True,
        xticklabels = True,
        cbar_kws    = {"label": "log2 FC", "shrink": 0.5,
                       "orientation": "horizontal"},
    )
    g.ax_heatmap.set_xlabel("Leiden cluster", fontsize=9)
    g.ax_heatmap.set_ylabel("Gene", fontsize=9)
    g.ax_heatmap.tick_params(axis="y", labelsize=6)
    g.ax_heatmap.tick_params(axis="x", labelsize=8)
    g.fig.suptitle(title, y=1.01, fontsize=11, fontweight="bold")

    g.fig.savefig(outfile, bbox_inches="tight", dpi=150)
    plt.close("all")
    print(f"    Saved: {os.path.basename(outfile)}")


# ── Curated-pathway expression heatmap ────────────────────────────────────────
def curated_expression_heatmap_pdf(outfile: str) -> None:
    expr_file    = os.path.join(DATADIR, "cluster_mean_expr.csv")
    curated_file = os.path.join(DATADIR, "curated_gene_sets.csv")

    if not os.path.exists(expr_file):
        print("  Skipping curated heatmap – run 00_extract_h5ad_data.py first.")
        return

    print("  Building curated pathway expression heatmap...")
    expr    = pd.read_csv(expr_file)
    curated = pd.read_csv(curated_file)

    gene_cols = [c for c in curated["gene"].unique() if c in expr.columns]
    if not gene_cols:
        print("  No curated genes found in expression matrix."); return

    # Column label: C<cluster>_<age2>_<cond4>
    expr = expr.copy()
    expr["col_label"] = (
        "C" + expr["leiden_harmony"].astype(str)
        + "_" + expr["age"].str[:2]
        + "_" + expr["condition"].str[:4]
    )

    mat = (
        expr.set_index("col_label")[gene_cols]
        .T   # genes × groups
        .astype(float)
    )

    # Z-score per gene across groups
    mat_z = mat.subtract(mat.mean(axis=1), axis=0).divide(
        mat.std(axis=1).replace(0, 1), axis=0
    )

    # ── Column annotations ─────────────────────────────────────────────────────
    col_meta = (
        expr.set_index("col_label")[["age", "condition"]]
        .loc[mat_z.columns]
    )
    age_colors  = col_meta["age"].map(AGE_PAL)
    cond_colors = col_meta["condition"].map(COND_PAL)
    col_colors  = pd.DataFrame(
        {"Age": age_colors, "Condition": cond_colors}
    )

    # ── Row annotations (pathway membership) ──────────────────────────────────
    pathway_list = sorted(curated["pathway"].unique())
    pw_palette   = {p: PATHWAY_PAL[i % len(PATHWAY_PAL)]
                    for i, p in enumerate(pathway_list)}
    gene_pathway = (
        curated[curated["gene"].isin(gene_cols)]
        .drop_duplicates("gene")
        .set_index("gene")["pathway"]
    )
    # reindex to match mat_z row order
    gene_pathway = gene_pathway.reindex(mat_z.index).fillna("Other")
    row_colors   = gene_pathway.map(pw_palette)

    vmax  = min(mat_z.abs().quantile(0.95).max(), 2.5)
    norm  = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    nrow, ncol = mat_z.shape
    fig_h = max(8,  nrow * 0.25 + 3)
    fig_w = max(12, ncol * 0.18 + 4)

    g = sns.clustermap(
        mat_z,
        cmap         = EXPR_CMAP,
        norm         = norm,
        figsize      = (fig_w, fig_h),
        row_cluster  = False,    # keep pathway grouping order
        col_cluster  = False,
        col_colors   = col_colors,
        row_colors   = row_colors,
        linewidths   = 0,
        yticklabels  = True,
        xticklabels  = True,
        cbar_kws     = {"label": "z-score", "shrink": 0.4,
                        "orientation": "horizontal"},
    )
    g.ax_heatmap.tick_params(axis="y", labelsize=6)
    g.ax_heatmap.tick_params(axis="x", labelsize=5.5, rotation=90)
    g.ax_heatmap.set_xlabel("Cluster × Age × Condition", fontsize=9)
    g.ax_heatmap.set_ylabel("")

    # ── Pathway legend ─────────────────────────────────────────────────────────
    handles = [mpatches.Patch(color=col, label=pw)
               for pw, col in pw_palette.items() if pw in gene_pathway.values]
    g.fig.legend(handles=handles, title="Pathway", loc="lower left",
                 fontsize=7, title_fontsize=8,
                 bbox_to_anchor=(0.0, 0.0), frameon=True)

    # Age / condition legend
    age_handles  = [mpatches.Patch(color=c, label=k) for k, c in AGE_PAL.items()]
    cond_handles = [mpatches.Patch(color=c, label=k) for k, c in COND_PAL.items()]
    leg2 = g.fig.legend(handles=age_handles + cond_handles,
                         title="Age / Condition", loc="lower right",
                         fontsize=7, title_fontsize=8,
                         bbox_to_anchor=(1.0, 0.0), frameon=True)

    g.fig.suptitle("Curated Pathway Gene Sets – Mean Expression (z-scored)",
                   y=1.01, fontsize=11, fontweight="bold")
    g.fig.savefig(outfile, bbox_inches="tight", dpi=150)
    plt.close("all")
    print(f"    Saved: {os.path.basename(outfile)}")


# ── Run ────────────────────────────────────────────────────────────────────────
print("=== LogFC heatmaps ===")
logfc_heatmap_pdf(
    "comp1_heatvsctrl",
    os.path.join(OUTDIR, "heatmap_logFC_heatvsctrl.pdf"),
    "Top DEGs – Heat Stress vs Control (pseudobulk logFC)",
)
logfc_heatmap_pdf(
    "comp1_agedvsadult",
    os.path.join(OUTDIR, "heatmap_logFC_agedvsadult.pdf"),
    "Top DEGs – Aged vs Adult (pseudobulk logFC)",
)

print("\n=== Curated expression heatmap ===")
curated_expression_heatmap_pdf(
    os.path.join(OUTDIR, "heatmap_expression_curated.pdf")
)

print("\nDone.")
