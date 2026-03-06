#!/usr/bin/env python3
# conda activate squidpy
"""
GO Biological Process over-representation analysis on pseudobulk DEGs.
Uses gseapy.enrichr (Enrichr API) with scipy.stats.fisher_exact as fallback.

Outputs (src/final_figures/output/figures/):
  pathway_dotplot_heatvsctrl.pdf   – per-cluster GO enrichment, heat vs ctrl
  pathway_dotplot_agedvsadult.pdf  – per-cluster GO enrichment, aged vs adult
  pathway_barplot_combined.pdf     – top GO terms, both comparisons overlaid

Saves intermediate results:
  src/final_figures/output/enrich_heatvsctrl.parquet
  src/final_figures/output/enrich_agedvsadult.parquet
"""

import os
import warnings
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

import gseapy as gp
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from plotnine import (
    ggplot, aes, geom_point, scale_size_continuous, scale_color_gradientn,
    scale_color_manual, scale_fill_manual, facet_wrap, labs,
    theme_bw, theme, element_text, element_blank, element_line,
    element_rect, coord_flip, position_dodge, geom_col, guides,
    guide_colorbar, guide_legend,
)

warnings.filterwarnings("ignore")

BASE    = "/dfs9/ucightf-lab/projects/ChatS/251006_0925Bio-06_ChatS_xenium"
DEG_PB  = os.path.join(BASE, "results/to_client/all/degs/pseudobulk")
DATADIR = os.path.join(BASE, "src/final_figures/output")
OUTDIR  = os.path.join(BASE, "src/final_figures/output/figures")
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

N_CLUSTERS = 17
P_GENE     = 0.05
FC_GENE    = 0.5
P_ENRICH   = 0.05
TOP_TERMS  = 5         # per cluster in faceted dot plot
TOP_GLOBAL = 25        # terms in combined bar plot
GENE_SET_DB = "GO_Biological_Process_2023"   # Enrichr library name


def shorten(s: str, n: int = 48) -> str:
    return s[:n - 1] + "…" if len(s) > n else s


# ── Load significant DEGs for one comparison ──────────────────────────────────
def load_sig(comparison: str) -> pd.DataFrame:
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


# ── Run Enrichr for one gene list ─────────────────────────────────────────────
def run_enrichr(genes: list[str], background: list[str] | None = None) -> pd.DataFrame | None:
    if len(genes) < 5:
        return None
    try:
        res = gp.enrichr(
            gene_list   = genes,
            gene_sets   = GENE_SET_DB,
            background  = background if background else "hsapiens_gene_ensembl",
            organism    = "mouse",
            outdir      = None,
            verbose     = False,
            cutoff      = 0.5,         # pre-filter; we reapply our threshold later
        )
        df = res.results.copy()
        df = df.rename(columns={
            "Term":            "term",
            "P-value":         "pvalue",
            "Adjusted P-value":"padj",
            "Genes":           "genes",
            "Overlap":         "overlap",
        })
        df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))
        df["gene_ratio"] = df["overlap"].apply(
            lambda x: int(x.split("/")[0]) / int(x.split("/")[1])
            if isinstance(x, str) and "/" in x else np.nan
        )
        return df[df["padj"] < P_ENRICH]
    except Exception as e:
        print(f"    Enrichr failed: {e}")
        return None


# ── Fallback: Fisher's exact ORA with local GO gene sets ──────────────────────
_CACHED_GENE_SETS: dict | None = None

def _load_go_gene_sets() -> dict:
    """Download GO BP mouse gene sets via gseapy if not cached."""
    global _CACHED_GENE_SETS
    if _CACHED_GENE_SETS is not None:
        return _CACHED_GENE_SETS
    cache_file = os.path.join(DATADIR, "_go_bp_mouse.json")
    if os.path.exists(cache_file):
        with open(cache_file) as fh:
            _CACHED_GENE_SETS = json.load(fh)
        return _CACHED_GENE_SETS
    try:
        print("    Downloading GO BP gene sets via gseapy…")
        sets = gp.get_library(name=GENE_SET_DB, organism="mouse")
        _CACHED_GENE_SETS = sets
        with open(cache_file, "w") as fh:
            json.dump(sets, fh)
        return sets
    except Exception as e:
        print(f"    Could not download gene sets: {e}")
        return {}


def run_fisher_ora(genes: list[str], background: list[str]) -> pd.DataFrame | None:
    gene_sets = _load_go_gene_sets()
    if not gene_sets:
        return None
    gene_set  = set(genes)
    bg_set    = set(background) if background else set()
    records   = []
    for term, members in gene_sets.items():
        term_set = set(members) & bg_set if bg_set else set(members)
        overlap  = gene_set & term_set
        if len(overlap) < 2:
            continue
        a = len(overlap)
        b = len(gene_set) - a
        c = len(term_set) - a
        d = len(bg_set) - a - b - c if bg_set else 10000 - a - b - c
        _, pval = fisher_exact([[a, b], [c, max(d, 0)]], alternative="greater")
        records.append({
            "term":       term,
            "pvalue":     pval,
            "genes":      "/".join(sorted(overlap)),
            "gene_ratio": a / len(gene_set),
            "n_overlap":  a,
        })
    if not records:
        return None
    df = pd.DataFrame(records)
    _, padj, _, _ = multipletests(df["pvalue"].values, method="fdr_bh")
    df["padj"]            = padj
    df["neg_log10_padj"]  = -np.log10(df["padj"].clip(lower=1e-300))
    return df[df["padj"] < P_ENRICH].sort_values("padj")


# ── Per-cluster enrichment ─────────────────────────────────────────────────────
def run_all_clusters(comparison: str) -> pd.DataFrame:
    sig_df   = load_sig(comparison)
    if sig_df.empty:
        return pd.DataFrame()

    panel_file = os.path.join(DATADIR, "panel_genes.csv")
    background = (pd.read_csv(panel_file)["gene"].tolist()
                  if os.path.exists(panel_file) else None)

    all_rows = []
    for cl in sorted(sig_df["cluster"].unique()):
        genes = (
            sig_df[(sig_df["cluster"] == cl)
                   & (sig_df["adj.P.Val"] < P_GENE)
                   & (sig_df["logFC"].abs() > FC_GENE)]["gene"]
            .dropna().unique().tolist()
        )
        if len(genes) < 5:
            continue

        print(f"    Cluster {cl}: {len(genes)} sig genes", flush=True)

        res = run_enrichr(genes, background)
        if res is None or res.empty:
            res = run_fisher_ora(genes, background or [])
        if res is None or res.empty:
            continue

        res["cluster"] = cl
        all_rows.append(res)

    return pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()


# ── Plotnine dot plot (per-cluster facets) ────────────────────────────────────
THEME_ENRICH = (
    theme_bw(base_size=8) +
    theme(
        panel_grid_minor   = element_blank(),
        panel_grid_major   = element_line(size=0.2, color="#EEEEEE"),
        strip_background   = element_rect(fill="#2C3E50", color="#2C3E50"),
        strip_text         = element_text(color="white", size=7),
        axis_text_y        = element_text(size=6.5),
        axis_text_x        = element_text(size=7),
        plot_title         = element_text(size=10, face="bold"),
    )
)


def dotplot_pdf(enrich_df: pd.DataFrame, title: str, outfile: str) -> None:
    if enrich_df.empty:
        print(f"  No enrichment data for: {title}"); return

    top_df = (
        enrich_df
        .sort_values("padj")
        .groupby("cluster", group_keys=False)
        .head(TOP_TERMS)
        .copy()
    )
    top_df["term"]    = top_df["term"].apply(lambda x: shorten(x))
    top_df["cluster"] = pd.Categorical(
        ["Cluster " + str(int(c)) for c in top_df["cluster"]],
        categories=["Cluster " + str(i)
                    for i in sorted(top_df["cluster"].astype(int).unique())],
        ordered=True,
    )

    n_clusters = top_df["cluster"].nunique()
    n_col      = min(5, n_clusters)
    n_row      = int(np.ceil(n_clusters / n_col))

    p = (
        ggplot(top_df, aes("gene_ratio", "term",
                           size="gene_ratio", color="neg_log10_padj"))
        + geom_point()
        + facet_wrap("cluster", ncol=n_col, scales="free_y")
        + scale_size_continuous(range=(1, 5), name="Gene ratio")
        + scale_color_gradientn(
            colors=["#FFF5F0", "#FC9272", "#CB181D"],
            name="-log₁₀(FDR)",
        )
        + labs(title=title, x="Gene ratio", y=None)
        + THEME_ENRICH
    )

    p.save(outfile,
           width=max(10, 4 * n_col),
           height=max(6, TOP_TERMS * n_row * 0.45 + 2),
           dpi=150, verbose=False)
    print(f"  Saved: {os.path.basename(outfile)}")


# ── Combined bubble plot (both comparisons) ───────────────────────────────────
def combined_bubble_pdf(heat_df: pd.DataFrame, age_df: pd.DataFrame,
                         outfile: str) -> None:
    def top_terms(df, label):
        if df.empty:
            return pd.DataFrame()
        return (
            df.groupby("term")["padj"].min().reset_index()
            .nsmallest(TOP_GLOBAL // 2, "padj")
            .assign(comparison=label)
            .merge(
                df.groupby("term").agg(
                    avg_nlp=("neg_log10_padj", "mean"),
                    n_clusters=("cluster", "nunique"),
                    gene_ratio=("gene_ratio", "mean"),
                ).reset_index(),
                on="term",
            )
        )

    combined = pd.concat([
        top_terms(heat_df, "Heat vs Ctrl"),
        top_terms(age_df,  "Aged vs Adult"),
    ], ignore_index=True)

    if combined.empty:
        print("  No data for combined bubble plot."); return

    # Order terms by total average significance
    term_order = (
        combined.groupby("term")["avg_nlp"].sum()
        .sort_values(ascending=True).index.tolist()
    )
    combined["term"] = pd.Categorical(
        combined["term"].apply(lambda x: shorten(x, 55)),
        categories=[shorten(t, 55) for t in term_order],
        ordered=True,
    )

    COMP_PAL = {"Heat vs Ctrl": "#FC8D62", "Aged vs Adult": "#8DA0CB"}

    p = (
        ggplot(combined, aes("avg_nlp", "term",
                             fill="comparison", size="n_clusters"))
        + geom_point(shape="o", alpha=0.85, color="#555555", stroke=0.3)
        + scale_fill_manual(values=COMP_PAL, name="Comparison")
        + scale_size_continuous(range=(2, 9), name="# Clusters")
        + labs(title="Top GO BP Terms – Main Comparisons",
               x="-log₁₀(FDR)", y=None)
        + theme_bw(base_size=9)
        + theme(
            panel_grid_minor  = element_blank(),
            panel_grid_major  = element_line(size=0.2, color="#EEEEEE"),
            axis_text_y       = element_text(size=7.5),
            plot_title        = element_text(size=11, face="bold"),
        )
    )

    n_terms = combined["term"].nunique()
    p.save(outfile, width=10, height=max(6, n_terms * 0.32 + 2),
           dpi=150, verbose=False)
    print(f"  Saved: {os.path.basename(outfile)}")


# ── Run ────────────────────────────────────────────────────────────────────────
for comparison, label in [
    ("comp1_heatvsctrl",  "Heat Stress vs Control"),
    ("comp1_agedvsadult", "Aged vs Adult"),
]:
    print(f"\n=== GO enrichment: {label} ===")
    parquet = os.path.join(DATADIR, f"enrich_{comparison.replace('comp1_','')}.parquet")
    enrich_df = run_all_clusters(comparison)

    if not enrich_df.empty:
        enrich_df.to_parquet(parquet, index=False)
        print(f"  Saved intermediate: {os.path.basename(parquet)}")
    else:
        print("  No enrichment results.")

    dotplot_pdf(
        enrich_df,
        f"GO BP Enrichment – {label}",
        os.path.join(OUTDIR, f"pathway_dotplot_{comparison.replace('comp1_','')}.pdf"),
    )

print("\n=== Combined comparison bubble plot ===")
heat_parquet = os.path.join(DATADIR, "enrich_heatvsctrl.parquet")
age_parquet  = os.path.join(DATADIR, "enrich_agedvsadult.parquet")
heat_df = pd.read_parquet(heat_parquet) if os.path.exists(heat_parquet) else pd.DataFrame()
age_df  = pd.read_parquet(age_parquet)  if os.path.exists(age_parquet)  else pd.DataFrame()
combined_bubble_pdf(
    heat_df, age_df,
    os.path.join(OUTDIR, "pathway_barplot_combined.pdf"),
)

print("\nDone.")
