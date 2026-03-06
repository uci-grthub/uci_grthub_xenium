#!/usr/bin/env python3
# conda activate squidpy
"""
Network visualizations for biologically relevant gene modules.
Requires enrichment parquet files from p04_pathway_enrichment.py.

Uses networkx for graph construction; matplotlib for layout/rendering
(plotnine has no graph geometry).

Outputs (src/final_figures/output/figures/):
  network_cnet_heatvsctrl.pdf    – pathway–gene concept network, heat vs ctrl
  network_cnet_agedvsadult.pdf   – pathway–gene concept network, aged vs adult
  network_coexpr_heatstress.pdf  – gene co-expression network with Louvain
                                    community detection
"""

import os
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize, to_rgba
from matplotlib.cm import ScalarMappable
import networkx as nx
from scipy.stats import pearsonr
from sklearn.preprocessing import scale as sklearn_scale

warnings.filterwarnings("ignore")

BASE    = "/dfs9/ucightf-lab/projects/ChatS/251006_0925Bio-06_ChatS_xenium"
DEG_PB  = os.path.join(BASE, "results/to_client/all/degs/pseudobulk")
DATADIR = os.path.join(BASE, "src/final_figures/output")
OUTDIR  = os.path.join(BASE, "src/final_figures/output/figures")
os.makedirs(OUTDIR, exist_ok=True)

N_CLUSTERS  = 17
N_TERMS     = 12    # top GO terms in concept network
N_GENES_PER = 5     # max genes shown per term
CORR_CUTOFF = 0.45  # |r| for co-expression edges
MAX_GENES   = 60    # genes in co-expression network

_PAIRED = ["#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",
           "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#66C2A5","#FC8D62",
           "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494"]


def shorten(s: str, n: int = 42) -> str:
    return s[:n - 1] + "…" if len(s) > n else s


# ── Concept (pathway–gene bipartite) network ──────────────────────────────────
def build_cnet(parquet_file: str, n_terms: int = N_TERMS,
               n_genes_per: int = N_GENES_PER) -> nx.Graph | None:
    if not os.path.exists(parquet_file):
        print(f"  Parquet not found – run p04 first: {parquet_file}")
        return None

    df = pd.read_parquet(parquet_file)
    if df.empty:
        return None

    # Pick top terms by best (minimum) padj across clusters
    top_term_ids = (
        df.groupby("term")["padj"].min()
        .nsmallest(n_terms).index.tolist()
    )
    sub = (
        df[df["term"].isin(top_term_ids)]
        .sort_values("padj")
        .drop_duplicates("term")          # keep best-padj row per term
    )

    G = nx.Graph()

    for _, row in sub.iterrows():
        term_label = shorten(row["term"])
        G.add_node(term_label, node_type="pathway",
                   neg_log10_padj=row.get("neg_log10_padj", 1.0))

        genes_raw = row.get("genes", "")
        if not isinstance(genes_raw, str):
            continue
        genes = genes_raw.split("/")[:n_genes_per]
        for g in genes:
            g = g.strip()
            if not g:
                continue
            if g not in G:
                G.add_node(g, node_type="gene")
            G.add_edge(term_label, g)

    # Keep only nodes with edges
    isolated = [n for n, d in G.degree() if d == 0]
    G.remove_nodes_from(isolated)
    return G if G.number_of_nodes() > 0 else None


def draw_cnet_pdf(G: nx.Graph, title: str, outfile: str) -> None:
    if G is None:
        print(f"  No graph to draw for: {title}"); return

    pos = nx.spring_layout(G, k=2.5 / np.sqrt(G.number_of_nodes()),
                            seed=42, iterations=80)

    pathway_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "pathway"]
    gene_nodes    = [n for n, d in G.nodes(data=True) if d.get("node_type") == "gene"]

    # Significance colour for pathway nodes
    nlp_vals = np.array([G.nodes[n].get("neg_log10_padj", 1.0)
                         for n in pathway_nodes])
    nlp_norm = Normalize(vmin=0, vmax=max(nlp_vals.max(), 1.0))
    cmap_pw  = plt.get_cmap("YlOrRd")
    pw_colors = [cmap_pw(nlp_norm(v)) for v in nlp_vals]

    fig, ax = plt.subplots(figsize=(13, 10))
    ax.set_axis_off()

    # Edges
    nx.draw_networkx_edges(G, pos, ax=ax,
                           alpha=0.3, edge_color="#888888", width=0.8)

    # Gene nodes
    nx.draw_networkx_nodes(G, pos, nodelist=gene_nodes, ax=ax,
                           node_color="#3498DB", node_size=120, alpha=0.85)

    # Pathway nodes
    nx.draw_networkx_nodes(G, pos, nodelist=pathway_nodes, ax=ax,
                           node_color=pw_colors, node_size=520,
                           edgecolors="#333333", linewidths=0.8)

    # Labels – genes small, pathways slightly larger
    nx.draw_networkx_labels(G, pos, ax=ax,
                             labels={n: n for n in gene_nodes},
                             font_size=6, font_color="#1A1A1A")
    nx.draw_networkx_labels(G, pos, ax=ax,
                             labels={n: n for n in pathway_nodes},
                             font_size=6.5, font_weight="bold",
                             font_color="#7B241C",
                             bbox=dict(boxstyle="round,pad=0.2",
                                       fc="white", ec="none", alpha=0.7))

    # Colorbar for pathway significance
    sm = ScalarMappable(cmap=cmap_pw, norm=nlp_norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.35, pad=0.01,
                        orientation="vertical")
    cbar.set_label("-log₁₀(FDR)", fontsize=9)

    # Legend
    handles = [
        mpatches.Patch(color="#3498DB", label="Gene"),
        mpatches.Patch(color=cmap_pw(0.85), label="GO term"),
    ]
    ax.legend(handles=handles, loc="lower left", fontsize=9,
              frameon=True, framealpha=0.9)

    ax.set_title(title, fontsize=13, fontweight="bold", pad=10)
    fig.tight_layout()
    fig.savefig(outfile, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved: {os.path.basename(outfile)}")


# ── Gene co-expression network ────────────────────────────────────────────────
def load_logfc_profiles(comparison: str) -> pd.DataFrame:
    """Build gene × cluster logFC matrix from pseudobulk DEG files."""
    rows = {}
    for cl in range(N_CLUSTERS):
        f = os.path.join(DEG_PB, f"DEG_pb_{comparison}_cluster{cl}.csv")
        if not os.path.exists(f):
            continue
        df = pd.read_csv(f, index_col=0)
        rows[cl] = df["logFC"]
    if not rows:
        return pd.DataFrame()
    mat = pd.DataFrame(rows).T    # clusters × genes
    mat.columns.name = "gene"
    return mat.fillna(0)


def draw_coexpr_pdf(comparison: str, outfile: str) -> None:
    print(f"  Building co-expression network ({comparison})…")
    mat = load_logfc_profiles(comparison)
    if mat.empty:
        print("  No data."); return

    # Keep most variable genes
    gene_var = mat.abs().mean(axis=0)
    top_genes = gene_var.nlargest(MAX_GENES).index.tolist()
    mat = mat[top_genes]

    # Pairwise Pearson correlation
    corr = mat.corr(method="pearson")

    # Build graph from edges above threshold
    G = nx.Graph()
    G.add_nodes_from(top_genes)
    for i, g1 in enumerate(top_genes):
        for j, g2 in enumerate(top_genes):
            if j <= i:
                continue
            r = corr.loc[g1, g2]
            if abs(r) >= CORR_CUTOFF:
                G.add_edge(g1, g2, weight=r)

    # Remove isolated nodes
    isolated = [n for n, d in G.degree() if d == 0]
    G.remove_nodes_from(isolated)

    if G.number_of_nodes() < 3:
        print(f"  Too few nodes after filtering (r ≥ {CORR_CUTOFF}). "
              "Lowering threshold to 0.3…")
        CORR_CUTOFF2 = 0.3
        for i, g1 in enumerate(top_genes):
            for j, g2 in enumerate(top_genes):
                if j <= i:
                    continue
                r = corr.loc[g1, g2]
                if abs(r) >= CORR_CUTOFF2 and not G.has_edge(g1, g2):
                    G.add_edge(g1, g2, weight=r)
        isolated = [n for n, d in G.degree() if d == 0]
        G.remove_nodes_from(isolated)

    if G.number_of_nodes() < 3:
        print("  Graph too sparse to visualise."); return

    # Louvain community detection
    try:
        communities = nx.community.louvain_communities(G, seed=42)
    except AttributeError:
        # older networkx
        communities = list(nx.community.greedy_modularity_communities(G))

    node_community = {}
    for idx, comm in enumerate(communities):
        for node in comm:
            node_community[node] = idx

    n_comms  = len(communities)
    comm_pal = {i: _PAIRED[i % len(_PAIRED)] for i in range(n_comms)}
    node_colors = [comm_pal[node_community[n]] for n in G.nodes()]

    # Node size ∝ degree (hub genes larger)
    degrees    = dict(G.degree())
    node_sizes = [60 + 25 * degrees[n] for n in G.nodes()]

    # Mean |logFC| across clusters as node attribute
    mean_fc = mat.abs().mean(axis=0)

    # Edge colors: positive=red, negative=blue
    edge_colors = ["#C0392B" if G[u][v]["weight"] > 0 else "#2980B9"
                   for u, v in G.edges()]
    edge_widths = [0.4 + abs(G[u][v]["weight"]) * 1.5
                   for u, v in G.edges()]

    pos = nx.spring_layout(G, weight="weight",
                            k=2.8 / np.sqrt(G.number_of_nodes()),
                            seed=7, iterations=120)

    fig, ax = plt.subplots(figsize=(13, 10))
    ax.set_axis_off()

    nx.draw_networkx_edges(G, pos, ax=ax,
                            edge_color=edge_colors,
                            width=edge_widths, alpha=0.55)
    nx.draw_networkx_nodes(G, pos, ax=ax,
                            node_color=node_colors,
                            node_size=node_sizes, alpha=0.88,
                            edgecolors="#333333", linewidths=0.5)
    nx.draw_networkx_labels(G, pos, ax=ax,
                             font_size=6.5, font_color="#1A1A1A")

    # Community legend (top communities only)
    top_comm = sorted(communities, key=len, reverse=True)[:6]
    legend_handles = []
    for idx, comm in enumerate(top_comm):
        hub = max(comm, key=lambda n: degrees.get(n, 0))
        legend_handles.append(
            mpatches.Patch(color=comm_pal[node_community[hub]],
                           label=f"Module {idx+1} (hub: {hub})")
        )
    ax.legend(handles=legend_handles, loc="lower left", fontsize=8,
              frameon=True, framealpha=0.9, title="Louvain modules",
              title_fontsize=9)

    # Edge type legend
    edge_handles = [
        mpatches.Patch(color="#C0392B", label="Positive corr."),
        mpatches.Patch(color="#2980B9", label="Negative corr."),
    ]
    ax.legend(handles=edge_handles, loc="lower right", fontsize=8,
              frameon=True, framealpha=0.9)

    label = comparison.replace("comp1_", "").replace("_", " ").title()
    ax.set_title(f"Gene Co-expression Network – {label}\n"
                 f"(|r| ≥ {CORR_CUTOFF}, n={G.number_of_nodes()} genes, "
                 f"Louvain modules)",
                 fontsize=12, fontweight="bold", pad=10)

    fig.tight_layout()
    fig.savefig(outfile, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved: {os.path.basename(outfile)}")


# ── Run ────────────────────────────────────────────────────────────────────────
print("=== Concept networks ===")
for comp, label in [
    ("heatvsctrl",  "Heat Stress vs Control"),
    ("agedvsadult", "Aged vs Adult"),
]:
    G = build_cnet(os.path.join(DATADIR, f"enrich_{comp}.parquet"))
    draw_cnet_pdf(
        G, f"Pathway–Gene Concept Network – {label}",
        os.path.join(OUTDIR, f"network_cnet_{comp}.pdf"),
    )

print("\n=== Co-expression network ===")
draw_coexpr_pdf(
    "comp1_heatvsctrl",
    os.path.join(OUTDIR, "network_coexpr_heatstress.pdf"),
)

print("\nDone.")
