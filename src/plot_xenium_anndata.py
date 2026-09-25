#!/usr/bin/env python3
"""
Standard plots for one processed Xenium AnnData (from xenium_to_anndata.py).

Writes into --outdir:

  qc_violin.png        genes, counts, cell and nucleus area per cell
  qc_spatial.png       counts, genes and cell area on the tissue
  pca_variance.png     variance explained per PC
  umap_qc.png          UMAP colored by counts, genes, cell area
  leiden.png           Leiden clusters on the UMAP and on the tissue
  spatial_leiden.png   Leiden clusters on the tissue, full size
  markers_dotplot.png  top marker genes per Leiden cluster (Wilcoxon)
  markers.csv          all marker statistics from rank_genes_groups

Usage:
  python src/plot_xenium_anndata.py --h5ad <sample>.h5ad --outdir <dir>
"""
import argparse
import os

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

QC_KEYS = ["n_genes_by_counts", "total_counts", "cell_area", "nucleus_area"]
QC_TITLES = {
    "n_genes_by_counts": "Genes per cell",
    "total_counts": "Transcripts per cell",
    "cell_area": "Cell area (µm²)",
    "nucleus_area": "Nucleus area (µm²)",
}
SEQ_CMAP = "viridis"


def spot_size(adata):
    """Point size that keeps dense and sparse sections readable."""
    return max(0.3, min(20.0, 120000 / adata.n_obs))


def style_spatial(ax):
    ax.set_aspect("equal")
    ax.invert_yaxis()  # Xenium y increases downward, as in the image
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")


def spatial(adata, color, ax, **kwargs):
    sc.pl.embedding(adata, basis="spatial", color=color, ax=ax, show=False,
                    size=spot_size(adata), frameon=True, **kwargs)
    style_spatial(ax)


def save(fig, outdir, name):
    fig.savefig(os.path.join(outdir, name), dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_qc(adata, sample, outdir):
    fig, axes = plt.subplots(1, len(QC_KEYS), figsize=(3.2 * len(QC_KEYS), 3.6))
    for ax, key in zip(axes, QC_KEYS):
        sc.pl.violin(adata, key, ax=ax, show=False, stripplot=False,
                     jitter=False, color="#8da0cb")
        ax.set_title(QC_TITLES[key])
        ax.set_ylabel("")
        ax.set_xticks([])
    fig.suptitle(f"{sample}: QC ({adata.n_obs:,} cells after filtering)")
    fig.tight_layout()
    save(fig, outdir, "qc_violin.png")

    keys = ["total_counts", "n_genes_by_counts", "cell_area"]
    fig, axes = plt.subplots(1, len(keys), figsize=(6 * len(keys), 5.5))
    for ax, key in zip(axes, keys):
        spatial(adata, key, ax, cmap=SEQ_CMAP, title=QC_TITLES[key])
    fig.suptitle(f"{sample}: QC on tissue")
    fig.tight_layout()
    save(fig, outdir, "qc_spatial.png")


def plot_embeddings(adata, sample, outdir):
    n_pcs = adata.obsm["X_pca"].shape[1]
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ratio = adata.uns["pca"]["variance_ratio"]
    ax.plot(range(1, n_pcs + 1), ratio, marker="o", markersize=4, linewidth=2,
            color="#3b6fb6")
    ax.set_xlabel("PC")
    ax.set_ylabel("Variance ratio")
    ax.set_title(f"{sample}: PCA variance explained")
    ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, outdir, "pca_variance.png")

    keys = ["total_counts", "n_genes_by_counts", "cell_area"]
    fig, axes = plt.subplots(1, len(keys), figsize=(5 * len(keys), 4.5))
    for ax, key in zip(axes, keys):
        sc.pl.umap(adata, color=key, ax=ax, show=False, cmap=SEQ_CMAP,
                   title=QC_TITLES[key])
    fig.suptitle(f"{sample}: UMAP")
    fig.tight_layout()
    save(fig, outdir, "umap_qc.png")

    n_clusters = adata.obs["leiden"].nunique()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sc.pl.umap(adata, color="leiden", ax=axes[0], show=False,
               legend_loc="on data", legend_fontsize=8, legend_fontoutline=2,
               title="Leiden on UMAP")
    spatial(adata, "leiden", axes[1], title="Leiden on tissue")
    fig.suptitle(f"{sample}: {n_clusters} Leiden clusters")
    fig.tight_layout()
    save(fig, outdir, "leiden.png")

    fig, ax = plt.subplots(figsize=(10, 10))
    spatial(adata, "leiden", ax, title=f"{sample}: Leiden on tissue")
    save(fig, outdir, "spatial_leiden.png")


def plot_markers(adata, sample, outdir, n_genes):
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=False,
                            pts=True)
    sc.get.rank_genes_groups_df(adata, group=None).to_csv(
        os.path.join(outdir, "markers.csv"), index=False)

    dp = sc.pl.rank_genes_groups_dotplot(adata, n_genes=n_genes, show=False,
                                         return_fig=True, standard_scale="var",
                                         title=f"{sample}: top {n_genes} markers per cluster")
    dp.savefig(os.path.join(outdir, "markers_dotplot.png"), dpi=150,
               bbox_inches="tight")
    plt.close("all")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--n-marker-genes", type=int, default=3)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    adata = ad.read_h5ad(args.h5ad)
    sample = adata.uns.get("sample_id", os.path.basename(args.h5ad).removesuffix(".h5ad"))

    plot_qc(adata, sample, args.outdir)
    plot_embeddings(adata, sample, args.outdir)
    plot_markers(adata, sample, args.outdir, args.n_marker_genes)
    print(f"{sample}: plots -> {args.outdir}")


if __name__ == "__main__":
    main()
