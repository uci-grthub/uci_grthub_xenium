#!/usr/bin/env python3
"""
Standard plots for one processed Xenium AnnData (from xenium_to_anndata.py).

Writes into --outdir:

  qc_violin.png        genes, counts, cell and nucleus area per cell
  qc_spatial.png       counts, genes and cell area on the tissue
  filtering_spatial.png  cells removed by QC filtering in red on the tissue,
                       next to the kept cells
  pca_variance.png     variance explained per PC
  umap_qc.png          UMAP colored by counts, genes, cell area
  leiden.png           Leiden clusters on the UMAP and on the tissue
  spatial_leiden.png   Leiden clusters on the tissue, full size
  markers_dotplot.png  top marker genes per Leiden cluster (Wilcoxon)
  markers.csv          all marker statistics from rank_genes_groups

With --merged (an object from merge_xenium_anndata.py; spatial coordinates
of different samples overlap, so no tissue plots):

  umap.png                 UMAP colored by sample and by Leiden cluster
  umap_by_sample.png       one panel per sample, its cells over all others in gray
  umap_qc.png              UMAP colored by counts, genes, cell area
  cluster_composition.png  fraction of each sample's cells per cluster
  markers_dotplot.png      top marker genes per Leiden cluster (Wilcoxon)
  markers.csv              all marker statistics from rank_genes_groups

Usage:
  python src/plot_xenium_anndata.py --h5ad <sample>.h5ad --outdir <dir> \
      --cells-csv <xenium_dir>/cells.csv.gz
  python src/plot_xenium_anndata.py --merged --h5ad <group>.h5ad --outdir <dir>
"""
import argparse
import os

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
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


def plot_filtering_spatial(adata, sample, cells_csv, outdir):
    """Cells removed by QC filtering highlighted on the tissue, next to the
    kept cells. Every detected cell comes from cells.csv.gz; kept cells are
    the ones that ended up in the h5ad."""
    cells = pd.read_csv(cells_csv, usecols=["cell_id", "x_centroid",
                                            "y_centroid"])
    cells["cell_id"] = cells["cell_id"].astype(str)
    kept_ids = set(adata.obs_names.str.removeprefix(f"{sample}_"))
    removed = ~cells["cell_id"].isin(kept_ids)
    n_all = len(cells)
    n_removed = int(removed.sum())
    size = max(0.3, min(20.0, 120000 / n_all))
    kept_cells = cells.loc[~removed]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))
    axes[0].scatter(kept_cells["x_centroid"], kept_cells["y_centroid"],
                    s=size, c="#bbbbbb", rasterized=True, linewidths=0)
    axes[0].scatter(cells.loc[removed, "x_centroid"],
                    cells.loc[removed, "y_centroid"],
                    s=size, c="#d62728", rasterized=True, linewidths=0)
    axes[0].set_title(f"all detected cells ({n_all:,})\nremoved by QC "
                      f"filtering in red ({n_removed:,}, "
                      f"{100 * n_removed / n_all:.1f}%)")
    axes[1].scatter(kept_cells["x_centroid"], kept_cells["y_centroid"],
                    s=size, c="#8da0cb", rasterized=True, linewidths=0)
    axes[1].set_title(f"after QC filtering ({len(kept_cells):,})")
    for ax in axes:
        style_spatial(ax)
    fig.suptitle(f"{sample}: QC filtering on tissue")
    fig.tight_layout()
    save(fig, outdir, "filtering_spatial.png")


def plot_umap_qc(adata, name, outdir):
    keys = ["total_counts", "n_genes_by_counts", "cell_area"]
    fig, axes = plt.subplots(1, len(keys), figsize=(5 * len(keys), 4.5))
    for ax, key in zip(axes, keys):
        sc.pl.umap(adata, color=key, ax=ax, show=False, cmap=SEQ_CMAP,
                   title=QC_TITLES[key])
    fig.suptitle(f"{name}: UMAP")
    fig.tight_layout()
    save(fig, outdir, "umap_qc.png")


def plot_merged(adata, name, outdir):
    samples = list(adata.obs["sample_id"].cat.categories)
    n_clusters = adata.obs["leiden"].nunique()

    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))
    sc.pl.umap(adata, color="sample_id", ax=axes[0], show=False,
               title="Sample")
    sc.pl.umap(adata, color="leiden", ax=axes[1], show=False,
               legend_loc="on data", legend_fontsize=8, legend_fontoutline=2,
               title=f"Leiden ({n_clusters} clusters)")
    correction = adata.uns.get("preprocessing", {}).get("batch_correction", "none")
    correction = "no batch correction" if correction == "none" else f"{correction} on sample_id"
    fig.suptitle(f"{name}: {len(samples)} samples, {adata.n_obs:,} cells, {correction}")
    fig.tight_layout()
    save(fig, outdir, "umap.png")

    ncols = min(4, len(samples))
    nrows = -(-len(samples) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows),
                             squeeze=False)
    umap = adata.obsm["X_umap"]
    colors = dict(zip(samples, adata.uns["sample_id_colors"]))
    for ax, sample in zip(axes.flat, samples):
        mine = (adata.obs["sample_id"] == sample).to_numpy()
        ax.scatter(umap[~mine, 0], umap[~mine, 1], s=0.5, c="#dddddd",
                   rasterized=True, linewidths=0)
        ax.scatter(umap[mine, 0], umap[mine, 1], s=1.5, c=colors[sample],
                   rasterized=True, linewidths=0)
        ax.set_title(f"{sample} ({mine.sum():,})", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in list(axes.flat)[len(samples):]:
        ax.axis("off")
    fig.suptitle(f"{name}: each sample over all cells")
    fig.tight_layout()
    save(fig, outdir, "umap_by_sample.png")

    plot_umap_qc(adata, name, outdir)

    frac = (adata.obs.groupby(["sample_id", "leiden"], observed=False).size()
            .unstack(fill_value=0))
    frac = frac.div(frac.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(9, 0.45 * len(samples) + 1.5))
    left = 0
    cluster_colors = adata.uns["leiden_colors"]
    for i, cluster in enumerate(frac.columns):
        ax.barh(frac.index, frac[cluster], left=left, color=cluster_colors[i],
                edgecolor="white", linewidth=0.5, label=cluster)
        left = left + frac[cluster]
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of sample's cells")
    ax.invert_yaxis()
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(title="Leiden", bbox_to_anchor=(1.01, 1), loc="upper left",
              ncol=2 if len(frac.columns) > 12 else 1, fontsize=8, frameon=False)
    ax.set_title(f"{name}: cluster composition per sample")
    save(fig, outdir, "cluster_composition.png")


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

    plot_umap_qc(adata, sample, outdir)

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
    parser.add_argument("--cells-csv",
                        help="raw cells.csv.gz (filtering_spatial.png; "
                             "required unless --merged)")
    parser.add_argument("--n-marker-genes", type=int, default=3)
    parser.add_argument("--merged", action="store_true",
                        help="input is a merged object; skip tissue plots")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    adata = ad.read_h5ad(args.h5ad)
    sample = adata.uns.get("sample_id", os.path.basename(args.h5ad).removesuffix(".h5ad"))

    if args.merged:
        plot_merged(adata, sample, args.outdir)
        plot_markers(adata, sample, args.outdir, args.n_marker_genes)
        print(f"{sample}: merged plots -> {args.outdir}")
        return

    if not args.cells_csv:
        parser.error("--cells-csv is required unless --merged")

    plot_qc(adata, sample, args.outdir)
    plot_filtering_spatial(adata, sample, args.cells_csv, args.outdir)
    plot_embeddings(adata, sample, args.outdir)
    plot_markers(adata, sample, args.outdir, args.n_marker_genes)
    print(f"{sample}: plots -> {args.outdir}")


if __name__ == "__main__":
    main()
