#!/usr/bin/env python3
"""
Process one Xenium output directory into an AnnData object.

Standard scanpy Xenium pipeline:

  1. load cell_feature_matrix.h5 (gene expression features only; control
     probes and codewords are dropped) and cell metadata from cells.csv.gz
  2. QC metrics; keep cells with >= --min-genes genes and genes detected in
     >= --min-cells cells
  3. normalize by cell size: counts / cell_area * --target-area (a fixed
     constant, default 100 µm², like Seurat's scale.factor), then log1p
  4. PCA, neighbors, UMAP, Leiden

Raw counts are kept in layers["counts"]; cell centroids in obsm["spatial"].

cells.csv.gz is read rather than cells.parquet because several delivered
parquet files fail to read ("Unexpected end of stream").

Usage:
  python src/xenium_to_anndata.py --xenium-dir <output-XETG...> \
      --sample-id <id> --out <sample>.h5ad
"""
import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc


def load_xenium(xenium_dir, sample_id):
    adata = sc.read_10x_h5(os.path.join(xenium_dir, "cell_feature_matrix.h5"))
    adata.var_names_make_unique()

    cells = pd.read_csv(os.path.join(xenium_dir, "cells.csv.gz")).set_index("cell_id")
    cells.index = cells.index.astype(str)
    missing = adata.obs_names.difference(cells.index)
    if len(missing):
        raise ValueError(f"{len(missing)} matrix cells missing from cells.csv.gz")
    # scanpy QC recomputes total_counts from genes only; keep Xenium's own
    cells = cells.rename(columns={"total_counts": "xenium_total_counts"})
    adata.obs = cells.loc[adata.obs_names].copy()
    adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy()

    adata.obs["sample_id"] = sample_id
    adata.obs_names = sample_id + "_" + adata.obs_names
    adata.uns["sample_id"] = sample_id
    return adata


def cell_size_normalize(adata, target_area):
    """
    counts / cell_area * target_area, then log1p.

    target_area is a fixed constant (µm²) shared by every sample, so values
    are comparable across per-sample and merged objects.
    """
    area = adata.obs["cell_area"].to_numpy(dtype=float)
    scale = target_area / area
    adata.obs["size_factor"] = scale
    adata.X = adata.X.multiply(scale[:, None]).tocsr()
    sc.pp.log1p(adata)


def harmony_correct(adata, batch_key):
    """Harmony on X_pca over obs[batch_key] -> obsm["X_pca_harmony"].

    Calls harmonypy directly: scanpy's harmony_integrate transposes Z_corr,
    which harmonypy >= 2 already returns as cells x PCs.
    """
    import harmonypy

    z = np.asarray(harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs,
                                         [batch_key], verbose=False).Z_corr)
    if z.shape[0] != adata.n_obs:
        z = z.T
    adata.obsm["X_pca_harmony"] = z


def embed_and_cluster(adata, n_pcs, resolution, batch_key=None):
    """PCA, neighbors, UMAP, Leiden; returns the number of PCs used.

    With batch_key, neighbors (and so UMAP and Leiden) use the Harmony-corrected
    PCs; X_pca stays uncorrected.
    """
    n_pcs = min(n_pcs, adata.n_obs - 1, adata.n_vars - 1)
    sc.pp.pca(adata, n_comps=n_pcs)
    if batch_key:
        harmony_correct(adata, batch_key)
        sc.pp.neighbors(adata, use_rep="X_pca_harmony")
    else:
        sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution, key_added="leiden",
                 flavor="igraph", n_iterations=2, directed=False)
    return n_pcs


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--xenium-dir", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--out", required=True, help="output .h5ad")
    parser.add_argument("--min-genes", type=int, default=10)
    parser.add_argument("--min-cells", type=int, default=5)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--resolution", type=float, default=0.3)
    parser.add_argument("--target-area", type=float, default=100.0,
                        help="cell area (µm²) every cell is scaled to")
    args = parser.parse_args()

    adata = load_xenium(args.xenium_dir, args.sample_id)
    n_cells_raw = adata.n_obs

    sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    adata = adata[adata.obs["cell_area"] > 0].copy()

    adata.layers["counts"] = adata.X.copy()
    cell_size_normalize(adata, args.target_area)
    n_pcs = embed_and_cluster(adata, args.n_pcs, args.resolution)

    adata.uns["preprocessing"] = {
        "n_cells_raw": n_cells_raw,
        "min_genes": args.min_genes,
        "min_cells": args.min_cells,
        "normalization": "counts / cell_area * target_area, log1p",
        "target_area": args.target_area,
        "n_pcs": n_pcs,
        "leiden_resolution": args.resolution,
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    adata.write_h5ad(args.out, compression="gzip")
    print(f"{args.sample_id}: {n_cells_raw} -> {adata.n_obs} cells, "
          f"{adata.n_vars} genes, {adata.obs['leiden'].nunique()} clusters -> {args.out}")


if __name__ == "__main__":
    main()
