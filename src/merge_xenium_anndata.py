#!/usr/bin/env python3
"""
Merge per-sample Xenium AnnData objects (from xenium_to_anndata.py) into one.

Each input is rebuilt from its raw counts (layers["counts"]), so the
per-sample normalization is discarded. The merged object is normalized by
cell size with the same fixed constant as the per-sample objects:

    counts / cell_area * target_area (default 100 µm²), then log1p

Genes are outer-joined (a gene filtered out of one sample is 0 there), then
kept if detected in >= --min-cells merged cells. PCA/UMAP/Leiden are run on
the merged object without batch correction; obs["sample_id"] holds the batch.
Per-sample Leiden labels are kept as obs["leiden_sample"].

Usage:
  python src/merge_xenium_anndata.py --h5ads output/xenium/anndata/*.h5ad \
      --out output/xenium/anndata/merged/merged.h5ad
"""
import argparse
import os

import anndata as ad
import scanpy as sc

from xenium_to_anndata import cell_size_normalize, embed_and_cluster

OBS_KEEP = ["sample_id", "x_centroid", "y_centroid", "cell_area", "nucleus_area",
            "nucleus_count", "xenium_total_counts", "n_genes_by_counts",
            "total_counts", "leiden"]


def counts_only(path):
    """A sample reduced to raw counts, core cell metadata and spatial coords."""
    adata = ad.read_h5ad(path)
    obs = adata.obs[[c for c in OBS_KEEP if c in adata.obs]].copy()
    obs = obs.rename(columns={"leiden": "leiden_sample"})
    out = ad.AnnData(X=adata.layers["counts"].copy(), obs=obs,
                     var=adata.var[["gene_ids"]].copy())
    out.obsm["spatial"] = adata.obsm["spatial"].copy()
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ads", nargs="+", required=True)
    parser.add_argument("--out", required=True, help="output .h5ad")
    parser.add_argument("--target-area", type=float, default=100.0,
                        help="cell area (µm²) every cell is scaled to")
    parser.add_argument("--min-cells", type=int, default=5)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--resolution", type=float, default=0.3)
    args = parser.parse_args()

    samples = [counts_only(p) for p in sorted(args.h5ads)]
    adata = ad.concat(samples, join="outer", fill_value=0, merge="first")
    adata.obs["sample_id"] = adata.obs["sample_id"].astype("category")
    del samples

    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    adata.layers["counts"] = adata.X.copy()

    cell_size_normalize(adata, args.target_area)
    n_pcs = embed_and_cluster(adata, args.n_pcs, args.resolution)

    per_sample = adata.obs.groupby("sample_id", observed=True)["cell_area"].median()
    adata.uns["preprocessing"] = {
        "samples": list(adata.obs["sample_id"].cat.categories),
        "normalization": "counts / cell_area * target_area, log1p",
        "target_area": args.target_area,
        "sample_median_cell_area": {k: float(v) for k, v in per_sample.items()},
        "min_cells": args.min_cells,
        "n_pcs": n_pcs,
        "leiden_resolution": args.resolution,
        "batch_correction": "none",
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    adata.write_h5ad(args.out, compression="gzip")
    print(f"merged {adata.obs['sample_id'].nunique()} samples: {adata.n_obs} cells, "
          f"{adata.n_vars} genes, target_area={args.target_area:g} µm², "
          f"{adata.obs['leiden'].nunique()} clusters -> {args.out}")
    print(per_sample.round(1).to_string())


if __name__ == "__main__":
    main()
