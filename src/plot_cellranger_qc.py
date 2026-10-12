#!/usr/bin/env python3
"""
Standard per-sample QC violin plots for one Cell Ranger sample.

Reads the sample's filtered_feature_bc_matrix.h5 and writes a single
figure with one violin panel per QC metric:

  genes per cell, transcripts per cell, % mitochondrial counts

Usage:
  python src/plot_cellranger_qc.py --matrix <sample>/outs/filtered_feature_bc_matrix.h5 \
      --sample-id <sample> --out output/qc/violin/<sample>.png
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

QC_KEYS = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
QC_TITLES = {
    "n_genes_by_counts": "Genes per cell",
    "total_counts": "Transcripts per cell",
    "pct_counts_mt": "% mitochondrial counts",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--matrix", required=True,
                        help="Cell Ranger filtered_feature_bc_matrix.h5")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--out", required=True, help="output PNG path")
    parser.add_argument("--mt-prefix", default="MT-",
                        help="mitochondrial gene name prefix (default: MT-)")
    args = parser.parse_args()

    adata = sc.read_10x_h5(args.matrix)
    adata.var_names_make_unique()

    adata.var["mt"] = adata.var_names.str.startswith(args.mt_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                               log1p=False, inplace=True)

    fig, axes = plt.subplots(1, len(QC_KEYS),
                             figsize=(3.2 * len(QC_KEYS), 3.6))
    for ax, key in zip(axes, QC_KEYS):
        sc.pl.violin(adata, key, ax=ax, show=False, stripplot=False,
                     jitter=False, color="#8da0cb")
        ax.set_title(QC_TITLES[key])
        ax.set_ylabel("")
        ax.set_xticks([])
    fig.suptitle(f"{args.sample_id}: QC "
                 f"({adata.n_obs:,} cells, unfiltered matrix)")
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"{args.sample_id}: QC violin plot -> {args.out}")


if __name__ == "__main__":
    main()
