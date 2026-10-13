#!/usr/bin/env python3
"""
Export raw counts from an .h5ad as a 10x-style Matrix Market directory for R
(Seurat::ReadMtx), the hand-off used by rules/azimuth.smk.

Writes into --outdir:

  matrix.mtx.gz    genes x cells, integer counts
  features.tsv.gz  gene names, one per line
  barcodes.tsv.gz  cell ids, one per line
  metadata.csv     --obs-columns per cell (index = cell id)

Two uses:

  query (a per-sample Xenium AnnData from xenium_to_anndata.py):
    python src/h5ad_to_mtx.py --h5ad <sample>.h5ad --layer counts --outdir <dir>

  reference (a CELLxGENE .h5ad: raw counts in raw/X, Ensembl ids as var index):
    python src/h5ad_to_mtx.py --h5ad reference.h5ad --layer raw \
        --panel <xenium_dir>/gene_panel.json --label-column cell_type \
        --keep-labels '["fibroblast", ...]' --max-cells-per-label 2000 \
        --obs-columns cell_type --outdir <dir>

With --panel, genes are matched to the Xenium panel by Ensembl id and renamed
to the panel's gene symbols, so reference and query share feature names. Cells
are selected (--keep-labels, --max-cells-per-label) before the count matrix is
read, so only those rows of a large atlas are loaded.
"""
import argparse
import gzip
import json
import os

import h5py
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from anndata.io import read_elem, sparse_dataset


def panel_genes(path):
    """Ensembl id -> gene symbol for the gene targets of a Xenium gene_panel.json."""
    with open(path) as handle:
        targets = json.load(handle)["payload"]["targets"]
    return {t["type"]["data"]["id"]: t["type"]["data"]["name"]
            for t in targets if t["type"]["descriptor"] == "gene"}


def select_cells(obs, label_column, keep_labels, max_per_label, seed):
    """Row positions to export: cells with a kept label, at most max_per_label each."""
    rows = np.arange(len(obs))
    if keep_labels:
        labels = obs[label_column].astype(str)
        missing = sorted(set(keep_labels) - set(labels))
        if missing:
            print(f"WARNING: --keep-labels not found in {label_column}: {missing}")
        rows = rows[labels.isin(keep_labels).to_numpy()]
    if max_per_label:
        rng = np.random.default_rng(seed)
        labels = obs[label_column].astype(str).to_numpy()[rows]
        picked = [rng.choice(rows[labels == lab], size=min(max_per_label, (labels == lab).sum()),
                             replace=False)
                  for lab in np.unique(labels)]
        rows = np.sort(np.concatenate(picked))
    return rows


def read_counts(f, layer, rows):
    """Cells x genes CSR counts for the given rows, read from disk."""
    key = {"raw": "raw/X", "X": "X"}.get(layer, f"layers/{layer}")
    elem = f[key]
    if isinstance(elem, h5py.Dataset):
        return sp.csr_matrix(elem[rows, :])
    matrix = sparse_dataset(elem)
    if matrix.format == "csc":
        return sp.csr_matrix(matrix.to_memory())[rows, :]
    return sp.csr_matrix(matrix[rows])


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--layer", default="counts",
                        help="'raw' (raw/X), 'X', or a key of layers (default counts)")
    parser.add_argument("--panel", help="Xenium gene_panel.json: keep and rename panel genes")
    parser.add_argument("--label-column", default="cell_type")
    parser.add_argument("--keep-labels", default="[]",
                        help="JSON list of --label-column values to keep (default all)")
    parser.add_argument("--max-cells-per-label", type=int, default=0,
                        help="random subsample per label; 0 keeps all")
    parser.add_argument("--obs-columns", nargs="*", default=[])
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    with h5py.File(args.h5ad, "r") as f:
        obs = read_elem(f["obs"])
        var_group = "raw/var" if args.layer == "raw" else "var"
        var_names = read_elem(f[var_group]).index.astype(str)

        rows = select_cells(obs, args.label_column, json.loads(args.keep_labels),
                            args.max_cells_per_label, args.seed)
        counts = read_counts(f, args.layer, rows)
    obs = obs.iloc[rows]

    if args.panel:
        panel = panel_genes(args.panel)
        cols = np.flatnonzero(var_names.isin(list(panel)))
        print(f"{len(cols)} of {len(panel)} panel genes found in {args.h5ad}")
        counts = counts[:, cols]
        var_names = pd.Index([panel[g] for g in var_names[cols]])

    # CELLxGENE raw counts are stored as float; R expects integers
    counts.data = np.rint(counts.data)

    os.makedirs(args.outdir, exist_ok=True)
    with gzip.open(os.path.join(args.outdir, "matrix.mtx.gz"), "wb") as handle:
        scipy.io.mmwrite(handle, counts.T.tocoo(), field="integer")
    pd.Series(var_names).to_csv(os.path.join(args.outdir, "features.tsv.gz"),
                                index=False, header=False)
    pd.Series(obs.index.astype(str)).to_csv(os.path.join(args.outdir, "barcodes.tsv.gz"),
                                            index=False, header=False)
    obs[args.obs_columns].to_csv(os.path.join(args.outdir, "metadata.csv"))

    print(f"{args.h5ad}: {counts.shape[0]} cells x {counts.shape[1]} genes -> {args.outdir}")
    if args.obs_columns:
        print(obs[args.obs_columns[0]].astype(str).value_counts().to_string())


if __name__ == "__main__":
    main()
