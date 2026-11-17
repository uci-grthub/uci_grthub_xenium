#!/usr/bin/env python3
"""
Add transferred cell types (Tangram or Seurat) to a per-sample or merged
AnnData's obs, in place (rules xenium_transfer_obs, xenium_transfer_obs_merged).

Reads --predictions, the merge group's predictions CSV from src/run_tangram.py
(--method tangram) or src/run_seurat_transfer.R (--method seurat_sct or
seurat), takes the object's cells and writes two obs columns (--column, e.g.
tangram_cell_type or seurat_celltype):

  <column>        label; score below --min-score is "unassigned", cells the
                  method skipped (< 50 transcripts) "low_counts"
  <column>_score  prediction score (NaN for low_counts)

Other obs columns of the same method family (tangram_* or seurat_*, e.g. from
an earlier run under another name) are dropped.

Only the obs group of the .h5ad is rewritten (h5py), not the matrix. The
file's mtime is restored afterwards so rules downstream of xenium_anndata do
not rerun. Rerunning xenium_anndata drops the columns; the --done flag is then
older than the .h5ad and Snakemake adds them again.

Usage:
  python src/add_transfer_obs.py --h5ad anndata/<sample>.h5ad --method seurat_sct \
      --predictions seurat_transfer/<group>_predictions_sct.csv \
      --column seurat_celltype --label cell_type_015 --min-score 0.5 \
      --done anndata/seurat_sct_obs/<sample>.done
"""
import argparse
import os

import h5py
import pandas as pd
from anndata.io import read_elem, write_elem

from plot_celltype_methods_spatial import assign


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ad", required=True, help="per-sample or merged AnnData (changed in place)")
    parser.add_argument("--method", required=True, choices=["tangram", "seurat_sct", "seurat"],
                        help="column prefix of the CSV")
    parser.add_argument("--predictions", required=True, help="predictions CSV of the merge group")
    parser.add_argument("--label", default="cell_type_015", help="label column of the CSV")
    parser.add_argument("--column", required=True, help="obs column to write")
    parser.add_argument("--min-score", type=float, default=0.5)
    parser.add_argument("--done", help="flag file touched when finished")
    args = parser.parse_args()

    key = f"{args.method}_{args.label}"
    family = args.method.split("_")[0] + "_"
    pred = pd.read_csv(args.predictions, index_col="cell_id", usecols=["cell_id", key, f"{key}_score"])
    stat = os.stat(args.h5ad)
    with h5py.File(args.h5ad, "r+") as f:
        obs = read_elem(f["obs"])
        missing = ~obs.index.isin(pred.index)
        if missing.any():
            raise ValueError(f"{args.predictions}: {missing.sum()} of {len(obs)} cells "
                             f"of {args.h5ad} have no prediction")
        pred = pred.reindex(obs.index)
        obs = obs.drop(columns=[c for c in obs.columns if c.startswith(family)])
        obs[args.column] = pd.Categorical(assign(pred[key], pred[f"{key}_score"], args.min_score))
        obs[f"{args.column}_score"] = pred[f"{key}_score"].to_numpy()
        del f["obs"]
        write_elem(f, "obs", obs)
    os.utime(args.h5ad, ns=(stat.st_atime_ns, stat.st_mtime_ns))

    counts = obs[args.column].value_counts()
    print(f"{args.h5ad}: {args.column} " + ", ".join(f"{k} {v}" for k, v in counts.items()))
    if args.done:
        os.makedirs(os.path.dirname(os.path.abspath(args.done)), exist_ok=True)
        open(args.done, "w").close()


if __name__ == "__main__":
    main()
