#!/usr/bin/env python3
"""
Add Tangram cell types to a per-sample or merged AnnData's obs, in place
(rules xenium_tangram_obs, xenium_tangram_obs_merged).

Reads --tangram, the merge group's predictions CSV from src/run_tangram.py,
takes the object's cells and writes two obs columns (--column, default
tangram_cell_type):

  <column>        label; score below --min-score is "unassigned", cells
                  Tangram skipped (< 50 transcripts) "low_counts"
  <column>_score  Tangram score (NaN for low_counts)

Other tangram_* obs columns (e.g. from an earlier run under another name) are
dropped.

Only the obs group of the .h5ad is rewritten (h5py), not the matrix. The
file's mtime is restored afterwards so rules downstream of xenium_anndata do
not rerun. Rerunning xenium_anndata drops the columns; the --done flag is then
older than the .h5ad and Snakemake adds them again.

Usage:
  python src/add_tangram_obs.py --h5ad anndata/<sample>.h5ad \
      --tangram tangram/<group>_predictions.csv --label cell_type_015 \
      --min-score 0.5 --done anndata/tangram_obs/<sample>.done
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
    parser.add_argument("--tangram", required=True, help="Tangram predictions CSV of the merge group")
    parser.add_argument("--label", default="cell_type_015", help="label column of the CSV")
    parser.add_argument("--column", default="tangram_cell_type", help="obs column to write")
    parser.add_argument("--min-score", type=float, default=0.5)
    parser.add_argument("--done", help="flag file touched when finished")
    args = parser.parse_args()

    key = f"tangram_{args.label}"
    pred = pd.read_csv(args.tangram, index_col="cell_id", usecols=["cell_id", key, f"{key}_score"])
    stat = os.stat(args.h5ad)
    with h5py.File(args.h5ad, "r+") as f:
        obs = read_elem(f["obs"])
        missing = ~obs.index.isin(pred.index)
        if missing.any():
            raise ValueError(f"{args.tangram}: {missing.sum()} of {len(obs)} cells "
                             f"of {args.h5ad} have no prediction")
        pred = pred.reindex(obs.index)
        obs = obs.drop(columns=[c for c in obs.columns if c.startswith("tangram_")])
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
