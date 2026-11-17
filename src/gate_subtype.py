#!/usr/bin/env python3
"""
Gate a subtype label transfer on a coarse one (src/run_tangram.py,
src/run_seurat_transfer.R predictions CSVs).

A subtype reference with only tumor-lineage labels would force every cell
into one of them, so only the cells the coarse run called tumor lineage are
mapped with it:

  cells    cell ids whose coarse label is one of --labels with score >=
           --min-score, one per line (input to --cells of run_tangram.py and
           h5ad_to_mtx.py)
  combine  full-slide predictions CSV: the gated cells from the subtype run,
           every other cell copied from the coarse run (label, score,
           total_counts; prob_ columns empty)

Usage:
  python src/gate_subtype.py cells --coarse <group>_predictions.csv \
      --labels Tumor Melanocyte Schwann_cell --min-score 0.5 --out <cells.txt>
  python src/gate_subtype.py combine --coarse <group>_predictions.csv \
      --subtype <group>_subtype_gated_predictions.csv --cells <cells.txt> \
      --out <group>_subtype_predictions.csv
"""
import argparse
import os

import pandas as pd


def label_key(columns, label):
    keys = [c for c in columns if c.endswith(f"_{label}")]
    if len(keys) != 1:
        raise ValueError(f"expected one column ending in _{label}, found {keys}")
    return keys[0]


def gate_cells(args):
    header = pd.read_csv(args.coarse, nrows=0).columns
    key = label_key(header, args.label)
    pred = pd.read_csv(args.coarse, usecols=["cell_id", key, f"{key}_score"])
    keep = pred[key].isin(args.labels) & (pred[f"{key}_score"] >= args.min_score)
    pred.loc[keep, "cell_id"].to_csv(args.out, index=False, header=False)
    print(f"{keep.sum()} of {len(pred)} cells with {key} in {args.labels} "
          f"(score >= {args.min_score:g}) -> {args.out}")
    print(pred.loc[keep, key].value_counts().to_string())


def combine(args):
    cells = pd.read_csv(args.cells, header=None).iloc[:, 0]
    sub = pd.read_csv(args.subtype, index_col="cell_id")
    key = label_key(sub.columns, args.label)
    header = pd.read_csv(args.coarse, nrows=0).columns
    if label_key(header, args.label) != key:
        raise ValueError(f"{args.coarse} and {args.subtype} are not from the same method")
    coarse = pd.read_csv(args.coarse, index_col="cell_id",
                         usecols=["cell_id", key, f"{key}_score", "total_counts"])
    missing = ~cells.isin(sub.index)
    if missing.any():
        raise ValueError(f"{missing.sum()} gated cells have no prediction in {args.subtype}")

    out = coarse.reindex(columns=sub.columns)
    out.loc[cells, :] = sub.loc[cells, :].to_numpy()
    out[key] = out[key].astype(str)
    out.to_csv(args.out)
    print(f"{len(cells)} gated of {len(out)} cells -> {args.out}")
    print(out.loc[cells, key].value_counts().to_string())


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    cells = sub.add_parser("cells", help="gated cell ids from the coarse run")
    cells.add_argument("--labels", nargs="+", default=["Tumor", "Melanocyte", "Schwann_cell"])
    cells.add_argument("--min-score", type=float, default=0.5)
    comb = sub.add_parser("combine", help="gated subtype + coarse predictions")
    comb.add_argument("--subtype", required=True, help="predictions of the gated cells")
    comb.add_argument("--cells", required=True, help="output of the cells command")
    for p in (cells, comb):
        p.add_argument("--coarse", required=True, help="coarse predictions CSV")
        p.add_argument("--label", default="cell_type_015", help="reference label column")
        p.add_argument("--out", required=True)
    args = parser.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    gate_cells(args) if args.command == "cells" else combine(args)


if __name__ == "__main__":
    main()
