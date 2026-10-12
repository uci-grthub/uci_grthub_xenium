#!/usr/bin/env python3
"""
Collect the numbers the project report needs from the processed AnnData
objects, so the report itself (src/generate_xenium_report.py) runs without
anndata/h5py.

Reads each .h5ad in backed mode (obs, var and uns only; X stays on disk) and
writes one JSON:

  {"samples": {<sample_id>: {...}}, "merged": {<group>: {...}}}

Per object: cells, genes, Leiden clusters, uns["preprocessing"]; merged
objects also carry cells per sample.

Usage:
  python src/xenium_report_stats.py \
      --h5ads output/xenium/anndata/*.h5ad \
      --merged output/xenium/anndata/merged/*.h5ad \
      --out output/xenium/report_stats.json
"""
import argparse
import json
import os

import anndata as ad
import numpy as np


def jsonable(value):
    """uns values (numpy scalars/arrays, nested dicts) as plain JSON types."""
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def object_stats(path):
    adata = ad.read_h5ad(path, backed="r")
    stats = {
        "path": path,
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "clusters": int(adata.obs["leiden"].nunique()) if "leiden" in adata.obs else None,
        "preprocessing": jsonable(dict(adata.uns.get("preprocessing", {}))),
    }
    if "sample_id" in adata.obs:
        counts = adata.obs["sample_id"].value_counts(sort=False)
        stats["cells_per_sample"] = {str(k): int(v) for k, v in counts.items()}
    adata.file.close()
    return stats


def name(path):
    return os.path.splitext(os.path.basename(path))[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ads", nargs="+", required=True, help="per-sample objects")
    parser.add_argument("--merged", nargs="*", default=[], help="merged objects")
    parser.add_argument("--out", required=True, help="output .json")
    args = parser.parse_args()

    stats = {
        "samples": {name(p): object_stats(p) for p in sorted(args.h5ads)},
        "merged": {name(p): object_stats(p) for p in sorted(args.merged)},
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(stats, fh, indent=2)
    print(f"{len(stats['samples'])} samples, {len(stats['merged'])} merged objects -> {args.out}")


if __name__ == "__main__":
    main()
