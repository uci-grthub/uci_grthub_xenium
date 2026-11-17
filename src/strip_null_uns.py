#!/usr/bin/env python3
"""
Remove None values from .h5ad files so anndata 0.11 can read them, in place.

anndata >= 0.12 writes a None in uns (e.g. uns/log1p/base, uns/pca/params/mask_var
from scanpy) with encoding-type "null", which anndata 0.11 has no reader for:

  IORegistryError: No read method registered for IOSpec(encoding_type='null', ...)

Every node with that encoding is deleted (h5py), which is what a missing key means
to scanpy anyway. The file's mtime is restored afterwards so Snakemake does not rerun
downstream rules. A file another process has open (HDF5 lock) is skipped and the
exit status is 1. New files are written without them (drop_none in
src/xenium_to_anndata.py).

Usage:
  python src/strip_null_uns.py output/xenium/anndata/*.h5ad output/xenium/anndata/merged*/*.h5ad
"""
import argparse
import os

import h5py


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("h5ads", nargs="+", help="AnnData files (changed in place)")
    args = parser.parse_args()

    locked = []
    for path in args.h5ads:
        stat = os.stat(path)
        try:
            f = h5py.File(path, "r+")
        except BlockingIOError:
            print(f"{path}: open in another process, skipped")
            locked.append(path)
            continue
        with f:
            nulls = []
            f.visititems(lambda name, obj: nulls.append(name)
                         if obj.attrs.get("encoding-type") == "null" else None)
            for name in nulls:
                del f[name]
        os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        print(f"{path}: removed {', '.join(nulls) if nulls else 'nothing'}")
    if locked:
        raise SystemExit(f"{len(locked)} file(s) skipped: {' '.join(locked)}")


if __name__ == "__main__":
    main()
