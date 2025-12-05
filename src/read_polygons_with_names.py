#!/usr/bin/env python3
"""
Utility to read a Xenium/coordinates.csv (skiprows=2), group points by the 'Selection'
column, optionally save per-selection CSVs, and return a dict mapping selection names to
Shapely geometries (Polygon, LineString, or raw point).

Usage example:
    from scripts.read_polygons_with_names import read_polygons_with_names
    polygons = read_polygons_with_names('data/coordinates.csv', pixel_size=0.2125, out_dir='out_csvs')

The function will raise a ValueError if expected columns are missing.
"""

from __future__ import annotations

import os
import re
from typing import Dict, Optional, Any

import pandas as pd

def clean_name(s: str) -> str:
    """Normalize a selection name to a safe, lowercase identifier.

    - remove any character not in a-z, 0-9, underscore or space
    - strip surrounding whitespace
    - remove internal spaces
    """
    s = str(s).lower()
    # remove any char that's not a-z, 0-9, underscore or space
    s = re.sub(r"[^a-z0-9_ ]", "", s)
    s = s.strip()
    # finally remove spaces
    s = s.replace(" ", "")
    return s


def _sanitize_filename(name: str) -> str:
    """Make a filesystem-safe filename from the selection name."""
    # Replace whitespace and any non safe chars with underscore
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(name))


def read_polygons_with_names(
    coordinates_csv: str,
    pixel_size: float = 1.0,
    out_dir: Optional[str] = None,
    save_csvs: bool = True,
) -> Dict[str, Any]:
    """
    Read a Xenium coordinates.csv and return a mapping from selection name to geometry.

    Inputs:
    - coordinates_csv: path to CSV file (file must contain columns 'Selection', 'X', 'Y').
      The CSV is read with pd.read_csv(..., skiprows=2) to accommodate the Xenium header.
    - pixel_size: divisor to convert pixel coordinates to physical units (default 1.0).
    - out_dir: optional directory where per-selection CSVs will be written; created if needed.
    - save_csvs: whether to save per-selection CSV files (default True).

    Returns:
    - dict mapping selection name -> geometry. Geometry will be either:
      - shapely.geometry.Polygon if >= 3 points
      - shapely.geometry.LineString if exactly 2 points
      - shapely.geometry.Point if exactly 1 point
      - None for empty selections

    Raises:
    - ValueError if required columns are missing.
    - If Shapely is not installed, the function will still work but will return
      plain Python lists/tuples of coordinates instead of Shapely geometry objects.
    """
    df = pd.read_csv(coordinates_csv, skiprows=2)

    # # Normalize the Selection column so names are cleaned and consistent
    # if "Selection" in df.columns:
    #     df["Selection"] = clean_name(df["Selection"].astype(str))

    # Import shapely lazily so other scripts can import this module even when
    # shapely isn't installed. If shapely is missing, we'll return raw coords.
    try:
        from shapely.geometry import Polygon, LineString, Point  # type: ignore

        _shapely_available = True
    except Exception:
        Polygon = LineString = Point = None  # type: ignore
        _shapely_available = False

    required = {"Selection", "X", "Y"}
    if not required.issubset(set(df.columns)):
        raise ValueError(f"coordinates CSV must contain columns: {', '.join(required)}")

    if out_dir is not None and save_csvs:
        os.makedirs(out_dir, exist_ok=True)

    result: Dict[str, Any] = {}

    for sel, group in df.groupby("Selection"):
        # preserve the order in file; convert to float and apply pixel_size
        pts = group[["X", "Y"]].to_numpy(dtype=float) / float(pixel_size)

        geom = None
        if len(pts) >= 3:
            if _shapely_available:
                geom = Polygon(pts)
                # Try to ensure valid polygon; buffer(0) can fix some invalid rings
                if not geom.is_valid:
                    fixed = geom.buffer(0)
                    if fixed.is_valid:
                        geom = fixed
                    else:
                        # fallback to convex hull
                        geom = geom.convex_hull
            else:
                # return plain coordinate list when shapely isn't present
                geom = pts.tolist()
        elif len(pts) == 2:
            if _shapely_available:
                geom = LineString(pts)
            else:
                geom = pts.tolist()
        elif len(pts) == 1:
            if _shapely_available:
                geom = Point(pts[0][0], pts[0][1])
            else:
                geom = (float(pts[0][0]), float(pts[0][1]))
        else:
            geom = None

        result[clean_name(sel)] = geom

        if save_csvs:
            fname = _sanitize_filename(sel) + ".csv"
            if out_dir:
                path = os.path.join(out_dir, fname)
            else:
                path = fname
            group.to_csv(path, index=False)

    return result


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description="Read Xenium coordinates.csv, save per-selection CSVs, and print a summary."
    )
    p.add_argument("coordinates_csv", help="Path to coordinates.csv (skiprows=2 will be used)")
    p.add_argument("--pixel-size", type=float, default=1.0, help="Pixel size to divide X/Y by")
    p.add_argument("--out-dir", default=None, help="Output directory to store per-selection CSVs")
    p.add_argument("--no-save", dest="save", action="store_false", help="Don't save per-selection CSVs")
    args = p.parse_args()

    polygons = read_polygons_with_names(
        args.coordinates_csv, pixel_size=args.pixel_size, out_dir=args.out_dir, save_csvs=args.save
    )

    # Print a short summary
    print(f"Processed {len(polygons)} selections from {args.coordinates_csv}")
    for name, geom in polygons.items():
        t = type(geom).__name__ if geom is not None else "None"
        print(f" - {name}: {t}")
