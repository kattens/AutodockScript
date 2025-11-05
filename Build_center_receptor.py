#!/usr/bin/env python3
# Build_center_receptor.py
import json
import math
from pathlib import Path

import pandas as pd
import numpy as np
from pandas.errors import EmptyDataError

# ---------- config ----------
BASE_DIR = Path(".").resolve()
CSV_DIR = BASE_DIR / "contacts_csv"      # folder with per-PDB contact CSVs
PDB_DIR = BASE_DIR / "malaria"           # folder with PDBs; keys are taken from here
OUTPUT_JSON = BASE_DIR / "Center_boxes.json"

PAD = 3.0
MIN_SIZE = 16.0
FORCE_CUBE = True
# ----------------------------

def parse_coord(val):
    """Accept '(x, y, z)' strings or tuples/lists; return (x, y, z) floats."""
    if val is None or (isinstance(val, float) and math.isnan(val)):
        raise ValueError("coord is NaN/None")

    if isinstance(val, (tuple, list)) and len(val) == 3:
        x, y, z = val
        return float(x), float(y), float(z)

    if not isinstance(val, str):
        raise ValueError(f"not a string/tuple/list coord: {val!r}")

    v = val.strip()
    if (v.startswith("(") and v.endswith(")")) or (v.startswith("[") and v.endswith("]")):
        v = v[1:-1]

    parts = [p.strip() for p in v.split(",")]
    if len(parts) != 3:
        raise ValueError(f"cannot parse coord: {val!r}")

    x, y, z = (float(p) for p in parts)
    return x, y, z

def process_csv(csv_path: Path):
    """Return dict with center/size/n_points from one CSV, or None if invalid/empty."""
    # quick skip for zero-byte files
    try:
        if csv_path.stat().st_size == 0:
            print(f"[warn] {csv_path.name}: empty file (0 bytes), skipping")
            return None
    except FileNotFoundError:
        print(f"[warn] missing file: {csv_path}")
        return None

    # safe read
    try:
        df = pd.read_csv(csv_path)
    except EmptyDataError:
        print(f"[warn] {csv_path.name}: EmptyDataError (no columns), skipping")
        return None

    if df is None or df.shape[1] == 0:
        print(f"[warn] {csv_path.name}: no columns after read, skipping")
        return None

    # Accept either column name; contact script emits 'prot_coord'
    coord_col = None
    for candidate in ("prot_coord", "protein_coord"):
        if candidate in df.columns:
            coord_col = candidate
            break
    if coord_col is None:
        print(f"[warn] {csv_path.name}: no prot_coord/protein_coord column, skipping")
        return None

    # extract coords
    coords = []
    for val in df[coord_col]:
        try:
            coords.append(parse_coord(val))
        except Exception:
            # ignore malformed rows
            continue

    if not coords:
        print(f"[warn] {csv_path.name}: no valid coords, skipping")
        return None

    # Deduplicate identical protein-atom coords to avoid overweighting
    coords = list({(round(x, 3), round(y, 3), round(z, 3)) for x, y, z in coords})
    arr = np.array(coords, dtype=float)  # (N, 3)

    center = arr.mean(axis=0)
    minc = arr.min(axis=0)
    maxc = arr.max(axis=0)
    size = (maxc - minc) + 2 * PAD
    size = np.maximum(size, MIN_SIZE)

    if FORCE_CUBE:
        m = float(size.max())
        size = np.array([m, m, m], dtype=float)

    return {
        "center": [float(center[0]), float(center[1]), float(center[2])],
        "size":   [float(size[0]),   float(size[1]),   float(size[2])],
        "n_points": int(arr.shape[0]),
    }

def main():
    results = {}

    # loop over PDBs so JSON keys match docking names
    for pdb_path in sorted(PDB_DIR.glob("*.pdb")):
        pdb_stem = pdb_path.stem

        # primary CSV name
        csv_path = CSV_DIR / f"{pdb_stem}.csv"
        if not csv_path.exists():
            # try alternate naming with/without _malaria
            alt = pdb_stem.removesuffix("_malaria") if pdb_stem.endswith("_malaria") else f"{pdb_stem}_malaria"
            alt_csv = CSV_DIR / f"{alt}.csv"
            if alt_csv.exists():
                csv_path = alt_csv
            else:
                print(f"[warn] CSV for {pdb_stem} not found (tried {csv_path.name}, {alt_csv.name})")
                continue

        data = process_csv(csv_path)
        if data is None:
            continue

        results[pdb_stem] = data
        print(f"[ok] {pdb_stem}: center={data['center']} size={data['size']} (N={data['n_points']})")

    with OUTPUT_JSON.open("w") as f:
        json.dump(results, f, indent=2)

    print(f"\nDone. Wrote {OUTPUT_JSON} with {len(results)} entries.")

if __name__ == "__main__":
    main()
