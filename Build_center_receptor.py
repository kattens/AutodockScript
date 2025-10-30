
import json
from pathlib import Path

import pandas as pd
import numpy as np

# config
BASE_DIR = Path(".").resolve()
CSV_DIR = BASE_DIR / "contact_output"     
PDB_DIR = BASE_DIR / "malaria"       
OUTPUT_JSON = BASE_DIR / "Center_boxes.json"

PAD = 3.0     
MIN_SIZE = 16.0  
FORCE_CUBE = True


def parse_coord(val: str):

    if not isinstance(val, str):
        raise ValueError(f"not a string coord: {val}")

    v = val.strip()
    # remove leading/trailing () or []
    if (v.startswith("(") and v.endswith(")")) or (v.startswith("[") and v.endswith("]")):
        v = v[1:-1]

    parts = v.split(",")
    if len(parts) != 3:
        raise ValueError(f"cannot parse coord: {val}")

    x, y, z = (float(p.strip()) for p in parts)
    return x, y, z


def process_csv(csv_path: Path):
    """Return dict with center/size/n_points from one CSV, or None if invalid."""
    df = pd.read_csv(csv_path)
    if "protein_coord" not in df.columns:
        print(f"[warn] {csv_path.name}: no 'protein_coord' column, skipping")
        return None

    coords = []
    for val in df["protein_coord"]:
        try:
            xyz = parse_coord(val)
            coords.append(xyz)
        except Exception as e:
            print(f"[warn] {csv_path.name}: bad coord {val!r} -> {e}")
            continue

    if not coords:
        print(f"[warn] {csv_path.name}: no valid coords, skipping")
        return None

    arr = np.array(coords, dtype=float)  # shape (N, 3)
    center = arr.mean(axis=0)

    # bounding box
    minc = arr.min(axis=0)
    maxc = arr.max(axis=0)
    size = (maxc - minc) + 2 * PAD

    # ensure min size
    size = np.maximum(size, MIN_SIZE)

    if FORCE_CUBE:
        m = float(size.max())
        size = np.array([m, m, m], dtype=float)

    return {
        "center": [float(center[0]), float(center[1]), float(center[2])],
        "size": [float(size[0]), float(size[1]), float(size[2])],
        "n_points": int(arr.shape[0]),
    }


def main():
    results = {}

    # loop over ALL pdbs
    for pdb_path in sorted(PDB_DIR.glob("*.pdb")):
        pdb_stem = pdb_path.stem  
        if pdb_stem.endswith("_malaria"):
            csv_stem = pdb_stem.removesuffix("_malaria")
        else:
            csv_stem = pdb_stem

        csv_path = CSV_DIR / f"{csv_stem}.csv"
        if not csv_path.exists():
            print(f"[warn] CSV for {pdb_stem} not found -> expected {csv_path.name}")
            continue

        data = process_csv(csv_path)
        if data is None:
            continue

        # store using the PDB name (with _malaria) so it matches docking files
        results[pdb_stem] = data
        print(f"[ok] {pdb_stem}: center={data['center']} size={data['size']} (N={data['n_points']})")

    # write once
    with OUTPUT_JSON.open("w") as f:
        json.dump(results, f, indent=2)

    print(f"\nDone. Wrote {OUTPUT_JSON} with {len(results)} entries.")


if __name__ == "__main__":
    main()
