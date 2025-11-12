import os
import sys
import subprocess
import pandas as pd
import numpy as np
import json
import random
from pathlib import Path

# ---------------- Config ----------------
#vina_executable = r"/usr/local/bin/vina"  # on Mac / Linux
vina_executable = r"C:\Program Files (x86)\PyRx\vina.exe"  # on Windows
RANDOM_SEED = 42              # set None for non-deterministic picks
NUM_RECEPTORS = 20            # pick up to this many random receptors
LIGANDS_PER_RECEPTOR = 2      # random ligands per receptor
AUTO_BOX_SIZE = (25.0, 25.0, 25.0)
HEAVY_ATOM_BUFFER = 5.0
# ---------------------------------------

base_dir = Path(__file__).resolve().parent
lig_folder = base_dir / "ligand_pdbqt"
prot_folder = base_dir / "malaria_pdbqt"
results_root = base_dir / "results"
results_root.mkdir(parents=True, exist_ok=True)

binding_centers_path = base_dir / "Center_boxes.json"

# Load precomputed centers if present
if binding_centers_path.exists():
    with open(binding_centers_path, "r") as f:
        BINDING_CENTERS = json.load(f)
    print(f"[info] loaded binding centers from {binding_centers_path} ({len(BINDING_CENTERS)} entries)")
else:
    BINDING_CENTERS = {}
    print("[info] Center_boxes.json not found, will use receptor-based auto box.")


def estimate_docking_box_from_pdbqt(pdbqt_path: Path, buffer: float = HEAVY_ATOM_BUFFER):
    xs, ys, zs = [], [], []
    try:
        with open(pdbqt_path, "r") as f:
            for line in f:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                atom_name = line[12:16].strip()
                if atom_name.startswith("H"):
                    continue
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    xs.append(x); ys.append(y); zs.append(z)
                except ValueError:
                    continue
    except FileNotFoundError:
        print(f"[warn] receptor file not found for auto box: {pdbqt_path}")
        return (0.0, 0.0, 0.0), (25.0, 25.0, 25.0)

    if not xs:
        print(f"[warn] no heavy-atom coords parsed; using (0,0,0) and 25 Å for {pdbqt_path.name}.")
        return (0.0, 0.0, 0.0), (25.0, 25.0, 25.0)

    coords = np.array(list(zip(xs, ys, zs)))
    center = coords.mean(axis=0)
    minc = coords.min(axis=0)
    maxc = coords.max(axis=0)
    size = (maxc - minc) + buffer
    return (tuple(center), tuple(size))


def build_box_for_receptor(receptor_path: Path):
    rec_stem = receptor_path.stem
    if rec_stem in BINDING_CENTERS:
        data = BINDING_CENTERS[rec_stem]
        center = tuple(data["center"])
        size   = tuple(data["size"])
        print(f"[site-box] {rec_stem}: center={center}, size={size}")
        return center, size
    center, _ = estimate_docking_box_from_pdbqt(receptor_path)
    size = AUTO_BOX_SIZE
    print(f"[auto-box] {rec_stem}: center={center}, size={size}")
    return center, size


def run_vina(size, center, log_path: Path, out_path: Path, receptor: Path, ligand: Path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(vina_executable),
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--out", str(out_path),
        "--cpu", "1",
    ]

    pretty = " ".join(f'"{c}"' if " " in c else c for c in cmd)
    print("\n[run] AutoDock Vina:\n", pretty, "\n")

    proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    with open(log_path, "w") as lf:
        lf.write(proc.stdout or "")

    if proc.stdout:
        print(proc.stdout)

    if proc.returncode != 0:
        raise RuntimeError(f"Vina failed (exit {proc.returncode}). See log: {log_path}")


def parse_affinities(log_path: Path):
    affinities = []
    if not log_path.exists():
        return affinities
    with open(log_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if parts and parts[0].isdigit() and len(parts) >= 2:
                try:
                    affinities.append(float(parts[1]))
                except ValueError:
                    pass
    return affinities[:10]


def save_to_csv(affinities, prefix: str, csv_dir: Path):
    df = pd.DataFrame({
        "Rank": list(range(1, len(affinities) + 1)),
        "Binding_Affinity_kcal/mol": affinities
    })
    csv_dir.mkdir(parents=True, exist_ok=True)
    csv_file = csv_dir / f"{prefix}.csv"
    df.to_csv(csv_file, index=False)
    print(f"[ok] saved CSV: {csv_file}")


def run_docking_pair(ligand_path: Path, receptor_path: Path, out_dir: Path, log_dir: Path, csv_dir: Path):
    lig_base = ligand_path.stem
    rec_base = receptor_path.stem
    prefix = f"{lig_base}_vs_{rec_base}"

    out_pdbqt = out_dir / f"{prefix}.pdbqt"
    log_path  = log_dir / f"{prefix}.log"

    # Skip if output already exists
    if out_pdbqt.exists():
        print(f"[skip] exists: {out_pdbqt.name}")
        return

    center, size = build_box_for_receptor(receptor_path)

    run_vina(size=size, center=center, log_path=log_path, out_path=out_pdbqt,
             receptor=receptor_path, ligand=ligand_path)

    affinities = parse_affinities(log_path)
    if affinities:
        save_to_csv(affinities, prefix, csv_dir)
    else:
        print(f"[warn] No affinities parsed. Check log: {log_path}")


def pick_random(paths, k):
    if k <= 0:
        return []
    if len(paths) <= k:
        return list(paths)
    return random.sample(paths, k)


def main():
    if RANDOM_SEED is not None:
        random.seed(RANDOM_SEED)

    lig_all = sorted(lig_folder.glob("*.pdbqt"))
    rec_all = sorted(prot_folder.glob("*.pdbqt"))
    print(f"[inventory] ligands: {len(lig_all)} | receptors: {len(rec_all)}")

    if len(lig_all) == 0 or len(rec_all) == 0:
        print("[fatal] Need both ligand_pdbqt and Malaria_Dataset_only with .pdbqt files.")
        sys.exit(1)

    receptors = pick_random(rec_all, NUM_RECEPTORS)
    print(f"[plan] receptors (random {len(receptors)}): {[p.stem for p in receptors]}")

    processed = 0
    errors = 0

    for rec_path in receptors:
        rec_stem = rec_path.stem
        rec_root = results_root / rec_stem
        out_dir  = rec_root / "out"
        log_dir  = rec_root / "logs"
        csv_dir  = rec_root / "csv"
        out_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)
        csv_dir.mkdir(parents=True, exist_ok=True)

        # Randomly pick two *different* ligands
        if len(lig_all) < LIGANDS_PER_RECEPTOR:
            ligands = lig_all
        else:
            ligands = random.sample(lig_all, LIGANDS_PER_RECEPTOR)

        print(f"[plan] {rec_stem}: ligands -> {[p.stem for p in ligands]}")

        for lig_path in ligands:
            pair = f"{lig_path.stem} vs {rec_stem}"
            try:
                run_docking_pair(lig_path, rec_path, out_dir, log_dir, csv_dir)
                processed += 1
            except Exception as e:
                print(f"[error] {pair}: {e}")
                errors += 1

    print(f"\nDone. Processed: {processed} | Errors: {errors}\n")


if __name__ == "__main__":
    main()
