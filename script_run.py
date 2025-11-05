import os
import sys
import json
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd

vina_executable = r"C:\Program Files (x86)\PyRx\vina.exe" #on windows
#vina_executable = r"/usr/local/bin/vina"   # Mac

base_dir = Path(__file__).resolve().parent
output_folder = base_dir / "docking_output"
log_folder = base_dir / "docking_logs"
lig_folder = base_dir / "ligand_pdbqt"
prot_folder = base_dir / "mraw_alaria_pdbqt"
csv_folder = base_dir / "top10_affinities"

binding_centers_path = base_dir / "Center_boxes.json"
# ------------------------------------------

# Load precomputed centers if present
if binding_centers_path.exists():
    with open(binding_centers_path, "r") as f:
        BINDING_CENTERS = json.load(f)
    print(f"[info] loaded binding centers from {binding_centers_path} ({len(BINDING_CENTERS)} entries)")
else:
    BINDING_CENTERS = {}
    print("[info] Center_boxes.json not found, will use receptor-based auto box.")

def estimate_docking_box_from_pdbqt(pdbqt_path: Path, buffer: float = 5.0):
    """
    Estimate center/size from receptor PDBQT (skips hydrogens).
    Returns (center_xyz, size_xyz) where center/size are tuples of 3 floats.
    """
    xs, ys, zs = [], [], []
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
    size = (25.0, 25.0, 25.0)
    print(f"[auto-box] {rec_stem}: center={center}, size={size}")
    return center, size

def run_vina(size, center, log_path: Path, out_path: Path, receptor: Path, ligand: Path):
    """
    Run vina 1.2.x: no --log option. Capture stdout and write it to log_path.
    """
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
    # Save stdout as our log file (since --log is not supported)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w") as lf:
        lf.write(proc.stdout or "")
    # Also echo to console
    print(proc.stdout)

    if proc.returncode != 0:
        raise RuntimeError(f"Vina failed (exit {proc.returncode}). See log: {log_path}")

def parse_affinities(log_path: Path):
    """
    Parse top-poses table from a Vina 1.2.x log (captured stdout).
    Returns up to top-10 affinities (kcal/mol).
    """
    affinities = []
    if not log_path.exists():
        return affinities
    with open(log_path, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            first = s.split()[0]
            if first.isdigit():
                parts = s.split()
                if len(parts) >= 2:
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

def run_docking_pair(ligand_path: Path, receptor_path: Path):
    lig_base = ligand_path.stem
    rec_base = receptor_path.stem
    prefix = f"{lig_base}_vs_{rec_base}"

    center, size = build_box_for_receptor(receptor_path)

    out_pdbqt = output_folder / f"{prefix}.pdbqt"
    log_path  = log_folder   / f"{prefix}.log"

    # Skip if output already exists
    if out_pdbqt.exists():
        print(f"[skip] {out_pdbqt.name} already exists, skipping.")
        return

    run_vina(size=size, center=center, log_path=log_path, out_path=out_pdbqt,
             receptor=receptor_path, ligand=ligand_path)

    affinities = parse_affinities(log_path)
    if affinities:
        save_to_csv(affinities, prefix, csv_folder)
    else:
        print(f"[warn] No affinities parsed. Check log: {log_path}")

if __name__ == "__main__":
    # Ensure output dirs exist
    output_folder.mkdir(parents=True, exist_ok=True)
    log_folder.mkdir(parents=True, exist_ok=True)
    csv_folder.mkdir(parents=True, exist_ok=True)

    # Build stem -> path maps
    lig_map = {p.stem: p for p in lig_folder.glob("*.pdbqt")}
    rec_map = {p.stem: p for p in prot_folder.glob("*.pdbqt")}

    print(f"[inventory] ligands: {len(lig_map)} | receptors: {len(rec_map)}")

    # Match by identical stem in both folders
    common = sorted(set(lig_map.keys()) & set(rec_map.keys()))
    print(f"[match] matched pairs by stem: {len(common)}")

    if not common:
        print("[fatal] No matched pairs found. Ensure both folders contain .pdbqt files with identical names (stems).")
        sys.exit(1)

    processed = 0
    errors = 0

    for stem in common:
        lig_path = lig_map[stem]
        rec_path = rec_map[stem]
        try:
            run_docking_pair(lig_path, rec_path)
            processed += 1
        except Exception as e:
            print(f"[error] {stem}: {e}")
            errors += 1

    print(f"\nDone. Processed: {processed} | Errors: {errors}\n")
