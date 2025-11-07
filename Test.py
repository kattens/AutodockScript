import os
import subprocess
import pandas as pd
import numpy as np
import json
from pathlib import Path

SELECTED_RECEPTOR_STEMS = [
    "4NCT_2K2_O77239"
]
SELECTED_LIGAND_STEMS = [
    "6VNE_1OB3_2TA","1CET_CLQ_Q27743"
]

#vina_executable = r"C:\Program Files (x86)\PyRx\vina.exe"  # on Windows
vina_executable = r"/usr/local/bin/vina"  # on Mac

base_dir = Path(__file__).resolve().parent
lig_folder = base_dir / "ligand_pdbqt"   
prot_folder = base_dir / "malaria_pdbqt"
results_root = base_dir / "results"
results_root.mkdir(parents=True, exist_ok=True)

binding_centers_path = base_dir / "Center_boxes.json"

# ------------------------------------------
# Load optional precomputed boxes
if binding_centers_path.exists():
    with open(binding_centers_path, "r") as f:
        BINDING_CENTERS = json.load(f)
    print(f"[info] loaded binding centers from {binding_centers_path} ({len(BINDING_CENTERS)} entries)")
else:
    BINDING_CENTERS = {}
    print("[info] binding_centers.json not found, will use receptor-based auto box.")

def estimate_docking_box_from_pdbqt(pdbqt_path: Path, buffer: float = 5.0):
    """
    Robust center/size estimation from PDBQT via fixed-width slicing.
    Skips hydrogens. Returns (center_xyz, size_xyz).
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
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                xs.append(x); ys.append(y); zs.append(z)
            except ValueError:
                continue

    if not xs:
        print(f"[warn] no heavy-atom coords parsed in {pdbqt_path.name}; using (0,0,0) and 25 Å.")
        return (0.0, 0.0, 0.0), (25.0, 25.0, 25.0)

    coords = np.array(list(zip(xs, ys, zs)))
    center = coords.mean(axis=0)
    minc = coords.min(axis=0)
    maxc = coords.max(axis=0)
    size = (maxc - minc) + buffer  # numpy broadcast, then tuple-ify below
    return (tuple(center), tuple(size))

def run_vina(size, center, log_path: Path, out_path: Path, receptor: Path, ligand: Path):
    """
    Vina 1.2.x: there is NO --log flag. Capture stdout and save it as the log.
    """
    # Make sure the directory exists for outputs
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
    # Save stdout to our log file
    with open(log_path, "w") as lf:
        lf.write(proc.stdout or "")
    print(proc.stdout)

    if proc.returncode != 0:
        raise RuntimeError(f"Vina failed (exit {proc.returncode}) for {ligand.name} vs {receptor.name}. See log: {log_path}")

def parse_affinities(log_path: Path):
    """Parse top-poses table from a Vina 1.2.x stdout log. Returns up to top-10 affinities (kcal/mol)."""
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

def save_affinities_csv(affinities, out_csv_path: Path):
    df = pd.DataFrame({
        "Rank": list(range(1, len(affinities) + 1)),
        "Binding_Affinity_kcal/mol": affinities
    })
    out_csv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv_path, index=False)
    print(f"[ok] saved CSV: {out_csv_path}")

def build_box_for_receptor(receptor_path: Path):
    rec_base = receptor_path.stem
    # Prefer precomputed binding site boxes if key matches the receptor stem
    if rec_base in BINDING_CENTERS:
        data = BINDING_CENTERS[rec_base]
        center = tuple(data["center"])
        size   = tuple(data["size"])
        print(f"[site-box] Using CSV-based box for {rec_base}: center={center}, size={size}")
        return center, size
    # Fallback: receptor-based automatic box
    center, _ = estimate_docking_box_from_pdbqt(receptor_path)
    size = (25.0, 25.0, 25.0)
    print(f"[auto-box] Using receptor-based box for {rec_base}: center={center}, size={size}")
    return center, size

def resolve_stems_to_paths(folder: Path, stems):
    """Return existing paths for the given stems (without extension), error on missing."""
    if not stems:
        raise ValueError("No names provided in SELECTED_*_STEMS.")
    available = {p.stem: p for p in folder.glob("*.pdbqt")}
    missing = [s for s in stems if s not in available]
    if missing:
        raise FileNotFoundError(
            f"The following stems were not found as .pdbqt in {folder}:\n  - " +
            "\n  - ".join(missing)
        )
    return [available[s] for s in stems]

def main():
    # 1) Resolve explicit selections
    receptors = resolve_stems_to_paths(prot_folder, SELECTED_RECEPTOR_STEMS)
    ligands   = resolve_stems_to_paths(lig_folder, SELECTED_LIGAND_STEMS)

    if len(ligands) != 2:
        print(f"[warn] You listed {len(ligands)} ligands; this script will dock ALL listed ligands per receptor.")
    print(f"[plan] receptors selected: {len(receptors)} -> {[p.stem for p in receptors]}")
    print(f"[plan] ligands selected:   {len(ligands)} -> {[p.stem for p in ligands]}")

    processed = 0
    errors = 0

    # 2) Run: each selected receptor with each selected ligand
    for rec_path in receptors:
        malaria_name = rec_path.stem
        malaria_out_dir = results_root / malaria_name
        malaria_out_dir.mkdir(parents=True, exist_ok=True)

        # Box once per receptor
        center, size = build_box_for_receptor(rec_path)

        for lig_path in ligands:
            lig_base = lig_path.stem
            rec_base = rec_path.stem
            prefix = f"{lig_base}_vs_{rec_base}"

            out_pdbqt = malaria_out_dir / f"{prefix}.pdbqt"
            log_path  = malaria_out_dir / f"{prefix}.log"
            out_csv   = malaria_out_dir / f"{prefix}.csv"

            try:
                run_vina(size=size, center=center, log_path=log_path, out_path=out_pdbqt,
                         receptor=rec_path, ligand=lig_path)

                affinities = parse_affinities(log_path)
                if affinities:
                    save_affinities_csv(affinities, out_csv)
                else:
                    print(f"[warn] no affinities parsed from {log_path.name}")

                processed += 1
            except Exception as e:
                print(f"[error] {prefix}: {e}")
                errors += 1

    print(f"\nDone. Dockings run: {processed} | Errors: {errors}")
    print(f"Results are under: {results_root}")

if __name__ == "__main__":
    main()
