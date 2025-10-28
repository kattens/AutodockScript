import os
import subprocess
import pandas as pd
import numpy as np

# --- CONFIG: set your ligand/receptor filenames here ---
LIGAND_FILENAME   = "24826799.pdbqt"  # e.g. in ./output/
RECEPTOR_FILENAME = "5DYK.pdbqt"                               # e.g. in ./NEW_malaria_pdbqt/

# Path to AutoDock Vina executable
vina_executable = r"C:\Program Files (x86)\PyRx\vina.exe"

# --- Paths ---
base_dir = os.path.dirname(os.path.abspath(__file__))
output_folder = os.path.join(base_dir, "final_docking_output")
log_folder = os.path.join(base_dir, "docking_logs")
lig_folder = os.path.join(base_dir, "output")
prot_folder = os.path.join(base_dir, "NEW_malaria_pdbqt")
csv_folder = os.path.join(base_dir, "top10_affinities")

# Ensure dirs exist
for folder in [output_folder, log_folder, csv_folder]:
    os.makedirs(folder, exist_ok=True)


def estimate_docking_box_from_pdbqt(pdbqt_path, buffer=5.0):
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
                # Skip malformed coordinate lines
                continue

    if not xs:
        print("Warning: no heavy-atom coords parsed; using (0,0,0) and 25 Å.")
        return (0.0, 0.0, 0.0), (25.0, 25.0, 25.0)

    coords = np.array(list(zip(xs, ys, zs)))
    center = coords.mean(axis=0)
    minc = coords.min(axis=0)
    maxc = coords.max(axis=0)
    size = (maxc - minc) + buffer
    return (tuple(center), tuple(size))

def run_vina(size, center, log_path, out_path, receptor, ligand):
    cmd = [
        vina_executable,
        "--receptor", receptor,
        "--ligand", ligand,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--out", out_path,
        "--log", log_path,
        "--cpu", "1",
    ]
    # show the exact command (quoted) for easy copy-paste
    pretty = " ".join(f'"{c}"' if " " in c else c for c in cmd)
    print("\nRunning AutoDock Vina with:\n", pretty, "\n")

    # capture Vina stdout+stderr
    proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(proc.stdout)  # show everything Vina said

    if proc.returncode != 0:
        raise RuntimeError(f"Vina failed (exit {proc.returncode}). See output above.")


def parse_affinities(log_path):
    """
    Read top-poses table from Vina log. Returns up to top-10 affinities (kcal/mol).
    """
    affinities = []
    if not os.path.exists(log_path):
        return affinities
    with open(log_path, "r") as f:
        for line in f:
            # Vina table lines begin with rank number 1..N
            s = line.strip()
            if not s:
                continue
            if s.split()[0].isdigit():
                parts = s.split()
                if len(parts) >= 2:
                    try:
                        affinities.append(float(parts[1]))
                    except ValueError:
                        pass
    return affinities[:10]


def save_to_csv(affinities, prefix, csv_folder):
    df = pd.DataFrame({
        "Rank": list(range(1, len(affinities) + 1)),
        "Binding_Affinity_kcal/mol": affinities
    })
    csv_file = os.path.join(csv_folder, f"{prefix}.csv")
    df.to_csv(csv_file, index=False)
    print(f"Top 10 affinities saved to {csv_file}")


def run_docking_with_filenames(ligand_filename, receptor_filename):
    # Build file paths
    ligand_path  = os.path.join(lig_folder, ligand_filename)
    receptor_path = os.path.join(prot_folder, receptor_filename)

    # Basic checks
    if not os.path.exists(ligand_path):
        raise FileNotFoundError(f"Ligand not found: {ligand_path}")
    if not os.path.exists(receptor_path):
        raise FileNotFoundError(f"Receptor not found: {receptor_path}")

    # Prefix for outputs
    lig_base = os.path.splitext(os.path.basename(ligand_path))[0]
    rec_base = os.path.splitext(os.path.basename(receptor_path))[0]
    prefix = f"{lig_base}_vs_{rec_base}"

    # Estimate box (center from receptor), use fixed size cube if desired
    center, _ = estimate_docking_box_from_pdbqt(receptor_path)
    size = (25.0, 25.0, 25.0)

    print(f"Docking box center: {center}, size: {size}")

    # Output paths
    out_pdbqt = os.path.join(output_folder, f"{prefix}.pdbqt")
    log_path  = os.path.join(log_folder,   f"{prefix}.log")

    # Run Vina
    run_vina(size, center, log_path, out_pdbqt, receptor_path, ligand_path)

    # Parse results
    affinities = parse_affinities(log_path)
    if affinities:
        save_to_csv(affinities, prefix, csv_folder)
    else:
        print("No binding affinities parsed. Check the Vina log:", log_path)


if __name__ == "__main__":
    run_docking_with_filenames(LIGAND_FILENAME, RECEPTOR_FILENAME)
