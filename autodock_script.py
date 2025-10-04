import os
import sys
import csv
import subprocess
import pandas as pd
import numpy as np
from collections import defaultdict

# ===================== CONFIG =====================
vina_executable = r"C:\Program Files (x86)\PyRx\vina.exe"
base_dir = os.path.dirname(os.path.abspath(__file__))

lig_folder    = os.path.join(base_dir, "ligands_pdbqt")
prot_folder   = os.path.join(base_dir, "malaria_pdbqt")
output_folder = os.path.join(base_dir, "output")
log_folder    = os.path.join(base_dir, "docking_logs")
csv_folder    = os.path.join(base_dir, "top10_affinities")

for d in (output_folder, log_folder, csv_folder):
    os.makedirs(d, exist_ok=True)

csv_path = os.path.join(base_dir, "Data.csv")

# ---- Box size: ~25 Å cube around the co-crystallized ligand ----
FIXED_EDGE = 25.0  # size_x = size_y = size_z

# Common solvents/ions/cofactors to exclude when auto-detecting bound ligand
EXCLUDE_RESN = {
    "HOH","WAT","DOD","TP3","H2O",
    "NA","K","CL","CA","MG","MN","ZN","FE","CU","CD","NI",
    "SO4","PO4","NO3","CO3","IOD","BR","I","F",
    "GOL","PG4","PEG","EDO","MPD","BME","TLA",
    "NAG","NDG","MAN","BMA","GAL","FUC","GLC","GLO","SUC","HEM","HEC","NAP",
}

# ==================================================


def ensure_pdbqt_name(name: str) -> str:
    base = os.path.splitext(str(name).strip())[0]
    return base + ".pdbqt"


def parse_pdbqt_het_groups(pdbqt_path):
    """
    Parse HETATM coordinates from a PDBQT receptor.
    Group atoms by (resname, chain, resseq). Return dict -> list of (x,y,z).
    """
    groups = defaultdict(list)
    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            # PDBQT uses the same fixed columns as PDB for ATOM/HETATM records
            if not line.startswith("HETATM"):
                continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            # residue name, chain, residue sequence (PDB fixed columns)
            resn  = line[17:21].strip().upper()
            chain = line[21:22].strip()
            resi  = line[22:26].strip()
            key = (resn, chain, resi)
            groups[key].append((x, y, z))
    return groups


def pick_ligand_group(groups, target_resn=None):
    """
    Choose which HETATM group to use for centering:
      - If target_resn is provided, pick the largest group with that resname.
      - Else, auto-detect: choose the largest non-excluded resname group.
    Returns (center_xyz tuple, (resn, chain, resi), n_atoms).
    Raises ValueError if not found.
    """
    if not groups:
        raise ValueError("No HETATM groups found in receptor PDBQT.")

    # Normalize for matching
    if target_resn:
        target_resn = str(target_resn).strip().upper()

    best = None
    if target_resn:
        # pick largest group with the requested resname
        candidates = [(k, v) for k, v in groups.items() if k[0] == target_resn]
        if not candidates:
            raise ValueError(f"Requested ligand residue '{target_resn}' not found in HETATM groups.")
        k, pts = max(candidates, key=lambda kv: len(kv[1]))
        coords = np.array(pts, dtype=np.float32)
        center = tuple(coords.mean(axis=0).tolist())
        return center, k, len(pts)

    # Auto-detect: filter out common solvents/ions, pick largest remaining
    filtered = [(k, v) for k, v in groups.items() if k[0] not in EXCLUDE_RESN]
    if not filtered:
        # fallback: pick absolute largest group even if excluded set (warn-ish)
        k, pts = max(groups.items(), key=lambda kv: len(kv[1]))
        coords = np.array(pts, dtype=np.float32)
        center = tuple(coords.mean(axis=0).tolist())
        return center, k, len(pts)

    k, pts = max(filtered, key=lambda kv: len(kv[1]))
    coords = np.array(pts, dtype=np.float32)
    center = tuple(coords.mean(axis=0).tolist())
    return center, k, len(pts)


def compute_site_box_from_receptor(receptor_pdbqt, lig_resname=None, fixed_edge=FIXED_EDGE):
    """
    Center the grid at the co-crystallized ligand (HETATM) center.
    If lig_resname is given, use that; otherwise auto-detect the largest non-excluded HET group.
    """
    groups = parse_pdbqt_het_groups(receptor_pdbqt)
    center, key, n = pick_ligand_group(groups, target_resn=lig_resname)
    # Fixed cube ~25 Å
    size = (fixed_edge, fixed_edge, fixed_edge)
    # key is (resn, chain, resi)
    return center, size, key, n


def run_vina(receptor, ligand, center, size, out_path, log_path):
    cmd = [
        vina_executable,
        "--receptor", receptor,
        "--ligand", ligand,
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", f"{size[0]:.3f}",
        "--size_y", f"{size[1]:.3f}",
        "--size_z", f"{size[2]:.3f}",
        "--out", out_path,
        "--log", log_path
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def parse_affinities(log_path, top_n=10):
    aff = []
    if not os.path.exists(log_path):
        return aff
    with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            parts = s.split()
            if parts and parts[0].isdigit() and len(parts) >= 2:
                try:
                    aff.append(float(parts[1]))
                except ValueError:
                    pass
            if len(aff) >= top_n:
                break
    return aff[:top_n]


def save_affinities_csv(affinities, out_csv_path):
    pd.DataFrame({
        "Rank": list(range(1, len(affinities) + 1)),
        "Binding_Affinity_kcal/mol": affinities
    }).to_csv(out_csv_path, index=False)


def run_from_csv():
    def pick(df, candidates):
        cols = {c.strip().lower(): c for c in df.columns}
        for name in candidates:
            if name in cols:
                return cols[name]
        return None

    df = pd.read_csv(csv_path)
    df.columns = [c.strip().lower() for c in df.columns]

    # Identify columns
    lig_col  = pick(df, ["pubchem_cid", "cid", "ligand_cid", "ligand", "ligand_id"])
    prot_col = pick(df, ["pdb_id", "pdbid", "pdb", "receptor", "protein", "malaria"])
    # Optional ligand residue name in receptor to center on
    ligres_col = pick(df, ["lig_resname", "site_ligand", "cocrystal_ligand", "co_ligand", "resname"])

    if lig_col is None or prot_col is None:
        raise RuntimeError(
            f"CSV columns = {list(df.columns)}\n"
            f"Need ligand→(pubchem_cid/cid/ligand) and receptor→(pdb_id/pdbid/pdb)."
        )

    print(f"[INFO] Using columns -> ligand: '{lig_col}', receptor: '{prot_col}', "
          f"ligand-resname: '{ligres_col or 'AUTO-DETECT'}'")

    pairs, missing = [], []
    for _, row in df.iterrows():
        lig_key  = str(row[lig_col]).strip()
        prot_key = str(row[prot_col]).strip()
        if not lig_key or not prot_key or lig_key.lower() == "nan" or prot_key.lower() == "nan":
            continue

        lig_resname = None
        if ligres_col:
            val = str(row[ligres_col]).strip()
            if val and val.lower() != "nan":
                lig_resname = val

        lig_path = os.path.join(lig_folder,  ensure_pdbqt_name(lig_key))
        rec_path = os.path.join(prot_folder, ensure_pdbqt_name(prot_key))

        if os.path.exists(lig_path) and os.path.exists(rec_path):
            pairs.append((rec_path, lig_path, os.path.basename(rec_path),
                          os.path.basename(lig_path), lig_resname))
        else:
            why = []
            if not os.path.exists(lig_path): why.append(f"ligand_missing:{os.path.basename(lig_path)}")
            if not os.path.exists(rec_path): why.append(f"receptor_missing:{os.path.basename(rec_path)}")
            missing.append((os.path.basename(rec_path), os.path.basename(lig_path), ";".join(why)))

    print(f"[CHECK] Valid pairs: {len(pairs)} | Missing/invalid: {len(missing)}")
    if missing:
        miss_path = os.path.join(base_dir, "missing_inputs.txt")
        with open(miss_path, "w", encoding="utf-8") as f:
            for prot_file, lig_file, why in missing:
                f.write(f"{prot_file}\t{lig_file}\t{why}\n")
        print(f"[CHECK] Missing file list saved to: {miss_path}")

    if not pairs:
        print("[ABORT] No valid {receptor, ligand} pairs found.")
        return

    failed, success, skipped = [], 0, 0
    for rec_path, lig_path, prot_file, lig_file, lig_resname in pairs:
        receptor_stem = os.path.splitext(prot_file)[0]
        ligand_stem   = os.path.splitext(lig_file)[0]
        pair_stem     = f"{receptor_stem}__{ligand_stem}"

        out_pdbqt = os.path.join(output_folder, pair_stem + ".pdbqt")
        out_log   = os.path.join(log_folder,   pair_stem + ".log")
        out_csv   = os.path.join(csv_folder,   pair_stem + ".csv")

        if os.path.exists(out_csv):
            skipped += 1
            continue

        try:
            center, size, group_key, n_atoms = compute_site_box_from_receptor(
                rec_path, lig_resname=lig_resname, fixed_edge=FIXED_EDGE
            )
            resn, chain, resi = group_key
            print(f"[{pair_stem}] site={resn}:{chain}{resi} (n={n_atoms}) "
                  f"center={tuple(round(x,3) for x in center)} size={size}")

            run_vina(
                receptor=rec_path,
                ligand=lig_path,
                center=center,
                size=size,
                out_path=out_pdbqt,
                log_path=out_log
            )

            affinities = parse_affinities(out_log, top_n=10)
            if affinities:
                save_affinities_csv(affinities, out_csv)
                success += 1
                print(f"[OK] {pair_stem} -> {out_csv}")
            else:
                failed.append((lig_file, prot_file, "no affinities in log"))
                print(f"[WARN] {pair_stem}: no affinities parsed")

        except subprocess.CalledProcessError as e:
            failed.append((lig_file, prot_file, f"vina error: {e.stderr.strip()}"))
            print(f"[FAIL] {pair_stem}: Vina error")
        except Exception as e:
            failed.append((lig_file, prot_file, f"{type(e).__name__}: {e}"))
            print(f"[FAIL] {pair_stem}: {e}")

    fail_path = os.path.join(base_dir, "failed_dockings.txt")
    if failed:
        with open(fail_path, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["ligand_file", "protein_file", "reason"])
            w.writerows(failed)
        print(f"\n[SUMMARY] {success} succeeded, {skipped} skipped, {len(failed)} failed.")
        print(f"[LOG] See: {fail_path}")
    else:
        print(f"\n[SUMMARY] {success} succeeded, {skipped} skipped, 0 failed. All pairs processed.")


if __name__ == "__main__":
    run_from_csv()
