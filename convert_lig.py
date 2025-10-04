"""
Batch-convert ligand SDFs in ./ligand to PDBQT in ./ligands_pdbqt.
Steps:
  1) SDF -> PDB (Open Babel: obabel)
  2) PDB -> PDBQT (prep_ligands.py / prepare_ligand4.py)
"""

import os
import sys
import subprocess
from pathlib import Path

base_dir = Path(__file__).resolve().parent
input_folder = base_dir / "ligand"
output_folder = base_dir / "ligands_pdbqt"
tmp_folder = base_dir / "_tmp_ligand_pdb"
output_folder.mkdir(parents=True, exist_ok=True)
tmp_folder.mkdir(parents=True, exist_ok=True)

prepare_script = base_dir / "prep_ligands.py"  # your prepare_ligand4 port
log_file_path = base_dir / "failed_conversions.txt"

failed = []
processed = 0
skipped = 0

def run(cmd, cwd=None):
    return subprocess.run(
        cmd,
        cwd=cwd,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        shell=False
    )

for filename in os.listdir(input_folder):
    if not filename.lower().endswith(".sdf"):
        continue

    sdf_path = input_folder / filename
    stem = sdf_path.stem
    pdb_tmp = tmp_folder / f"{stem}.pdb"
    out_pdbqt = output_folder / f"{stem}.pdbqt"

    # Skip if already converted
    if out_pdbqt.exists():
        skipped += 1
        continue

    try:
        # 1) SDF -> PDB using Open Babel (requires obabel in PATH)
        # If you prefer to include hydrogens at this step: add "--addhs"
        run(["obabel", str(sdf_path), "-O", str(pdb_tmp)])

        # 2) PDB -> PDBQT using your ligand prep script
        cmd = [
            sys.executable,
            str(prepare_script),
            "-l", str(pdb_tmp),
            "-o", str(out_pdbqt),
            "-v"
        ]
        run(cmd, cwd=base_dir)

        processed += 1

    except subprocess.CalledProcessError as e:
        failed.append(
            f"{filename}\nCMD: {' '.join(e.cmd)}\nSTDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}\n"
        )
    except FileNotFoundError as e:
        # Likely 'obabel' not found
        failed.append(f"{filename}\nERROR: {e}\n")

# Write failures (if any)
if failed:
    with open(log_file_path, "w", encoding="utf-8") as log:
        log.write("Failed ligand conversions:\n\n")
        for entry in failed:
            log.write(entry + "\n")
    print(f"\n[SUMMARY] {processed} converted, {skipped} skipped, {len(failed)} failed.")
    print(f"[LOG] See: {log_file_path}")
else:
    print(f"\n[SUMMARY] {processed} converted, {skipped} skipped, 0 failed. All ligands processed successfully.")
