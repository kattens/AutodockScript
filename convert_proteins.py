"""
Batch-convert protein PDBs in ./malaria to PDBQT in ./malaria_pdbqt
using prepare_receptor4.py (AutoDockTools style).
"""

import os
import sys
import subprocess

base_dir = os.path.dirname(os.path.abspath(__file__))

input_folder = os.path.join(base_dir, "malaria")
output_folder = os.path.join(base_dir, "malaria_pdbqt")
os.makedirs(output_folder, exist_ok=True)

prepare_script = os.path.join(base_dir, "prepare_receptor4.py")
log_file_path = os.path.join(base_dir, "failed_conversions_proteins.txt")

failed = []
processed = 0
skipped = 0

for filename in os.listdir(input_folder):
    if not filename.lower().endswith(".pdb"):
        continue

    protein_path = os.path.join(input_folder, filename)
    output_name = os.path.splitext(filename)[0] + ".pdbqt"
    output_path = os.path.join(output_folder, output_name)

    # Skip if already converted
    if os.path.exists(output_path):
        skipped += 1
        continue

    cmd = [
        sys.executable,            # use current Python interpreter
        prepare_script,
        "-r", protein_path,
        "-o", output_path,
        "-v"
    ]

    try:
        subprocess.run(
            cmd,
            cwd=base_dir,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        processed += 1
    except subprocess.CalledProcessError as e:
        print(f"[FAIL] {filename}")
        # Optional: print e.stderr for debugging
        failed.append(f"{filename}\nSTDERR:\n{e.stderr}\n")

# Write failures (if any)
if failed:
    with open(log_file_path, "w", encoding="utf-8") as log:
        log.write("Failed protein conversions:\n\n")
        for entry in failed:
            log.write(entry + "\n")
    print(f"\n[SUMMARY] {processed} converted, {skipped} skipped, {len(failed)} failed.")
    print(f"[LOG] See: {log_file_path}")
else:
    print(f"\n[SUMMARY] {processed} converted, {skipped} skipped, 0 failed. All proteins processed successfully.")
