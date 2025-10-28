import os
import subprocess

# === CONFIG ===
base_dir = os.path.dirname(os.path.abspath(__file__))
input_folder = os.path.join(base_dir, "ligsdf")    # folder with .sdf files
output_folder = os.path.join(base_dir, "ligpdb")   # where to save .pdb files

os.makedirs(output_folder, exist_ok=True)

# === LOOP THROUGH ALL SDF FILES ===
for filename in os.listdir(input_folder):
    if not filename.endswith(".sdf"):
        continue

    sdf_path = os.path.join(input_folder, filename)
    pdb_path = os.path.join(output_folder, os.path.splitext(filename)[0] + ".pdb")

    print(f"Converting {filename} -> {os.path.basename(pdb_path)}")
    cmd = [
        "obabel",
        "-isdf", sdf_path,    # input format and file
        "-opdb",              # output format
        "-O", pdb_path,       # output file
        "--gen3d"             # generate 3D coordinates if missing
    ]
    subprocess.run(cmd, check=True)

print("\n All SDF files converted successfully.")
