import os
import pandas as pd
import requests

# === CONFIGURATION ===
csv_path = "Data.csv"       # Path to your CSV file
output_dir = "malaria"       # Folder to save downloaded PDBs
pdb_column = "pdb_id"          # Column name in CSV

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Load CSV
df = pd.read_csv(csv_path)

# Iterate over PDB IDs
for idx, pdb_id in enumerate(df[pdb_column].dropna().unique(), start=1):
    pdb_id = str(pdb_id).strip().upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    save_path = os.path.join(output_dir, f"{pdb_id}.pdb")

    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            with open(save_path, "wb") as f:
                f.write(response.content)
            print(f"[{idx}] Downloaded {pdb_id}")
        else:
            print(f"[{idx}] Failed to download {pdb_id} (HTTP {response.status_code})")
    except Exception as e:
        print(f"[{idx}] Error downloading {pdb_id}: {e}")
