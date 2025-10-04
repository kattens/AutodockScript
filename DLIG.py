import os
import pandas as pd
import requests
from tqdm import tqdm

# === config ===
csv_path = "Data.csv"  
save_dir = "ligand"      # output folder
os.makedirs(save_dir, exist_ok=True)

# load csv
df = pd.read_csv(csv_path)

# iterate rows and download 3D SDF
for cid in tqdm(df["pubchem_cid"].dropna().astype(str), desc="Downloading SDFs"):
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/conformers/"
        f"0091B5F30000000E/SDF?response_type=save&response_basename=Conformer3D_COMPOUND_CID_{cid}"
    )

    out_file = os.path.join(save_dir, f"{cid}.sdf")

    if os.path.exists(out_file):
        continue  # skip if already downloaded

    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200 and len(r.content) > 0:
            with open(out_file, "wb") as f:
                f.write(r.content)
        else:
            print(f"[WARN] CID {cid}: no 3D SDF available (HTTP {r.status_code}).")
    except Exception as e:
        print(f"[ERROR] CID {cid}: {e}")