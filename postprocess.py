import pandas as pd
import os

#expdata = pd.read_csv('EXPFULL')
#compdata = pd.read_csv('COMPFULL')

csv_folder = 'top10_affinities'

"""
for file in os.listdir(csv_folder):
    if "_vs" in file:
        new_name = file.split("_vs")[0] + os.path.splitext(file)[1]
        os.rename(os.path.join(csv_folder, file),
                  os.path.join(csv_folder, new_name))
        print(f"{file} -> {new_name}")
"""


import os
import pandas as pd

# ===== CONFIG =====
main_csv = "Data.csv"     # your main CSV
results_folder = csv_folder           # folder containing files like 8TZE_1V0B_T3X.csv
key_col = "3D Interaction"              # match this column (values like ... .pdb)
# Source column(s) present in result CSVs; we'll pick the first that exists
source_aff_cols = ["Binding_Affinity", "Binding_Affinity_kcal/mol"]
# Column to write into main CSV (use kcal/mol since that's what's in your files)
out_col = "Binding_Affinity_kcal/mol"
out_csv = "Data_updated.csv"            # output file
# ==================
"""
def normalize(name: str) -> str:
    
    s = str(name).strip()
    for ext in (".pdbqt", ".pdb", ".csv"):
        if s.lower().endswith(ext):
            s = s[: -len(ext)]
            break
    # also strip any remaining extension, if present
    base, _ = os.path.splitext(s)
    return base

# Load main CSV
df = pd.read_csv(main_csv)

# Ensure output column exists
if out_col not in df.columns:
    df[out_col] = pd.NA

# Build mapping: normalized base name -> csv file path
file_map = {}
for f in os.listdir(results_folder):
    if f.lower().endswith(".csv"):
        base = normalize(f)
        file_map[base] = os.path.join(results_folder, f)

matched = 0
for idx, val in df[key_col].items():
    key = normalize(val)
    path = file_map.get(key)
    if not path:
        continue

    try:
        sub = pd.read_csv(path)
        # pick the first source affinity column that exists
        src_col = next((c for c in source_aff_cols if c in sub.columns), None)
        if src_col is None or sub.empty:
            continue
        affinity_val = sub[src_col].iloc[0]
        df.at[idx, out_col] = affinity_val
        matched += 1
    except Exception as e:
        print(f"[warn] failed reading {path}: {e}")

df.to_csv(out_csv, index=False)
print(f"[done] wrote {out_csv} | rows updated: {matched}")
"""



#check the column "Binding_Affinity_kcal/mol" if value >= -6 (higher than), remove row

df = pd.read_csv("Comp_Data_updated.csv")  # your file
df = df[df["Binding_Affinity_kcal/mol"] < -6]  # keep only rows with value < -6
df.to_csv("Comp_Data_updated.csv", index=False)
