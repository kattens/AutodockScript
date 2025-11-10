#seperate the csv by the drug name, cleanup, add them back together
import pandas as pd
import os

# ===== CONFIG =====
input_csv = "Comp_Data_updated.csv"     # your main CSV file
output_folder = "out"      # folder to save separated CSVs
drug_col = "drug ligand name"     # column that contains drug names
# ==================
"""
# Make sure output folder exists
os.makedirs(output_folder, exist_ok=True)

# Read the main CSV
df = pd.read_csv(input_csv)

# Check column
if drug_col not in df.columns:
    raise ValueError(f"Column '{drug_col}' not found in CSV.")

# Group by drug name and save each group
for drug, group in df.groupby(drug_col):
    safe_name = "".join(c if c.isalnum() or c in ('-', '_') else '_' for c in str(drug))
    out_path = os.path.join(output_folder, f"{safe_name}.csv")
    group.to_csv(out_path, index=False)
    print(f"[OK] Saved {out_path} ({len(group)} rows)")
"""

#merge all the cleaned csvs back together
import pandas as pd
import os
from pathlib import Path

# ===== CONFIG =====
input_folder = Path("out")        # folder with cleaned CSV files
output_csv = "merged.csv"         # final merged CSV
# ==================

# Collect all CSV files in the folder
csv_files = list(input_folder.glob("*.csv"))

# Merge all into one DataFrame
dfs = []
for file in csv_files:
    try:
        df = pd.read_csv(file)
        df["Source_File"] = file.name  # optional: track which file it came from
        dfs.append(df)
        print(f"[merged] {file.name} ({len(df)} rows)")
    except Exception as e:
        print(f"[skip] {file.name} ({e})")

# Concatenate all DataFrames and save
if dfs:
    merged_df = pd.concat(dfs, ignore_index=True)
    merged_df.to_csv(output_csv, index=False)
    print(f"\nMerged {len(csv_files)} files into '{output_csv}' ({len(merged_df)} total rows)")
else:
    print("No valid CSV files found.")
