import pandas as pd
from pathlib import Path

# --- Config ---
result_folder = Path('results')
csv_file = 'CompData.csv'
aff_col = "Binding_Affinity_kcal/mol"   # column inside the small CSVs
id_col = "3D Interaction"               # column in CompData to match folder names
n_examples = 3                          # number of example affinity columns
# -------------------------------------------------------
"""
def normalize(name: str) -> str:
    name = str(name).strip()
    for ext in [".pdb", ".pdbqt", ".csv"]:
        if name.endswith(ext):
            name = name[:-len(ext)]
    return name

# Load main table
df = pd.read_csv(csv_file)

# Make sure ex1, ex2, ex3 exist
for i in range(1, n_examples + 1):
    col_name = f"ex{i}"
    if col_name not in df.columns:
        df[col_name] = pd.NA

# Process each row
for idx, row in df.iterrows():

    raw_id = str(row[id_col])
    folder_name = normalize(raw_id)

    # Path: results/<interaction>/csvs/
    folder = result_folder / folder_name / "csvs"

    if not folder.is_dir():
        print(f"[skip] Missing folder: {folder}")
        continue

    # Collect CSV files inside csvs/
    csv_paths = sorted([p for p in folder.iterdir() if p.suffix.lower() == ".csv"])

    if not csv_paths:
        print(f"[skip] No CSV files in {folder}")
        continue

    # Read first N csv files
    for j, csv_path in enumerate(csv_paths[:n_examples]):
        try:
            df_small = pd.read_csv(csv_path)
            if aff_col not in df_small.columns or df_small.empty:
                print(f"[warn] Missing column or empty file: {csv_path}")
                continue

            # First row affinity
            affinity = df_small.loc[0, aff_col]

            # Write into ex1/ex2/ex3
            df.at[idx, f"ex{j+1}"] = affinity

        except Exception as e:
            print(f"[error] Failed reading {csv_path}: {e}")

# Save output
out_file = "CompData_with_examples.csv"
df.to_csv(out_file, index=False)

print(f"[DONE] Saved updated CSV → {out_file}")


import pandas as pd
import matplotlib.pyplot as plt

# Load your updated CSV
df = pd.read_csv("CompData_with_examples.csv")

# Keep only rows where ex2 has a value
df = df.dropna(subset=["ex2"])

# (Optional but recommended) ensure numeric
df["autodock_binding_affinity"] = pd.to_numeric(df["autodock_binding_affinity"], errors="coerce")
df["ex2"] = pd.to_numeric(df["ex2"], errors="coerce")

# Drop rows where conversion failed
df = df.dropna(subset=["autodock_binding_affinity", "ex2"])

# Prepare data for boxplot
data = [
    df["autodock_binding_affinity"],  # positive (native)
    df["ex2"],                        # negative (non-native)
]
labels = ["Native (AutoDock)", "Non-native (ex2)"]

# Make boxplot
plt.figure(figsize=(8, 6))
plt.boxplot(
    data,
    labels=labels,
    showmeans=True,
    meanline=True,
)

plt.ylabel("Binding Affinity (kcal/mol)")
plt.title("Native vs Non-Native Docking Affinities (Filtered by ex2 presence)")
plt.grid(axis="y", linestyle="--", alpha=0.4)

plt.tight_layout()
plt.show()
"""

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("CompData_with_examples.csv")

# Keep only rows where ex2 exists and both cols are numeric
df = df.dropna(subset=["ex2", "autodock_binding_affinity"]).copy()
df["autodock_binding_affinity"] = pd.to_numeric(df["autodock_binding_affinity"], errors="coerce")
df["ex2"] = pd.to_numeric(df["ex2"], errors="coerce")
df = df.dropna(subset=["autodock_binding_affinity", "ex2"])

native = df["autodock_binding_affinity"]
nonnative = df["ex2"]

print("=== Native (AutoDock) stats ===")
print(native.describe())

print("\n=== Non-native (ex2) stats ===")
print(nonnative.describe())

print("\nMean difference (Native - Non-native):")
print(native.mean() - nonnative.mean())


# True if non-native is better (more negative)
mask_nonnative_better = nonnative < native

count_nonnative_better = mask_nonnative_better.sum()
total = len(df)

print(f"Non-native better in {count_nonnative_better}/{total} "
      f"({100*count_nonnative_better/total:.1f}%) of cases")

# Optional: inspect these "problematic" ones
problem_cases = df[mask_nonnative_better]
print("\nExamples where non-native is better:")
print(problem_cases.head())

from scipy.stats import wilcoxon

# Diff > 0  ⇒ native is stronger on that pair (because ex2 - native > 0)
diff = nonnative - native

stat, p = wilcoxon(diff)
print("Wilcoxon signed-rank test:")
print("  statistic =", stat)
print("  p-value  =", p)


plt.figure(figsize=(8, 5))
plt.hist(native, bins=20, alpha=0.5, label="Native (AutoDock)")
plt.hist(nonnative, bins=20, alpha=0.5, label="Non-native (ex2)")

plt.xlabel("Binding Affinity (kcal/mol)")
plt.ylabel("Count")
plt.title("Histogram of Native vs Non-native Docking Affinities")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
