import os
import shutil
from pathlib import Path
import pandas as pd

# --- Config ---
csv_path = 'CompData.csv'            # file OR folder of CSVs
top10_affinities = 'top10_affinities'
docking_log = 'docking_logs'
docking_output = 'docking_output'
Out_folder = 'OutFolder'
col_name = "3D Interaction"
# --------------

base = Path('.').resolve()
src_dirs = [base / top10_affinities, base / docking_log, base / docking_output]
out_dir = base / Out_folder
out_dir.mkdir(parents=True, exist_ok=True)

def load_ids(csv_path: str, column: str) -> set:
    p = Path(csv_path)
    frames = []
    if p.is_dir():
        csvs = sorted(p.glob("*.csv"))
        if not csvs:
            raise FileNotFoundError(f"No CSV files found in folder: {p}")
        for f in csvs:
            try:
                df = pd.read_csv(f)
                if column in df.columns:
                    frames.append(df[[column]])
            except Exception:
                pass
        if not frames:
            raise KeyError(f"Column '{column}' not found in any CSV in {p}")
        df_all = pd.concat(frames, ignore_index=True)
    else:
        if not p.exists():
            raise FileNotFoundError(f"CSV file not found: {p}")
        df_all = pd.read_csv(p)
        if column not in df_all.columns:
            raise KeyError(f"Column '{column}' not in {p}")
        df_all = df_all[[column]]

    vals = (
        df_all[column]
        .dropna()
        .astype(str)
        .str.strip()
        .replace({"": None})
        .dropna()
        .unique()
        .tolist()
    )
    return set(vals)

def strip_ext(name: str) -> str:
    p = Path(name)
    stem = p.stem
    if stem.endswith(".pdb") or stem.endswith(".pdbqt") or stem.endswith(".csv") or stem.endswith(".log"):
        stem = Path(stem).stem
    return stem

def filename_matches(target: str, fname: str) -> bool:
    # match base equality OR contains either way
    t = strip_ext(target)
    f = strip_ext(fname)
    if not t or not f:
        return False
    return (t == f) or (t in f) or (f in t)

def safe_copy(src: Path, dst_dir: Path):
    dst = dst_dir / src.name
    if not dst.exists():
        shutil.copy2(src, dst)
        return dst
    stem, suffix = dst.stem, dst.suffix
    i = 1
    while True:
        candidate = dst_dir / f"{stem}_{i}{suffix}"
        if not candidate.exists():
            shutil.copy2(src, candidate)
            return candidate
        i += 1

def main():
    targets = load_ids(csv_path, col_name)
    if not targets:
        print("[info] No IDs found to match.")
        return

    # pre-normalize targets once
    norm_targets = {strip_ext(t) for t in targets if strip_ext(t)}

    copied = 0
    matched_targets = set()
    missing_targets = set(norm_targets)

    for sdir in src_dirs:
        if not sdir.exists():
            print(f"[warn] Source folder missing: {sdir}")
            continue
        for f in sdir.rglob("*"):
            if not f.is_file() or f.name.startswith("."):
                continue
            # check against ALL targets (do NOT discard after first match)
            for t in norm_targets:
                if filename_matches(t, f.name):
                    safe_copy(f, out_dir)
                    copied += 1
                    matched_targets.add(t)
                    # note: do NOT stop looking for t in other folders/files
                    # we want ALL matches across ALL three dirs

    # report unmatched targets (those that never matched any file)
    missing_targets = norm_targets - matched_targets

    print(f"[done] Copied {copied} file(s) into '{out_dir}'.")
    print(f"[info] Matched {len(matched_targets)} unique target(s) across all folders.")
    if missing_targets:
        print(f"[note] {len(missing_targets)} target(s) had no file match (showing up to 20):")
        for t in list(missing_targets)[:20]:
            print("   -", t)
        if len(missing_targets) > 20:
            print("   ...")

if __name__ == "__main__":
    main()
