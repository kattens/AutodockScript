import os, sys, re, subprocess
from pathlib import Path

base_dir = Path(__file__).resolve().parent
input_folder = base_dir / "malaria"
output_folder = base_dir / "malaria_pdbqt"
output_folder.mkdir(parents=True, exist_ok=True)
log_file_path = base_dir / "failed_conversions_proteins.txt"
failed = []

def run(cmd):
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, shell=False)
    return proc.returncode, proc.stdout

def has_adt():
    try:
        import AutoDockTools_py3  # noqa
        return True
    except Exception:
        return False

def adt_prepare_receptor(pdb_in: Path, pdbqt_out: Path):
    cmd = [
        sys.executable, "-m", "AutoDockTools_py3.Utilities24.prepare_receptor4",
        "-r", str(pdb_in),
        "-o", str(pdbqt_out),
        "-A", "hydrogens",
        "-U", "nphs_lps_waters_nonstdres",
        "-v",
    ]
    return run(cmd)

def obabel_prepare_receptor(pdb_in: Path, pdbqt_out: Path):
    # Explicit -i/-o is more reliable on Windows
    cmd = [
        "obabel",
        "-ipdb", str(pdb_in),
        "-opdbqt", "-O", str(pdbqt_out),
        "-xr", "-xh",
    ]
    return run(cmd)

def _fix_atom_line_columns(line: str) -> str:
    line = line.replace("\t", " ")
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return line
    if len(line) < 80:
        line = line.rstrip("\n") + " " * (80 - len(line)) + "\n"

    def try_cols(s):
        try:
            return float(s[30:38]), float(s[38:46]), float(s[46:54])
        except ValueError:
            return None

    xyz = try_cols(line)
    if xyz is None:
        m = None
        for m in re.finditer(r'(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)', line):
            pass
        if not m:
            return line
        try:
            xyz = tuple(float(v) for v in m.groups())
        except ValueError:
            return line

    x, y, z = xyz
    coord = f"{x:8.3f}{y:8.3f}{z:8.3f}"
    line_list = list(line)
    line_list[30:54] = list(coord)
    return "".join(line_list)

def sanitize_and_reflow(src: Path, dst: Path):
    with src.open("r") as fin, dst.open("w") as fout:
        for raw in fin:
            if raw.startswith(("ATOM", "HETATM")):
                raw = _fix_atom_line_columns(raw)
            else:
                raw = raw.replace("\t", " ")
            fout.write(raw)

def validate_pdbqt_coords(path: Path, max_lines: int = 200):
    bad = []
    with path.open("r") as f:
        seen = 0
        for i, line in enumerate(f, start=1):
            if line.startswith(("ATOM", "HETATM")):
                seen += 1
                try:
                    float(line[30:38]); float(line[38:46]); float(line[46:54])
                except ValueError:
                    bad.append((i, line.rstrip("\n")))
                if seen >= max_lines:
                    break
    return bad

print(f"Starting conversion of PDB files from: {input_folder}")

for filename in os.listdir(input_folder):
    if not filename.lower().endswith(".pdb"):
        continue

    pdb_in = (input_folder / filename).resolve()
    pdbqt_out = (output_folder / (pdb_in.stem + ".pdbqt")).resolve()
    tmp_out = pdbqt_out.with_suffix(".pdbqt.tmp")
    clean_out = pdbqt_out.with_suffix(".pdbqt.clean")

    print(f"\nProcessing: {pdb_in.name}")

    ok = False
    if has_adt():
        code, out = adt_prepare_receptor(pdb_in, tmp_out)
        ok = (code == 0)
        if not ok:
            print("  ADT failed, will try Open Babel.")
            # print(out)  # uncomment to see ADT details
    else:
        print("  ADT not installed in this env, will try Open Babel.")

    if not ok:
        code2, out2 = obabel_prepare_receptor(pdb_in, tmp_out)
        ok = (code2 == 0)
        if not ok:
            print("  Open Babel failed too. Check that openbabel is installed.")
            # print(out2)  # uncomment to see OBabel details
            failed.append(pdb_in.name)
            continue

    # sanitize/reflow -> fixed columns
    sanitize_and_reflow(tmp_out, clean_out)
    bad = validate_pdbqt_coords(clean_out, max_lines=200)
    if bad:
        print("  Validation found malformed lines (showing up to 5):")
        for i, L in bad[:5]:
            print(f"    line {i}: {L}")
        failed.append(pdb_in.name)
        continue

    try:
        if pdbqt_out.exists():
            pdbqt_out.unlink()
        clean_out.rename(pdbqt_out)
        tmp_out.unlink(missing_ok=True)
        print(f"  OK -> {pdbqt_out}")
    except Exception as e:
        print(f"  Error finalizing {pdbqt_out.name}: {e}")
        failed.append(pdb_in.name)

# Report
if failed:
    with open(log_file_path, "w") as log:
        log.write("Failed protein conversions:\n")
        for name in failed:
            log.write(name + "\n")
    print(f"\nLogged {len(failed)} failures to {log_file_path}")
else:
    print("\nAll receptor PDB files processed successfully.")
