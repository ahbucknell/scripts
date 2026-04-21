#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter PDB files by average pLDDT score."
    )
    parser.add_argument(
        "directory", type=Path, help="Path to the directory containing PDB files."
    )
    parser.add_argument(
        "threshold",
        type=float,
        help="pLDDT threshold (0-100). Files with average pLDDT >= threshold pass.",
    )
    parser.add_argument(
        "-d",
        "--delete",
        action="store_true",
        help="Delete files that fail the threshold instead of moving them.",
    )
    return parser.parse_args()


def get_average_plddt(pdb_path: Path):
    """
    Parse a PDB file and return the average pLDDT score.

    pLDDT scores are stored in the B-factor column (column 61-66) of ATOM/HETATM records.
    AlphaFold PDB files store per-residue pLDDT in the B-factor field.
    We average over all ATOM records (one per atom), which reflects the
    per-atom pLDDT values as written by AlphaFold.
    """
    scores = []
    try:
        for line in pdb_path.read_text().splitlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    scores.append(float(line[60:66].strip()))
                except (ValueError, IndexError):
                    continue
    except OSError as e:
        print(f"  [WARNING] Could not read '{pdb_path}': {e}", file=sys.stderr)
        return None

    if not scores:
        print(
            f"  [WARNING] No ATOM/HETATM records found in '{pdb_path}'.",
            file=sys.stderr,
        )
        return None

    return sum(scores) / len(scores)


def main():
    args = parse_args()

    directory: Path = args.directory.resolve()
    threshold: float = args.threshold
    delete_mode: bool = args.delete

    if not directory.is_dir():
        print(f"Error: '{directory}' is not a valid directory.", file=sys.stderr)
        sys.exit(1)

    if not (0 <= threshold <= 100):
        print(
            "Warning: threshold is outside the typical pLDDT range of 0-100.",
            file=sys.stderr,
        )

    # Collect PDB files (non-recursive)
    pdb_files = sorted(
        p for p in directory.iterdir() if p.suffix.lower() == ".pdb" and p.is_file()
    )

    if not pdb_files:
        print(f"No PDB files found in '{directory}'.")
        sys.exit(0)

    print(f"Found {len(pdb_files)} PDB file(s) in '{directory}'.")
    print(
        f"Threshold: {threshold}  |  Mode: {'delete failures' if delete_mode else 'sort into subdirs'}\n"
    )

    # Threshold label for directory names (strip trailing zeros for cleanliness)
    t_label = f"{threshold:g}"

    if not delete_mode:
        over_dir = directory / f"over{t_label}"
        under_dir = directory / f"under{t_label}"
        over_dir.mkdir(exist_ok=True)
        under_dir.mkdir(exist_ok=True)

    passed = []
    failed = []
    skipped = []

    for fpath in pdb_files:
        avg = get_average_plddt(fpath)
        print(fpath)
        print(avg)

        if avg is None:
            skipped.append(fpath)
            continue

        if avg >= threshold:
            passed.append((fpath, avg))
        else:
            failed.append((fpath, avg))

    # --- Process passing files ---
    if not delete_mode:
        for fpath, avg in passed:
            fpath.rename(over_dir / fpath.name)
            print(f"  [PASS  {avg:6.2f}]  {fpath.name}  ->  over{t_label}/")
    else:
        for fpath, avg in passed:
            print(f"  [PASS  {avg:6.2f}]  {fpath.name}  (kept)")

    # --- Process failing files ---
    for fpath, avg in failed:
        if delete_mode:
            fpath.unlink()
            print(f"  [FAIL  {avg:6.2f}]  {fpath.name}  (deleted)")
        else:
            fpath.rename(under_dir / fpath.name)
            print(f"  [FAIL  {avg:6.2f}]  {fpath.name}  ->  under{t_label}/")

    # --- Skipped files ---
    for fpath in skipped:
        print(f"  [SKIP         ]  {fpath.name}  (could not parse pLDDT)")

    # --- Summary ---
    print(f"\nDone.")
    print(f"  Passed  : {len(passed)}")
    print(f"  Failed  : {len(failed)}")
    if skipped:
        print(f"  Skipped : {len(skipped)}")


if __name__ == "__main__":
    main()
