import argparse
import json
import pathlib
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Move ranked_0.pdb and companion .pkl files up one level, "
            "renaming them after their parent subdirectory."
        )
    )
    parser.add_argument(
        "base_dir",
        type=pathlib.Path,
        help="Path to the base directory containing subdirectories.",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Overwrite existing destination files (default is to skip).",
    )
    parser.add_argument(
        "-r",
        "--rename",
        action="store_true",
        help="Strip the leading 'N_' prefix from subdirectory names (e.g. '42_my_protein' → 'my_protein').",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned moves without making any changes.",
    )
    return parser.parse_args()


def resolve_stem(subdir_name, rename):
    """Return the destination stem for a given subdirectory name."""
    if rename:
        stem = re.sub(r"^\d+_", "", subdir_name)
        if stem == subdir_name:
            print(
                f"  [WARNING] Could not strip prefix from '{subdir_name}', using full name."
            )
        return stem
    return subdir_name


def get_pkl_source(subdir, subdir_name):
    """
    Parse ranking_debug.json in subdir and return the Path of the .pkl file
    named by the first element of the 'order' array.
    Returns None and prints a warning if anything goes wrong.
    """
    json_file = subdir / "ranking_debug.json"
    if not json_file.exists():
        print(
            f"  [WARNING] Missing ranking_debug.json in '{subdir_name}', skipping both files."
        )
        return None

    try:
        with json_file.open() as f:
            data = json.load(f)
        first_entry = data["order"][0]
    except (json.JSONDecodeError, KeyError, IndexError) as e:
        print(
            f"  [WARNING] Could not read 'order' from ranking_debug.json in '{subdir_name}' ({e}), skipping both files."
        )
        return None

    pkl_file = subdir / f"result_{first_entry}.pkl"
    if not pkl_file.exists():
        print(
            f"  [WARNING] PKL file '{pkl_file.name}' not found in '{subdir_name}', skipping both files."
        )
        return None

    return pkl_file


def move_file(src, dst, force, dry_run):
    """
    Attempt to move src to dst, respecting force and dry_run flags.
    Returns one of: 'moved', 'overwritten', 'skipped'.
    """
    if dst.exists():
        if force:
            print(f"  [WARNING] Destination exists, overwriting: {dst.name}")
            if not dry_run:
                src.rename(dst)
            return "overwritten"
        else:
            print(f"  [WARNING] Destination exists, skipping: {dst.name}")
            return "skipped"

    if dry_run:
        print(f"  (dry run) {src.name} → {dst}")
    else:
        src.rename(dst)
    return "moved"


def main():
    args = parse_args()
    base_dir = args.base_dir.resolve()

    if not base_dir.is_dir():
        print(f"Error: '{base_dir}' is not a valid directory.", file=sys.stderr)
        sys.exit(1)

    pdb_files = sorted(base_dir.glob("*/ranked_0.pdb"))

    if not pdb_files:
        print("No 'ranked_0.pdb' files found in any subdirectory.")
        sys.exit(0)

    if args.dry_run:
        print("=== DRY RUN — no files will be moved ===\n")

    counts = {"moved": 0, "overwritten": 0, "skipped": 0}

    for pdb_file in pdb_files:
        subdir = pdb_file.parent
        subdir_name = subdir.name
        stem = resolve_stem(subdir_name, args.rename)

        # Resolve companion .pkl before moving anything — skip both if missing
        pkl_file = get_pkl_source(subdir, subdir_name)
        if pkl_file is None:
            continue

        pdb_dest = base_dir / f"{stem}.pdb"
        pkl_dest = base_dir / f"{stem}.pkl"

        print(f"[{subdir_name}]")
        counts[move_file(pdb_file, pdb_dest, args.force, args.dry_run)] += 1
        counts[move_file(pkl_file, pkl_dest, args.force, args.dry_run)] += 1

    label = "Would move" if args.dry_run else "Moved"
    print(f"\n--- Summary (counts include both .pdb and .pkl files) ---")
    print(f"  {label}:       {counts['moved']}")
    print(f"  Overwritten: {counts['overwritten']}")
    print(f"  Skipped:     {counts['skipped']}")


if __name__ == "__main__":
    main()
