import argparse
import pathlib
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Move ranked_0.pdb files up one level, renaming them after their parent subdirectory."
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
        help="Strip the leading 'N_' prefix from subdirectory names (e.g. '42_my_protein' → 'my_protein.pdb').",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned moves without making any changes.",
    )
    return parser.parse_args()


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

    moved = 0
    skipped = 0
    overwritten = 0

    for pdb_file in pdb_files:
        subdir_name = pdb_file.parent.name
        if args.rename:
            stem = re.sub(r"^\d+_", "", subdir_name)
            if stem == subdir_name:
                print(
                    f"[WARNING] Could not strip prefix from '{subdir_name}', using full name."
                )
        else:
            stem = subdir_name
        destination = base_dir / f"{stem}.pdb"

        if destination.exists():
            if args.force:
                print(f"[WARNING] Destination exists, overwriting: {destination}")
                if not args.dry_run:
                    pdb_file.rename(destination)
                else:
                    print(f"  (dry run) {pdb_file} → {destination}")
                overwritten += 1
            else:
                print(f"[WARNING] Destination exists, skipping: {destination}")
                skipped += 1
            continue

        if args.dry_run:
            print(f"  (dry run) {pdb_file} → {destination}")
        else:
            pdb_file.rename(destination)

        moved += 1

    print(f"\n--- Summary ---")
    if args.dry_run:
        print(f"  Would move:      {moved}")
        print(f"  Would overwrite: {overwritten}")
        print(f"  Would skip:      {skipped}")
    else:
        print(f"  Moved:      {moved}")
        print(f"  Overwritten: {overwritten}")
        print(f"  Skipped:    {skipped}")


if __name__ == "__main__":
    main()
