import argparse
import json
from pathlib import Path


# Given the path of a sub-directory, return the confidences of all models.
# Return is a dictionary: {modelID: confidence_score}.
def extractConfidences(subDir: Path) -> dict[str, float]:
    confidences: dict[str, float] = {}
    pdbFiles = list(subDir.glob("*.pdb"))
    for f in pdbFiles:
        id = f.stem
        confPath = subDir / f"{id}_confidence.json"
        try:
            with confPath.open("r") as f:
                data = json.load(f)
            confidences[id] = data["confidence_score"]
        except json.JSONDecodeError:
            print(f"Invalid JSON: {confPath}")
        except FileNotFoundError:
            print(f"Missing file: {confPath}")
        except KeyError:
            print(f"Missing 'confidence_score' key: {confPath.name}")
        except Exception as e:
            print(f"Error processing {confPath.name}: {e}")
    return confidences


def moveBestModel(dirPath: Path, seed: str, modelID: str, dry_run: bool) -> None:
    files = ["_confidence.json", "_pae.npz", "_pde.npz", "_plddt.npz", ".pdb"]
    for f in files:
        oldFile = dirPath / seed / (modelID + f)
        newFile = dirPath / (modelID + "_" + seed + f)
        if not oldFile.is_file():
            print(f"Cannot find file {oldFile.name}")
            continue
        try:
            if dry_run:
                print(f"{oldFile} to {newFile}")
                continue
            if newFile.exists():
                print(f"Overwriting existing file: {newFile.name}")
                newFile.unlink()
            _ = oldFile.rename(newFile)
            print(f"Moved {oldFile.name} to {newFile.name}")
        except PermissionError:
            print(f"Permission denied: {oldFile.name}")
        except FileNotFoundError:
            print(f"Disappeared during move: {oldFile.name}")
        except OSError as e:
            print(f"OS error moving {oldFile.name}: {e}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Find best Boltz-2 models by confidence score and move to parent directory."
    )
    parser.add_argument(
        "target_dir", type=Path, help="Directory contain seed sub-directories."
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print moves without actually moving them.",
    )

    args = parser.parse_args()

    targetDir = Path(args.target_dir).absolute()
    dryRunBool: bool = args.dry_run

    if dryRunBool:
        print("Dry-run mode enabled.")

    # bestSeed stores reference models and highest confidence_score.
    # Stored as a dictionary: {modelID: (seed, confidence_score)}.
    bestSeed: dict[str, tuple[Path, float]] = {}
    subDirs = [p for p in targetDir.iterdir() if p.is_dir()]

    for s in subDirs:
        tmp_conf = extractConfidences(s)
        for k in tmp_conf:
            if k in bestSeed:
                if bestSeed[k][1] < tmp_conf[k]:
                    bestSeed[k] = (s, tmp_conf[k])
            else:
                bestSeed[k] = (s, tmp_conf[k])

    for m in bestSeed:
        tmpModel = m
        tmpSeed = bestSeed[m][0].name
        moveBestModel(targetDir, tmpSeed, tmpModel, dryRunBool)


if __name__ == "__main__":
    main()
