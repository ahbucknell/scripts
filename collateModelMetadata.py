import csv
import json
import subprocess
import sys
from multiprocessing import Array
from pathlib import Path

import requests
from Bio.PDB import MMCIFIO, PDBParser


# Given a path, check if ipsae.py is present.
# If not, download it from GitHub and return the path.
def downloadIPSAE(dir: Path) -> Path:
    targetPath = dir / "ipsae.py"
    if targetPath.exists():
        print(f"ipSAE script found at {targetPath}")
        return targetPath
    print("ipSAE script not found. Downloading...")
    targetURL = (
        "https://raw.githubusercontent.com/DunbrackLab/IPSAE/refs/heads/main/ipsae.py"
    )
    resp = requests.get(targetURL, timeout=30)
    resp.raise_for_status()
    _ = targetPath.write_bytes(resp.content)
    return targetPath


# Given the path of a PDB file, convert it to mmCIF.
# Return the path of the mmCIF file.
def convertToCIF(inPath: Path, outPath: Path) -> Path:
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("structure", str(inPath))
    io = MMCIFIO()
    io.set_structure(struct)
    io.save(str(outPath))
    return outPath


def runIPSAE(cifPath: Path, paePath: Path, ipsaePath: Path) -> None:
    _ = subprocess.run(
        [sys.executable, str(ipsaePath), str(paePath), str(cifPath), "10", "10"],
        check=True,
    )
    # Confirm creation of all ipSAE-related files.
    targetPattern = cifPath.parent / (cifPath.stem + "_10_10")
    fileSfxs = [".txt", ".pml", "_byres.txt"]
    ipsaeFiles = [(str(targetPattern) + s) for s in fileSfxs]
    if not all(Path(p) for p in ipsaeFiles):
        print(f"Not all ipSAE-related files created for {cifPath.stem}")


def readJsonData(path: Path, keys: list[str]) -> dict[str, str | None]:
    metrics: dict[str, str | None] = {}
    with open(path) as f:
        data = json.load(f)
    for k in keys:
        metrics[k] = data.get(k, None)
    return metrics


def readIpsaeData(path: Path, keys: list[str]) -> dict[str, str | None] | None:
    metrics: dict[str, str | None] = {}
    with open(path) as f:
        next(f)
        reader = csv.DictReader(f, delimiter=" ", skipinitialspace=True)
        for row in reader:
            if row["Type"] == "max":
                for k in keys:
                    metrics[k] = row.get(k, None)
                return metrics
    return None


def mergeData(modelID: str) -> dict[str, str | None]:
    jsonFile = Path(targetDir / (modelID + "_confidence.json"))
    tsvFile = Path(tmpDir / (modelID + "_10_10.txt"))
    jsonData = readJsonData(jsonFile, JSON_KEYS)
    tsvData = readIpsaeData(tsvFile, TSV_KEYS) or {}
    row = {"sample": modelID} | jsonData | tsvData
    print(row)
    return row


targetDir = Path("/Users/doz23per/Documents/Meps-HMA")
outPath = Path("/Users/doz23per/Documents/Meps-HMA/data.csv")
tmpDir = targetDir / "tmp"
ipsaeFile = downloadIPSAE(targetDir)
pdbFiles = list(targetDir.glob("*.pdb"))

JSON_KEYS = [
    "confidence_score",
    "ptm",
    "iptm",
    "complex_plddt",
    "complex_iplddt",
    "complex_pde",
    "complex_ipde",
]
TSV_KEYS = ["ipSAE", "pDockQ", "pDockQ2", "LIS"]

for p in pdbFiles:
    cifPath = tmpDir / p.with_suffix(".cif").name
    if cifPath.exists():
        print(cifPath.name, " already exists. Skipping...")
        continue
    _ = convertToCIF(p, cifPath)

dataRows: list[dict[str, str | None] | None] = []
cifFiles = list(tmpDir.glob("*.cif"))
for p in cifFiles:
    paeFile = targetDir / (p.stem + "_pae.npz")
    if not paeFile.exists():
        print(f"Cannot find file: {paeFile.name}")
    elif (tmpDir / (p.stem + "_10_10.txt")).exists():
        print("IPSAE file already exists. Skipping...")
    else:
        runIPSAE(p, paeFile, ipsaeFile)
    row = mergeData(p.stem)
    dataRows.append(row)

fields = ["sample"] + JSON_KEYS + TSV_KEYS
with open(outPath, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fields)
    writer.writeheader()
    writer.writerows(dataRows)
