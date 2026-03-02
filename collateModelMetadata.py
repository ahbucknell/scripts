import subprocess
import sys
from genericpath import exists
from pathlib import Path

import Bio
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


targetDir = Path("/Users/doz23per/Documents/controls-HMA")
tmpDir = targetDir / "tmp"
ipsaeFile = downloadIPSAE(targetDir)
pdbFiles = list(targetDir.glob("*.pdb"))

for p in pdbFiles:
    cifPath = tmpDir / p.with_suffix(".cif").name
    if cifPath.exists():
        print(cifPath.name, " already exists. Skipping...")
        continue
    _ = convertToCIF(p, cifPath)

cifFiles = list(tmpDir.glob("*.cif"))
for p in cifFiles:
    paeFile = targetDir / (p.stem + "_pae.npz")
    if not paeFile.exists():
        print(f"Cannot find file: {paeFile.name}")
        continue
    if (tmpDir / (p.stem + "_10_10.txt")).exists():
        print("IPSAE file already exists. Skipping...")
        continue
    runIPSAE(p, paeFile, ipsaeFile)
