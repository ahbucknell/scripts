#!/usr/bin/env python3
import sys
from pathlib import Path

if len(sys.argv) != 3:
    sys.exit("usage: convert_topology_to_gff3.py input.txt output.gff3")

inp = Path(sys.argv[1])
out = Path(sys.argv[2])

seq_lengths = {}
features = {}
current_seqid = None
feature_counters = {}

with inp.open() as fh:
    for raw in fh:
        line = raw.strip()
        if not line or line == "//":
            continue
        if line.startswith("##gff-version"):
            continue
        if line.startswith("#"):
            parts = line[1:].strip().split()
            if len(parts) >= 3 and parts[1] == "Length:":
                seq_lengths[parts[0]] = parts[2]
            continue

        cols = line.split("\t")
        if len(cols) < 4:
            continue

        seqid, ftype, start, end = cols[0], cols[1], cols[2], cols[3]

        feature_counters.setdefault(seqid, 0)
        feature_counters[seqid] += 1

        if seqid not in features:
            features[seqid] = []

        attrs = f"ID={seqid}.{feature_counters[seqid]};Name={ftype}"
        features[seqid].append([
            seqid,
            ".",
            ftype,
            start,
            end,
            ".",
            ".",
            ".",
            attrs
        ])

with out.open("w") as fh:
    fh.write("##gff-version 3\n")
    for seqid, rows in features.items():
        if seqid in seq_lengths:
            fh.write(f"##sequence-region {seqid} 1 {seq_lengths[seqid]}\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
