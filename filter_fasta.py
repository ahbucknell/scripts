#!/usr/bin/env python3
"""
Filter a FASTA file to retain sequences strictly under a given length threshold.
Usage: python filter_fasta.py input.fasta output.fasta 1000
"""

import sys
from Bio import SeqIO

def filter_fasta(input_file, output_file, max_length):
    all_seqs = list(SeqIO.parse(input_file, "fasta"))
    filtered = [rec for rec in all_seqs if len(rec.seq) < max_length]
    removed = len(all_seqs) - len(filtered)

    SeqIO.write(filtered, output_file, "fasta")

    print(f"Total sequences:   {len(all_seqs)}")
    print(f"Sequences retained: {len(filtered)} (< {max_length} aa)")
    print(f"Sequences removed:  {removed} (>= {max_length} aa)")
    print(f"Output written to:  {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_fasta.py input.fasta output.fasta max_length")
        print("Example: python filter_fasta.py proteins.fasta filtered.fasta 1000")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    try:
        max_length = int(sys.argv[3])
    except ValueError:
        print("Error: max_length must be an integer")
        sys.exit(1)

    filter_fasta(input_file, output_file, max_length)
