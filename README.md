# scripts

Assorted one-off bioinformatics scripts.

## ko_verify.sh

Checks whether one CRISPR + cassette knock-in strain is correct, from paired-end
WGS. One strain, one target locus, one submission.

Maps R1/R2 as a pair against the genome and cassette together, then reports three
things from that single alignment: coverage dropout at the target, cassette-to-genome
junction reads at each flank, and ectopic insertions with a copy-number estimate.
`PASS` needs all four — target deleted, both flanks supported, copy number ~1x,
nothing off-target.

```bash
# Build the reference once. It is just a concatenation.
cat MGGv8_genome.fasta cassette.fasta > chimeric.fasta
bwa index chimeric.fasta && samtools faidx chimeric.fasta
grep '>' cassette.fasta          # this name is --cassette-contig

sbatch ko_verify.sh --sample KO_A --r1 A_1.fq.gz --r2 A_2.fq.gz \
       --ref chimeric.fasta --cassette-contig HPH_CASSETTE \
       --target contig1:50000-52000 --outdir results
```

Writes `results/KO_A/KO_A.ko_report.tsv` — one row, ending in
`PASS` / `REVIEW` / `FAIL`. `--help` for the rest, including how to read a
`REVIEW`.
