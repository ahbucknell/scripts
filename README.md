# scripts

Assorted one-off bioinformatics scripts.

## KO strain verification

Checks whether CRISPR + cassette knock-in strains are correct, from paired-end WGS.

| Script | Purpose |
|---|---|
| `ko_prepare_ref.sh` | Build + index the chimeric reference (genome + cassette). Run once. |
| `ko_verify.sh` | Verify one strain. |
| `ko_verify_array.sbatch` | Run `ko_verify.sh` across a samplesheet as a SLURM array. |

Maps each strain once to genome + cassette combined, then reports three things
from that alignment: coverage dropout at the target, cassette-to-genome junction
reads at each flank, and ectopic insertions with a copy-number estimate. A strain
passes only if the target is deleted, both flanks have junction support, copy
number is ~1x, and nothing landed off-target.

```bash
mkdir -p logs
./ko_prepare_ref.sh --genome MGGv8_genome.fasta --cassette cassette.fasta \
                    --outdir chimeric_ref

# samples.tsv: header + sample<TAB>r1<TAB>r2<TAB>contig:start-end[,contig:start-end]
sbatch --array=1-$(($(wc -l < samples.tsv) - 1)) ko_verify_array.sbatch \
       samples.tsv chimeric_ref/chimeric.fasta results
```

Results land in `results/<sample>/<sample>.ko_report.tsv`, one row per target,
ending in a `PASS` / `REVIEW` / `FAIL` verdict. `--help` on any script for options.

Single strain, without SLURM:

```bash
./ko_verify.sh --sample KO_A --r1 A_1.fq.gz --r2 A_2.fq.gz \
               --ref chimeric_ref/chimeric.fasta \
               --target contig1:50000-52000 --outdir results
```
