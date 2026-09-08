# KO strain verification from WGS

Confirm that CRISPR + cassette knock-in strains are what they are supposed to be,
using paired-end whole-genome sequencing. Reusable across strains and targets.

## Scripts

| Script | Runs | Does |
|---|---|---|
| `ko_prepare_ref.sh` | once per genome+cassette | builds and indexes the chimeric reference |
| `ko_verify.sh` | once per strain | aligns, measures, writes the report |
| `ko_verify_array.sbatch` | once per batch | SLURM array over a samplesheet |

## Run

```bash
mkdir -p logs
./ko_prepare_ref.sh --genome MGGv8_genome.fasta --cassette cassette.fasta \
                    --outdir chimeric_ref

sbatch --array=1-$(($(wc -l < samples.tsv) - 1)) ko_verify_array.sbatch \
       samples.tsv chimeric_ref/chimeric.fasta results
```

`samples.tsv` is TAB-separated with a header: `sample`, `r1`, `r2`, `targets`
(one or more `contig:start-end`, comma-separated).

## Architecture

One chimeric reference — host genome with the cassette appended as a contig
named `KO_CASSETTE` — and one paired alignment per strain. All three checks read
from that single BAM:

1. **Dropout.** Depth across the target locus against genome-wide median depth.
   A correct replacement leaves the target at ~0x.
2. **Junction.** Read pairs with one mate on the cassette and the other on unique
   genome sequence, plus `SA:Z` split alignments crossing the boundary. Binned at
   500 bp. Bins inside the target flank window are the expected integration
   junctions.
3. **Ectopic.** The same clusters landing anywhere else, plus cassette median
   depth ÷ genome median depth as a copy-number estimate.

A strain passes only if the target is deleted, **both** flanks carry junction
support, cassette copy number is ~1x, and no ectopic cluster is found.

## Decisions worth keeping

- **Cassette goes into the reference, not a separate alignment.** Reads decide
  between genome and cassette in one competitive mapping. Aligning to each
  separately makes every arm-derived read map confidently to both.
- **Cassette contigs are renamed to `KO_CASSETTE*` at prep time.** Downstream
  code finds the cassette without being told what its FASTA header said. Prep
  aborts if the genome already uses that name.
- **Reference prep is a separate script, not a step inside the verify script.**
  Parallel array tasks would otherwise race to build the same bwa index.
- **Junction evidence comes from discordant pairs, not from reads on the
  cassette's homology arms.** The arms are identical to the genome flanks, so
  reads there are multi-mappers at MAPQ 0 and get filtered. The informative read
  is the one whose mate sits in unique genome sequence.
- **Depth is measured twice per target, at `--min-mapq` and at MAPQ 0.** Zero
  depth at MAPQ 20 with healthy depth at MAPQ 0 is an unmappable repeat, not a
  deletion. Reported as `ambiguous_repeat`, never as a pass.
- **Median depth comes from a capped histogram**, not from sorting ~43M lines.
  Use a strict `>` on the running count: `>=` returns the lower median, which on
  a genome half-covered by zeros collapses the denominator to 0.
- **Module names are hints; the version assert is the check.** Sourcing
  `bwa-0.7.7` on this cluster has delivered bwa 0.7.17. Both scripts log the
  resolved path and version. Override with `BWA_MODULE` / `SAMTOOLS_MODULE`.
- **Targets are arguments, never inlined.** A typo'd contig name aborts against
  the `.fai` rather than silently yielding zero coverage and a confident
  "knockout confirmed".

## Known limits

- **No sequenced WT control.** A pre-existing repeat or a segmental duplication
  can produce an ectopic-looking cluster. Clusters are reported as
  `ectopic_candidate` and force a `REVIEW` verdict; they are not called.
- **Copy number is a depth ratio**, so it is sensitive to GC bias and to
  coverage that is uneven across the run. Treat ~1x as "consistent with single
  copy", not as proof of it.
- **Small CRISPR indels are out of scope.** This detects cassette replacement.
  A markerless edit needs variant calling instead.
- Bins are fixed 500 bp, so a junction falling on a bin boundary splits across
  two clusters. Both still classify correctly; the per-cluster counts halve.

## Superseded

`*_convert_fastq_to_bam.bash` in the Novogene run directories globbed `*.fq.gz`
and passed one file per `bwa mem` call, producing a separate single-end BAM for
R1 and for R2. No proper pairs, no insert sizes, and no way to compute the
mate-based junction evidence above. `ko_verify.sh` passes R1 and R2 to a single
`bwa mem` call and aborts if the properly-paired count is zero.
