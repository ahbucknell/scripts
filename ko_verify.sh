#!/bin/bash
#SBATCH --job-name=ko_verify
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=ko_verify_%j.out
#SBATCH --error=ko_verify_%j.err
#
# Verify one CRISPR + cassette knock-in strain at one target locus, from
# paired-end WGS. Maps R1/R2 as a pair against the genome and the cassette
# together, then reports three things from that single alignment:
#
#   dropout   depth across the target vs genome-wide median (a correct
#             replacement leaves the target at ~0x)
#   junction  read pairs with one mate on the cassette and the other in unique
#             genome sequence, clustered by position -- clusters at the target
#             flanks are the expected integration junctions
#   ectopic   the same clusters landing anywhere else, plus cassette depth
#             divided by genome depth as a copy-number estimate
#
# PASS requires all four: target deleted, BOTH flanks supported, copy number
# ~1x, no ectopic cluster.
#
# Build the reference once -- it is just a concatenation:
#   cat MGGv8_genome.fasta cassette.fasta > chimeric.fasta
#   bwa index chimeric.fasta && samtools faidx chimeric.fasta
#   grep '>' cassette.fasta        # this name is --cassette-contig
#
# Usage:
#   ko_verify.sh --sample NAME --r1 R1.fq.gz --r2 R2.fq.gz \
#                --ref chimeric.fasta --cassette-contig NAME \
#                --target contig:start-end --outdir DIR \
#                [--threads N] [--keep-bam]
#
# On the cluster, one strain per submission:
#   sbatch ko_verify.sh --sample KO_A --r1 ... --r2 ... --ref chimeric.fasta \
#          --cassette-contig HPH_CASSETTE --target contig1:50000-52000 \
#          --outdir results
#
# Reading the result:
#   - dropout=ambiguous_repeat means the locus is unmappable, not deleted.
#   - ectopic clusters force REVIEW, they are not called: with no sequenced WT
#     control a pre-existing repeat looks the same.
#   - copy number is a depth ratio, so ~1x means "consistent with single copy",
#     not proof of it.
#   - small markerless CRISPR indels are out of scope; that needs variant calling.
#
# Supersedes the *_convert_fastq_to_bam.bash scripts, which globbed *.fq.gz and
# passed one file per bwa mem call, making a separate single-end BAM for R1 and
# for R2 -- no proper pairs, so no mate-based junction evidence at all.

set -e

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*" >&2; }

# --- tuning constants -------------------------------------------------------
# Every threshold lives here, including the ones the verdict logic applies.
# Edit these rather than adding flags.

MINMAPQ=20      # MAPQ 20 = 1% chance the read is misplaced. bwa mem gives 0 to
                # equal-best multi-mappers and scales to 60. Field standard.

GAP=500         # Junction clustering: positions further apart than this start a
                # new cluster. Mates from one junction scatter over about one
                # insert (~350 bp for Novogene PE150), so this is ~1.5 inserts.

FLANK=1000      # How far from the target edge a cluster still counts as that
                # target's junction. Derived, not guessed: a mate cannot land
                # further out than homology arm (70 bp) + insert (~350 bp), so
                # 1000 is that with slack. Raising this hides real ectopic
                # insertions by absorbing them into "expected".

MINJUNC_FLOOR=5 # Junction reads needed before a cluster is believed. Scaled
MINJUNC_FRAC=0.15  # against measured depth: at 100x expect 20-35 pairs per
                # flank, so a fixed floor of 3 would accept noise. On a poor
                # library the floor takes over. Computed once GMED is known.

DROP_RATIO=0.10 # Target counts as deleted below this fraction of genome median
DROP_ZERO=0.90  # ...and with at least this fraction of its bases at zero depth
PART_ZERO=0.30  # Above this but short of DROP_ZERO is a partial deletion
REPEAT_FRAC=0.30   # MAPQ-0 depth above this fraction of genome median means the
                # "dropout" is an unmappable repeat, not a deletion
CN_HIGH=1.5     # Cassette depth ratio above this suggests an extra copy
CN_LOW=0.60     # ...and below this a partial or truncated integration
CN_ABSENT=0.10  # ...and below this there is no cassette in the strain at all

SAMPLE=""; R1=""; R2=""; REF=""; CAS=""; TARGET=""; OUTDIR=""
THREADS="${SLURM_CPUS_PER_TASK:-4}"; KEEPBAM=0

while [ $# -gt 0 ]; do
  case "$1" in
    --sample)          SAMPLE="$2";  shift 2 ;;
    --r1)              R1="$2";      shift 2 ;;
    --r2)              R2="$2";      shift 2 ;;
    --ref)             REF="$2";     shift 2 ;;
    --cassette-contig) CAS="$2";     shift 2 ;;
    --target)          TARGET="$2";  shift 2 ;;
    --outdir)          OUTDIR="$2";  shift 2 ;;
    --threads)         THREADS="$2"; shift 2 ;;
    --keep-bam)        KEEPBAM=1;    shift 1 ;;
    -h|--help)         sed -n '11,45p' "$0"; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

# --- validate everything before asking for real work ------------------------
# A tool that exits 0 having produced nothing useful is the failure mode that
# reaches results unnoticed, so inputs are checked here and outputs asserted
# at the end.
[ -n "$SAMPLE" ] || die "--sample is required"
[ -n "$OUTDIR" ] || die "--outdir is required"
[ -s "$R1" ] || die "R1 missing or empty: $R1"
[ -s "$R2" ] || die "R2 missing or empty: $R2"
[ "$R1" != "$R2" ] || die "--r1 and --r2 are the same file: $R1"
[ -s "$REF" ] || die "reference missing or empty: $REF"
[ -s "${REF}.bwt" ] || die "reference is not bwa-indexed: ${REF}.bwt absent"
[ -s "${REF}.fai" ] || die "reference is not faidx-indexed: ${REF}.fai absent"
[ -n "$CAS" ] || die "--cassette-contig is required (see: grep '>' cassette.fasta)"
[ -n "$TARGET" ] || die "--target contig:start-end is required"
case "$THREADS" in (*[!0-9]*|"") die "--threads must be an integer: $THREADS" ;; esac

awk -v c="$CAS" '$1 == c { found = 1 } END { exit !found }' "${REF}.fai" \
  || die "cassette contig '$CAS' is not in ${REF}.fai -- check: grep '>' your cassette FASTA"

# A typo'd contig must stop the job, not quietly return zero coverage and a
# confident "knockout confirmed".
echo "$TARGET" | grep -Eq '^[^:]+:[0-9]+-[0-9]+$' \
  || die "--target must look like contig:start-end, got: $TARGET"
TC="${TARGET%%:*}"; SPAN="${TARGET#*:}"; TS="${SPAN%%-*}"; TE="${SPAN##*-}"
[ "$TS" -ge 1 ] || die "target start must be >= 1: $TARGET"
[ "$TE" -ge "$TS" ] || die "target end is before its start: $TARGET"
[ "$TC" != "$CAS" ] || die "--target points at the cassette; it must be the genome locus being replaced"
CLEN=$(awk -v c="$TC" '$1 == c { print $2; exit }' "${REF}.fai")
[ -n "$CLEN" ] || die "contig '$TC' is not in ${REF}.fai (typo? wrong reference?)"
[ "$TE" -le "$CLEN" ] || die "target end $TE exceeds length of $TC ($CLEN)"

# --- tools ------------------------------------------------------------------
# Module names are a hint, not a guarantee: sourcing bwa-0.7.7 on this cluster
# has delivered bwa 0.7.17. Check activation by its outcome and log what we got.
BWA_MODULE="${BWA_MODULE:-bwa-0.7.17}"
SAMTOOLS_MODULE="${SAMTOOLS_MODULE:-samtools-1.10}"

# Activation scripts are not written for strict mode: the TSL wrapper reads
# $WRAPPER_DEBUG unset and contains a harmlessly-failing command. "source X ||
# true" does not contain this -- errexit still fires inside the sourced file.
set +eu
source "$BWA_MODULE"
source "$SAMTOOLS_MODULE"
set -e
command -v bwa      >/dev/null || die "bwa not on PATH after sourcing $BWA_MODULE"
command -v samtools >/dev/null || die "samtools not on PATH after sourcing $SAMTOOLS_MODULE"
log "bwa:      $(command -v bwa)      $(bwa 2>&1 | awk '/^Version/ { print $2; exit }')"
log "samtools: $(command -v samtools) $(samtools --version | awk 'NR == 1 { print $2 }')"

WORK="$OUTDIR/$SAMPLE"; mkdir -p "$WORK"
TMP="$WORK/tmp.$$"; mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT
BAM="$WORK/${SAMPLE}.sorted.bam"

# --- 1. align, as a pair ----------------------------------------------------
if [ -s "$BAM" ] && samtools quickcheck "$BAM" 2>/dev/null; then
  log "reusing existing BAM: $BAM"
else
  log "aligning $SAMPLE (paired, $THREADS threads)"
  # pipefail only here: without it a bwa failure mid-stream leaves a truncated
  # but validly-sorted BAM and samtools still exits 0.
  set -o pipefail
  bwa mem -t "$THREADS" -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
      "$REF" "$R1" "$R2" \
    | samtools sort -@ "$THREADS" -T "$TMP/sort" -o "$BAM" -
  set +o pipefail
  samtools index -@ "$THREADS" "$BAM"
fi
samtools quickcheck "$BAM" || die "BAM failed samtools quickcheck: $BAM"
[ -s "${BAM}.bai" ] || samtools index -@ "$THREADS" "$BAM"

samtools flagstat -@ "$THREADS" "$BAM" > "$WORK/${SAMPLE}.flagstat.txt"
PAIRED=$(awk '/properly paired/ { print $1; exit }' "$WORK/${SAMPLE}.flagstat.txt")
# If this is zero the run was effectively single-end and every junction number
# below is meaningless -- which is exactly how the old scripts failed.
[ "${PAIRED:-0}" -gt 0 ] || die "zero properly-paired reads -- R1/R2 were not aligned as a pair"
log "properly paired: $PAIRED"

# --- 2. depth ---------------------------------------------------------------
# Median from a capped histogram rather than sorting ~43M lines. The comparison
# must be strict >: with >= a genome half-covered by zeros returns median 0 and
# collapses every ratio below.
depth_stats() {  # `samtools depth -a` on stdin -> "mean median n zero_frac"
  awk '
    { d = $3; sum += d; n++; if (d > 1000) d = 1000; h[d]++ }
    END {
      if (n == 0) { print "0 0 0 1"; exit }
      half = n / 2; c = 0; med = 0
      for (d = 0; d <= 1000; d++) { c += h[d]; if (c > half) { med = d; break } }
      printf "%.3f %d %d %.4f\n", sum / n, med, n, h[0] / n
    }'
}

log "computing genome-wide depth"
GMED=$(samtools depth -a -Q "$MINMAPQ" "$BAM" \
  | awk -v c="$CAS" '$1 != c' | depth_stats | awk '{ print $2 }')
[ "$GMED" -gt 0 ] || die "genome-wide median depth is 0 -- alignment or reference is wrong"

CMED=$(samtools depth -a -Q "$MINMAPQ" -r "$CAS" "$BAM" | depth_stats | awk '{ print $2 }')
CRATIO=$(awk -v c="$CMED" -v g="$GMED" 'BEGIN { printf "%.2f", c / g }')

TSTAT=$(samtools depth -a -Q "$MINMAPQ" -r "$TARGET" "$BAM" | depth_stats)
TMEAN=$(echo "$TSTAT" | awk '{ print $1 }')
TZERO=$(echo "$TSTAT" | awk '{ print $4 }')
# Also at MAPQ 0: healthy depth here but none above suggests an unmappable
# repeat rather than a deletion.
TMEAN_Q0=$(samtools depth -a -Q 0 -r "$TARGET" "$BAM" | depth_stats | awk '{ print $1 }')
TRATIO=$(awk -v a="$TMEAN" -v g="$GMED" 'BEGIN { printf "%.3f", a / g }')
log "genome median ${GMED}x  target mean ${TMEAN}x  cassette ratio ${CRATIO}"

# Junction threshold scales with the depth actually observed, so the script does
# not assume a coverage level.
MINJUNC=$(awk -v g="$GMED" -v f="$MINJUNC_FRAC" -v fl="$MINJUNC_FLOOR" \
  'BEGIN { m = int(g * f); print (m > fl ? m : fl) }')
log "junction threshold: $MINJUNC reads (floor $MINJUNC_FLOOR, ${MINJUNC_FRAC} x depth)"

# --- 3. cassette-to-genome junctions ----------------------------------------
# Two signals pooled: discordant pairs (read on the cassette, mate on the
# genome) and split reads (an SA:Z alignment crossing over). 0x90C drops
# secondary, supplementary, unmapped and mate-unmapped records.
#
# The 70 bp homology arms are shorter than one read, so a read overlapping an
# arm almost always extends into unique genome or into the marker and is
# uniquely placeable. Long arms would break this: with arms longer than the
# insert size no pair can reach from the marker into unique genome, correct
# targeting creates no novel genome junction at all, and this whole check
# becomes unsatisfiable. Re-derive FLANK below if the construct changes.
log "extracting junction reads"
{
  samtools view -q "$MINMAPQ" -F 0x90C "$BAM" "$CAS" \
    | awk -v OFS='\t' -v c="$CAS" '$7 != "=" && $7 != "*" && $7 != c { print $7, $8, "pair" }'
  samtools view -q "$MINMAPQ" -F 0x100 "$BAM" "$CAS" \
    | awk -v OFS='\t' -v c="$CAS" '
        {
          for (i = 12; i <= NF; i++)
            if ($i ~ /^SA:Z:/) {
              n = split(substr($i, 6), recs, ";")
              for (j = 1; j <= n; j++) {
                if (recs[j] == "") continue
                split(recs[j], f, ",")
                if (f[1] != c) print f[1], f[2], "split"
              }
            }
        }'
} > "$TMP/links.tsv"

# Single-linkage clustering on a gap, not fixed bins: a fixed grid splits any
# junction that lands on a boundary across two clusters, halving both counts and
# potentially dropping both below threshold. A gap has no boundary to land on.
JUNC="$WORK/${SAMPLE}.junctions.tsv"
sort -k1,1 -k2,2n "$TMP/links.tsv" \
  | awk -v OFS='\t' -v gap="$GAP" -v minj="$MINJUNC" -v flank="$FLANK" \
        -v tc="$TC" -v ts="$TS" -v te="$TE" '
      function flush(   cls) {
        if (n < minj) return
        cls = "ectopic_candidate"
        if (c == tc && end >= ts - flank && start <= te + flank) cls = "expected_junction"
        print c, start, end, n, np + 0, ns + 0, cls
      }
      BEGIN { print "contig", "start", "end", "n_reads", "n_pair", "n_split", "class" }
      # n/np/ns must be initialised here too, not only in the reset branch: an
      # unset awk variable prints as "" rather than 0, which drops a field and
      # shifts every column of the first cluster left by one.
      NR == 1 { c = $1; start = $2; end = $2; n = 0; np = 0; ns = 0 }
      $1 != c || $2 - end > gap { flush(); c = $1; start = $2; n = 0; np = 0; ns = 0 }
      { end = $2; n++; if ($3 == "pair") np++; else ns++ }
      END { if (NR > 0) flush() }' \
  | { read -r h; echo "$h"; sort -k4,4nr; } > "$JUNC"

# A cluster counts for the LEFT flank only if it ends at or before the target
# start, and RIGHT only if it starts at or after the target end. Overlapping
# windows would let a one-sided integration claim support on both flanks --
# exactly the false pass this script exists to catch.
JL=$(awk -v c="$TC" -v s="$TS" 'NR > 1 && $7 == "expected_junction" && $1 == c && $3 <= s { n += $4 } END { print n + 0 }' "$JUNC")
JR=$(awk -v c="$TC" -v e="$TE" 'NR > 1 && $7 == "expected_junction" && $1 == c && $2 >= e { n += $4 } END { print n + 0 }' "$JUNC")
NECT=$(awk 'NR > 1 && $7 == "ectopic_candidate"' "$JUNC" | wc -l | tr -d ' ')
log "junction reads: left=$JL right=$JR  ectopic clusters=$NECT"

# --- 4. verdict -------------------------------------------------------------
DROPOUT="no"
if awk -v r="$TRATIO" -v z="$TZERO" -v dr="$DROP_RATIO" -v dz="$DROP_ZERO" \
     'BEGIN { exit !(r < dr && z > dz) }'; then
  DROPOUT="yes"
  # Healthy depth at MAPQ 0 where there is none above it means the locus is
  # unmappable, not deleted. Never report that as a knockout.
  if awk -v q="$TMEAN_Q0" -v g="$GMED" -v f="$REPEAT_FRAC" \
       'BEGIN { exit !(q > f * g) }'; then
    DROPOUT="ambiguous_repeat"
  fi
elif awk -v z="$TZERO" -v pz="$PART_ZERO" 'BEGIN { exit !(z > pz) }'; then
  # Part of the ORF is gone but not all of it. Distinguishing this from "no
  # deletion at all" matters -- both used to report the same FAIL.
  DROPOUT="partial"
fi

INTEGRATION="no"
if [ "$JL" -ge "$MINJUNC" ] && [ "$JR" -ge "$MINJUNC" ]; then
  INTEGRATION="both_flanks"
elif [ "$JL" -ge "$MINJUNC" ] || [ "$JR" -ge "$MINJUNC" ]; then
  INTEGRATION="one_flank_only"
fi

COPYNUM="single"
if awk -v r="$CRATIO" -v h="$CN_HIGH" 'BEGIN { exit !(r > h) }'; then
  COPYNUM="multi_copy"
elif awk -v r="$CRATIO" -v a="$CN_ABSENT" 'BEGIN { exit !(r < a) }'; then
  COPYNUM="absent"          # no cassette at all: untransformed, or wrong strain
elif awk -v r="$CRATIO" -v l="$CN_LOW" 'BEGIN { exit !(r < l) }'; then
  COPYNUM="low_or_partial"
fi

VERDICT="PASS"
case "$DROPOUT" in
  yes)     ;;
  partial) VERDICT="REVIEW" ;;
  *)       VERDICT="FAIL" ;;
esac
[ "$INTEGRATION" = "both_flanks" ] || VERDICT="FAIL"
if [ "$VERDICT" = "PASS" ]; then
  [ "$COPYNUM" = "single" ] || VERDICT="REVIEW"
  [ "$NECT" -eq 0 ] || VERDICT="REVIEW"
fi

REPORT="$WORK/${SAMPLE}.ko_report.tsv"
{
  printf 'sample\ttarget\tgenome_median_depth\ttarget_mean_depth\ttarget_mean_depth_mapq0\t'
  printf 'target_zero_frac\tdepth_ratio\tcassette_depth_ratio\tjunctions_left\t'
  printf 'junctions_right\tectopic_clusters\tdropout\tintegration\tcopy_number\tverdict\n'
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n' \
    "$SAMPLE" "$TARGET" "$GMED" "$TMEAN" "$TMEAN_Q0" "$TZERO" "$TRATIO" "$CRATIO" \
    "$JL" "$JR" "$NECT" "$DROPOUT" "$INTEGRATION" "$COPYNUM" "$VERDICT"
} > "$REPORT"

samtools depth -a -Q "$MINMAPQ" -r "$TARGET" "$BAM" > "$WORK/${SAMPLE}.target_depth.tsv"

# --- assert on the outputs --------------------------------------------------
[ "$(wc -l < "$REPORT")" -eq 2 ] || die "report is malformed: $REPORT"
[ -s "$WORK/${SAMPLE}.target_depth.tsv" ] || die "target depth file is empty -- was the region valid?"

[ "$KEEPBAM" -eq 1 ] || { rm -f "$BAM" "${BAM}.bai"; log "removed BAM (--keep-bam to retain)"; }

log "verdict for $SAMPLE: $VERDICT"
column -t "$REPORT" >&2 2>/dev/null || cat "$REPORT" >&2
echo "$REPORT"
