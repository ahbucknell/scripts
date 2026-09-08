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
#     not proof. It also reads low if the cassette's homology arms are long,
#     since arm reads multi-map and are filtered.
#   - small markerless CRISPR indels are out of scope; that needs variant calling.
#
# Supersedes the *_convert_fastq_to_bam.bash scripts, which globbed *.fq.gz and
# passed one file per bwa mem call, making a separate single-end BAM for R1 and
# for R2 -- no proper pairs, so no mate-based junction evidence at all.

set -e

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*" >&2; }

# Tuning constants. Edit here rather than adding flags for them.
BIN=500        # junction clustering window, bp
FLANK=5000     # how far from the target a cluster still counts as its junction
MINJUNC=3      # reads needed before a cluster is believed
MINMAPQ=20     # below this a read is treated as ambiguously placed

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

# --- 3. cassette-to-genome junctions ----------------------------------------
# Two signals pooled: discordant pairs (read on the cassette, mate on the
# genome) and split reads (an SA:Z alignment crossing over). 0x90C drops
# secondary, supplementary, unmapped and mate-unmapped records.
#
# Note this deliberately does not use reads sitting on the cassette's homology
# arms: those are identical to the genome flanks, so they multi-map at MAPQ 0.
# The informative read is the one whose mate is in unique genome sequence.
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

JUNC="$WORK/${SAMPLE}.junctions.tsv"
awk -v OFS='\t' -v bin="$BIN" -v minj="$MINJUNC" -v flank="$FLANK" \
    -v tc="$TC" -v ts="$TS" -v te="$TE" '
  {
    b = int($2 / bin) * bin; k = $1 SUBSEP b
    total[k]++; contig[k] = $1; pos[k] = b
    if ($3 == "pair") npair[k]++; else nsplit[k]++
  }
  END {
    print "contig", "bin_start", "bin_end", "n_reads", "n_pair", "n_split", "class"
    for (k in total) {
      if (total[k] < minj) continue
      cls = "ectopic_candidate"
      if (contig[k] == tc && pos[k] + bin >= ts - flank && pos[k] <= te + flank)
        cls = "expected_junction"
      print contig[k], pos[k], pos[k] + bin, total[k], npair[k] + 0, nsplit[k] + 0, cls
    }
  }' "$TMP/links.tsv" \
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
if awk -v r="$TRATIO" -v z="$TZERO" 'BEGIN { exit !(r < 0.10 && z > 0.90) }'; then
  DROPOUT="yes"
  if awk -v q="$TMEAN_Q0" -v g="$GMED" 'BEGIN { exit !(q > 0.30 * g) }'; then
    DROPOUT="ambiguous_repeat"
  fi
fi

INTEGRATION="no"
if [ "$JL" -ge "$MINJUNC" ] && [ "$JR" -ge "$MINJUNC" ]; then
  INTEGRATION="both_flanks"
elif [ "$JL" -ge "$MINJUNC" ] || [ "$JR" -ge "$MINJUNC" ]; then
  INTEGRATION="one_flank_only"
fi

COPYNUM="single"
if awk -v r="$CRATIO" 'BEGIN { exit !(r > 1.5) }'; then
  COPYNUM="multi_copy"
elif awk -v r="$CRATIO" 'BEGIN { exit !(r < 0.60) }'; then
  COPYNUM="low_or_partial"
fi

VERDICT="PASS"
[ "$DROPOUT" = "yes" ] || VERDICT="FAIL"
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
