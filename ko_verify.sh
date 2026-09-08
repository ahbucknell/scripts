#!/bin/bash
# Verify a CRISPR + cassette knock-in strain from paired-end WGS.
#
# Maps R1/R2 as a PAIR to a chimeric reference (genome + cassette contig), then
# reports three independent lines of evidence:
#
#   1. dropout   -- depth across the target locus vs genome-wide median
#   2. junction  -- read pairs bridging the cassette to unique genome sequence,
#                   clustered by position; clusters at the target flanks are the
#                   expected integration junctions
#   3. ectopic   -- the same clusters anywhere else, plus cassette depth ratio
#                   as a copy-number estimate
#
# Usage:
# Build the chimeric reference first -- it is just a concatenation:
#   cat MGGv8_genome.fasta cassette.fasta > chimeric.fasta
#   bwa index chimeric.fasta && samtools faidx chimeric.fasta
#   grep '>' cassette.fasta        # note the name(s) for --cassette-contig
#
# Usage:
#   ko_verify.sh --sample NAME --r1 R1.fq.gz --r2 R2.fq.gz \
#                --ref chimeric.fasta --cassette-contig NAME \
#                --target contig:start-end --outdir DIR \
#                [--target contig:start-end ...] \
#                [--threads 4] [--min-mapq 20] [--flank-window 5000] \
#                [--min-junction-reads 3] [--arms arms.bed] [--keep-bam]
#
# Full run, for a batch of strains:
#   mkdir -p logs
#   cat MGGv8_genome.fasta cassette.fasta > chimeric.fasta
#   bwa index chimeric.fasta && samtools faidx chimeric.fasta
#   sbatch --array=1-$(($(wc -l < samples.tsv) - 1)) ko_verify_array.sbatch \
#          samples.tsv chimeric.fasta results HPH_CASSETTE
#
# A strain passes only if the target is deleted, BOTH flanks carry junction
# support, cassette copy number is ~1x, and no ectopic cluster is found.
#
# Limits, so the report is not over-read:
#   - No sequenced WT control, so a pre-existing repeat or segmental duplication
#     can look ectopic. Such clusters force REVIEW; they are not called.
#   - Copy number is a depth ratio and is sensitive to GC and coverage bias.
#     ~1x means "consistent with single copy", not proof of it.
#   - Small markerless CRISPR indels are out of scope; that needs variant
#     calling, not coverage.
#   - Bins are a fixed 500 bp, so a junction on a bin boundary splits across two
#     clusters. Both still classify correctly, but per-cluster counts halve.
#
# Supersedes the *_convert_fastq_to_bam.bash scripts in the Novogene run
# directories, which globbed *.fq.gz and passed ONE file per bwa mem call,
# producing a separate single-end BAM for R1 and for R2.

set -e

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*" >&2; }

SAMPLE=""; R1=""; R2=""; REF=""; OUTDIR=""; ARMS=""; CAS_ARG=""
THREADS=4; MINMAPQ=20; FLANK=5000; MINJUNC=3; BIN=500; KEEPBAM=0
TARGETS=()

while [ $# -gt 0 ]; do
  case "$1" in
    --sample)             SAMPLE="$2";  shift 2 ;;
    --r1)                 R1="$2";      shift 2 ;;
    --r2)                 R2="$2";      shift 2 ;;
    --ref)                REF="$2";     shift 2 ;;
    --outdir)             OUTDIR="$2";  shift 2 ;;
    --arms)               ARMS="$2";    shift 2 ;;
    --cassette-contig)    CAS_ARG="$2"; shift 2 ;;
    --target)             TARGETS+=("$2"); shift 2 ;;
    --threads)            THREADS="$2"; shift 2 ;;
    --min-mapq)           MINMAPQ="$2"; shift 2 ;;
    --flank-window)       FLANK="$2";   shift 2 ;;
    --min-junction-reads) MINJUNC="$2"; shift 2 ;;
    --bin)                BIN="$2";     shift 2 ;;
    --keep-bam)           KEEPBAM=1;    shift 1 ;;
    -h|--help)            sed -n '2,40p' "$0"; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

# --- validate every input BEFORE requesting real work ----------------------
# A tool that exits 0 having produced nothing useful is the failure mode that
# reaches results unnoticed, so the checks are front-loaded and the asserts at
# the end are on outputs.
[ -n "$SAMPLE" ] || die "--sample is required"
[ -n "$OUTDIR" ] || die "--outdir is required"
[ -s "$R1" ] || die "R1 missing or empty: $R1"
[ -s "$R2" ] || die "R2 missing or empty: $R2"
[ "$R1" != "$R2" ] || die "--r1 and --r2 are the same file: $R1"
[ -s "$REF" ] || die "chimeric reference missing or empty: $REF"
[ -s "${REF}.bwt" ] || die "reference is not bwa-indexed: ${REF}.bwt absent"
[ -s "${REF}.fai" ] || die "reference is not faidx-indexed: ${REF}.fai absent"
[ "${#TARGETS[@]}" -ge 1 ] || die "at least one --target contig:start-end is required"
[ -n "$CAS_ARG" ] || die "--cassette-contig is required (its name in the reference; see: grep '>' cassette.fasta)"

# Cassette contigs are matched by exact name, never by pattern: real cassette
# headers carry dots and dashes that would be live regex metacharacters.
CAS_NAMES=()
IFS=',' read -r -a CAS_SPLIT <<< "$CAS_ARG"
for c in "${CAS_SPLIT[@]}"; do
  [ -n "$c" ] || continue
  awk -v c="$c" '$1 == c { found = 1 } END { exit !found }' "${REF}.fai" \
    || die "cassette contig '$c' is not in ${REF}.fai -- check: grep '>' your cassette FASTA"
  CAS_NAMES+=("$c")
done
[ "${#CAS_NAMES[@]}" -ge 1 ] || die "--cassette-contig parsed to nothing: $CAS_ARG"
[ -z "$ARMS" ] || [ -s "$ARMS" ] || die "arms BED missing or empty: $ARMS"

case "$THREADS" in (*[!0-9]*|"") die "--threads must be an integer: $THREADS" ;; esac
case "$MINMAPQ" in (*[!0-9]*|"") die "--min-mapq must be an integer: $MINMAPQ" ;; esac
case "$MINJUNC" in (*[!0-9]*|"") die "--min-junction-reads must be an integer: $MINJUNC" ;; esac

# Reject a target whose contig is not in the reference, or whose coordinates are
# reversed / out of range. mmseqs-style silent garbage-in is the thing to avoid:
# a typo'd contig name must stop the job, not quietly yield zero coverage and a
# confident "knockout confirmed".
for t in "${TARGETS[@]}"; do
  echo "$t" | grep -Eq '^[^:]+:[0-9]+-[0-9]+$' \
    || die "target must look like contig:start-end, got: $t"
  tc="${t%%:*}"; span="${t#*:}"; ts="${span%%-*}"; te="${span##*-}"
  [ "$ts" -ge 1 ]   || die "target start must be >= 1: $t"
  [ "$te" -ge "$ts" ] || die "target end is before its start: $t"
  clen=$(awk -v c="$tc" '$1 == c { print $2; exit }' "${REF}.fai")
  [ -n "$clen" ] || die "contig '$tc' is not in ${REF}.fai (typo? wrong reference?)"
  [ "$te" -le "$clen" ] || die "target end $te exceeds length of $tc ($clen): $t"
  for c in "${CAS_NAMES[@]}"; do
    [ "$tc" != "$c" ] || die "--target points at the cassette contig '$c'; it must be the genome locus being replaced"
  done
done

BWA_MODULE="${BWA_MODULE:-bwa-0.7.17}"
SAMTOOLS_MODULE="${SAMTOOLS_MODULE:-samtools-1.10}"

set +eu
source "$BWA_MODULE"
source "$SAMTOOLS_MODULE"
set -e
command -v bwa      >/dev/null || die "bwa not on PATH after sourcing $BWA_MODULE"
command -v samtools >/dev/null || die "samtools not on PATH after sourcing $SAMTOOLS_MODULE"
log "bwa:      $(command -v bwa)      $(bwa 2>&1 | awk '/^Version/ { print $2; exit }')"
log "samtools: $(command -v samtools) $(samtools --version | awk 'NR == 1 { print $2 }')"

mkdir -p "$OUTDIR"
WORK="$OUTDIR/$SAMPLE"; mkdir -p "$WORK"
TMP="${SLURM_TMPDIR:-$WORK}/tmp.$$"; mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

BAM="$WORK/${SAMPLE}.sorted.bam"

CAS_LIST="$TMP/cassette_contigs.txt"
printf '%s\n' "${CAS_NAMES[@]}" > "$CAS_LIST"
CAS_CONTIGS=$(cat "$CAS_LIST")
log "cassette contigs: $(printf '%s,' "${CAS_NAMES[@]}" | sed 's/,$//')"

# ---------------------------------------------------------------------------
# 1. Align R1 and R2 together, as a pair.
# ---------------------------------------------------------------------------
# This is the bug in the original convert_fastq_to_bam scripts: they globbed
# *.fq.gz and passed ONE file per bwa call, so R1 and R2 became two single-end
# BAMs. No proper-pair flags, no insert sizes -- and the mate-based junction
# evidence below is impossible to compute from those.
if [ -s "$BAM" ] && samtools quickcheck "$BAM" 2>/dev/null; then
  log "reusing existing BAM: $BAM"
else
  log "aligning $SAMPLE (paired) with $THREADS threads"
  # pipefail is on ONLY here: without it a bwa failure mid-stream leaves a
  # truncated-but-sorted BAM and samtools still exits 0.
  set -o pipefail
  bwa mem -t "$THREADS" \
      -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
      "$REF" "$R1" "$R2" \
    | samtools sort -@ "$THREADS" -T "$TMP/sort" -o "$BAM" -
  set +o pipefail
  samtools index -@ "$THREADS" "$BAM"
fi

samtools quickcheck "$BAM" || die "BAM failed samtools quickcheck: $BAM"
[ -s "${BAM}.bai" ] || samtools index -@ "$THREADS" "$BAM"

samtools flagstat -@ "$THREADS" "$BAM" > "$WORK/${SAMPLE}.flagstat.txt"

# Assert the pairing actually happened. If this is 0 the alignment was
# effectively single-end and every junction number below is meaningless.
PAIRED=$(awk '/properly paired/ { gsub(/[()]/, "", $1); print $1; exit }' "$WORK/${SAMPLE}.flagstat.txt")
MAPPED=$(awk '/ mapped \(/ { print $1; exit }' "$WORK/${SAMPLE}.flagstat.txt")
[ "${PAIRED:-0}" -gt 0 ] || die "zero properly-paired reads -- R1/R2 were not aligned as a pair"
log "mapped=$MAPPED properly_paired=$PAIRED"

# ---------------------------------------------------------------------------
# 2. Depth statistics
# ---------------------------------------------------------------------------
# Median comes from a depth histogram rather than sorting ~43M lines.
depth_stats() {  # reads `samtools depth -a` on stdin -> "mean median n zero_frac"
  awk '
    { d = $3; sum += d; n++; if (d > 1000) d = 1000; h[d]++ }
    END {
      if (n == 0) { print "0 0 0 1"; exit }
      half = n / 2; c = 0; med = 0
      # strict > gives the upper median; >= would return the lower one, which
      # for a genome half-covered by zeros collapses the denominator to 0.
      for (d = 0; d <= 1000; d++) { c += h[d]; if (c > half) { med = d; break } }
      printf "%.3f %d %d %.4f\n", sum / n, med, n, h[0] / n
    }'
}

log "computing genome-wide depth (cassette contigs excluded)"
GENOME_STATS=$(samtools depth -a -Q "$MINMAPQ" "$BAM" \
  | awk 'NR == FNR { cas[$1]; next } !($1 in cas)' "$CAS_LIST" - | depth_stats)
GENOME_MEDIAN=$(echo "$GENOME_STATS" | awk '{ print $2 }')
GENOME_MEAN=$(echo "$GENOME_STATS" | awk '{ print $1 }')
[ "$GENOME_MEDIAN" -gt 0 ] || die "genome-wide median depth is 0 -- alignment or reference is wrong"
log "genome median depth = ${GENOME_MEDIAN}x (mean ${GENOME_MEAN}x)"

# Cassette depth, optionally restricted to the marker core so that homology
# arms -- which are identical to the genome flanks and therefore multi-map --
# do not drag the copy-number estimate down.
if [ -n "$ARMS" ]; then
  log "cassette depth restricted to regions in $ARMS"
  CAS_STATS=$(samtools depth -a -Q "$MINMAPQ" -b "$ARMS" "$BAM" \
    | awk 'NR == FNR { cas[$1]; next } ($1 in cas)' "$CAS_LIST" - | depth_stats)
else
  # Query the cassette contigs by region so this is an indexed lookup rather
  # than a second scan of the whole BAM.
  CAS_STATS=$(for c in $CAS_CONTIGS; do
                samtools depth -a -Q "$MINMAPQ" -r "$c" "$BAM"
              done | depth_stats)
fi
CAS_MEDIAN=$(echo "$CAS_STATS" | awk '{ print $2 }')
CAS_RATIO=$(awk -v c="$CAS_MEDIAN" -v g="$GENOME_MEDIAN" 'BEGIN { printf "%.2f", c / g }')
log "cassette median depth = ${CAS_MEDIAN}x  ratio to genome = ${CAS_RATIO}"

# ---------------------------------------------------------------------------
# 3. Junction evidence: cassette-to-genome links
# ---------------------------------------------------------------------------
# Two independent signals, pooled:
#   (a) discordant pairs -- read on the cassette, mate on a genome contig
#   (b) split reads      -- an SA:Z supplementary alignment crossing over
# 0x90C drops secondary, supplementary, unmapped and mate-unmapped records so
# (a) counts only real cross-contig pairs.
log "extracting cassette-to-genome junction reads (MAPQ >= $MINMAPQ)"
: > "$TMP/links.tsv"

samtools view -q "$MINMAPQ" -F 0x90C "$BAM" $CAS_CONTIGS \
  | awk -v OFS='\t' '
      NR == FNR { cas[$1]; next }
      $7 != "=" && $7 != "*" && !($7 in cas) { print $7, $8, "pair" }' "$CAS_LIST" - \
  >> "$TMP/links.tsv"

samtools view -q "$MINMAPQ" -F 0x100 "$BAM" $CAS_CONTIGS \
  | awk -v OFS='\t' '
      NR == FNR { cas[$1]; next }
      {
        for (i = 12; i <= NF; i++) {
          if ($i ~ /^SA:Z:/) {
            n = split(substr($i, 6), recs, ";")
            for (j = 1; j <= n; j++) {
              if (recs[j] == "") continue
              split(recs[j], f, ",")
              if (!(f[1] in cas)) print f[1], f[2], "split"
            }
          }
        }
      }' "$CAS_LIST" - \
  >> "$TMP/links.tsv"

NLINKS=$(wc -l < "$TMP/links.tsv" | tr -d ' ')
log "cassette-to-genome links: $NLINKS"

# Cluster into fixed bins and classify against the target flanks.
printf '%s\n' "${TARGETS[@]}" \
  | awk -F'[:-]' -v OFS='\t' '{ print $1, $2, $3 }' > "$TMP/targets.tsv"

awk -v OFS='\t' -v bin="$BIN" -v flank="$FLANK" -v minj="$MINJUNC" '
  # first file: targets
  NR == FNR { tc[FNR] = $1; ts[FNR] = $2; te[FNR] = $3; nt = FNR; next }
  # second file: links -> bin them
  {
    b = int($2 / bin) * bin
    key = $1 SUBSEP b
    total[key]++
    if ($3 == "pair") npair[key]++; else nsplit[key]++
    contig[key] = $1; pos[key] = b
  }
  END {
    print "contig", "bin_start", "bin_end", "n_reads", "n_pair", "n_split", "class"
    for (k in total) {
      if (total[k] < minj) continue
      cls = "ectopic_candidate"
      for (i = 1; i <= nt; i++) {
        if (contig[k] == tc[i] && pos[k] + bin >= ts[i] - flank && pos[k] <= te[i] + flank) {
          cls = "expected_junction"
          break
        }
      }
      print contig[k], pos[k], pos[k] + bin, total[k], npair[k] + 0, nsplit[k] + 0, cls
    }
  }' "$TMP/targets.tsv" "$TMP/links.tsv" \
  | { read -r hdr; echo "$hdr"; sort -k4,4nr; } > "$WORK/${SAMPLE}.junctions.tsv"

N_EXPECTED=$(awk 'NR > 1 && $7 == "expected_junction"' "$WORK/${SAMPLE}.junctions.tsv" | wc -l | tr -d ' ')
N_ECTOPIC=$(awk 'NR > 1 && $7 == "ectopic_candidate"'  "$WORK/${SAMPLE}.junctions.tsv" | wc -l | tr -d ' ')
log "junction clusters: expected=$N_EXPECTED ectopic_candidates=$N_ECTOPIC"

# ---------------------------------------------------------------------------
# 4. Per-target report
# ---------------------------------------------------------------------------
REPORT="$WORK/${SAMPLE}.ko_report.tsv"
{
  printf 'sample\ttarget\ttarget_len\tgenome_median_depth\ttarget_mean_depth\t'
  printf 'target_mean_depth_mapq0\ttarget_zero_frac\tdepth_ratio\tcassette_depth_ratio\t'
  printf 'junctions_left\tjunctions_right\tectopic_clusters\tdropout\tintegration\tcopy_number\tverdict\n'
} > "$REPORT"

OVERALL="PASS"
for t in "${TARGETS[@]}"; do
  tc="${t%%:*}"; span="${t#*:}"; ts="${span%%-*}"; te="${span##*-}"
  TLEN=$((te - ts + 1))

  T_STATS=$(samtools depth -a -Q "$MINMAPQ" -r "$t" "$BAM" | depth_stats)
  T_MEAN=$(echo "$T_STATS" | awk '{ print $1 }')
  T_ZERO=$(echo "$T_STATS" | awk '{ print $4 }')

  # MAPQ 0 depth too: if the locus looks empty at MAPQ 20 but full at MAPQ 0,
  # it is a repeat/mappability artefact, not a deletion.
  T_MEAN_Q0=$(samtools depth -a -Q 0 -r "$t" "$BAM" | depth_stats | awk '{ print $1 }')

  RATIO=$(awk -v a="$T_MEAN" -v b="$GENOME_MEDIAN" 'BEGIN { printf "%.3f", a / b }')

  # Junction support on each side of the target separately -- a correct
  # single-locus integration should show both, not just one.
  # A cluster counts for the LEFT flank only if it ends at or before the target
  # start, and for the RIGHT flank only if it begins at or after the target end.
  # Overlapping windows here would let a single-sided integration report support
  # on both flanks -- exactly the false pass this script exists to catch. The
  # flank distance limit is already applied by the expected_junction class.
  JL=$(awk -v c="$tc" -v s="$ts" \
        'NR > 1 && $7 == "expected_junction" && $1 == c && $3 <= s { n += $4 } END { print n + 0 }' \
        "$WORK/${SAMPLE}.junctions.tsv")
  JR=$(awk -v c="$tc" -v e="$te" \
        'NR > 1 && $7 == "expected_junction" && $1 == c && $2 >= e { n += $4 } END { print n + 0 }' \
        "$WORK/${SAMPLE}.junctions.tsv")

  # Every test below is written as `if`, not `test && VAR=x`. Under set -e a
  # false `a && b` chain is a failing command and would kill the run silently.
  DROPOUT="no"
  if awk -v r="$RATIO" -v z="$T_ZERO" 'BEGIN { exit !(r < 0.10 && z > 0.90) }'; then
    DROPOUT="yes"
    # Distinguish a real deletion from a locus that is merely unmappable: if
    # MAPQ-0 depth is healthy where MAPQ-filtered depth is zero, the "dropout"
    # is a repeat artefact and must not be reported as a knockout.
    if awk -v q0="$T_MEAN_Q0" -v g="$GENOME_MEDIAN" 'BEGIN { exit !(q0 > 0.30 * g) }'; then
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
  if awk -v r="$CAS_RATIO" 'BEGIN { exit !(r > 1.5) }'; then
    COPYNUM="multi_copy"
  elif awk -v r="$CAS_RATIO" 'BEGIN { exit !(r < 0.60) }'; then
    COPYNUM="low_or_partial"
  fi

  VERDICT="PASS"
  [ "$DROPOUT" = "yes" ] || VERDICT="FAIL"
  [ "$INTEGRATION" = "both_flanks" ] || VERDICT="FAIL"
  [ "$COPYNUM" = "single" ] || VERDICT="REVIEW"
  [ "$N_ECTOPIC" -eq 0 ] || VERDICT="REVIEW"
  [ "$VERDICT" = "PASS" ] || OVERALL="$VERDICT"

  printf '%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n' \
    "$SAMPLE" "$t" "$TLEN" "$GENOME_MEDIAN" "$T_MEAN" "$T_MEAN_Q0" "$T_ZERO" \
    "$RATIO" "$CAS_RATIO" "$JL" "$JR" "$N_ECTOPIC" \
    "$DROPOUT" "$INTEGRATION" "$COPYNUM" "$VERDICT" >> "$REPORT"
done

# Per-base depth over each target, for plotting.
: > "$WORK/${SAMPLE}.target_depth.tsv"
for t in "${TARGETS[@]}"; do
  samtools depth -a -Q "$MINMAPQ" -r "$t" "$BAM" >> "$WORK/${SAMPLE}.target_depth.tsv"
done

# --- assert the outputs exist and are non-trivial ---------------------------
[ -s "$REPORT" ] || die "report was not written: $REPORT"
[ "$(wc -l < "$REPORT")" -gt 1 ] || die "report has a header but no rows: $REPORT"

[ "$KEEPBAM" -eq 1 ] || { rm -f "$BAM" "${BAM}.bai"; log "removed BAM (pass --keep-bam to keep it)"; }

log "verdict for $SAMPLE: $OVERALL"
column -t "$REPORT" >&2 || cat "$REPORT" >&2
echo "$REPORT"
