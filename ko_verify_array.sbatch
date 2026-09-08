#!/bin/bash
#SBATCH --job-name=ko_verify
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/ko_verify_%A_%a.out
#SBATCH --error=logs/ko_verify_%A_%a.err

# Run ko_verify.sh across a samplesheet, one strain per array task.
#
#   mkdir -p logs
#   ko_prepare_ref.sh --genome MGGv8_genome.fasta --cassette cassette.fasta \
#                     --outdir /path/to/chimeric_ref
#   sbatch --array=1-$(($(wc -l < samples.tsv) - 1)) ko_verify_array.sbatch \
#          samples.tsv /path/to/chimeric_ref/chimeric.fasta results
#
# samples.tsv is TAB-separated WITH a header line:
#   sample<TAB>r1<TAB>r2<TAB>targets
# where `targets` is one or more contig:start-end, comma-separated.
#
# The resource requests above are guesses, not measurements. Check a finished
# job with `seff $SLURM_JOB_ID` and tighten them.

set -e

die() { echo "ERROR: $*" >&2; exit 1; }

SAMPLESHEET="${1:?usage: sbatch --array=N-M ko_verify_array.sbatch samples.tsv chimeric.fasta outdir}"
REF="${2:?chimeric.fasta is required (see ko_prepare_ref.sh)}"
OUTDIR="${3:?output directory is required}"

[ -s "$SAMPLESHEET" ] || die "samplesheet missing or empty: $SAMPLESHEET"
[ -s "$REF" ]         || die "chimeric reference missing or empty: $REF"
[ -n "${SLURM_ARRAY_TASK_ID:-}" ] || die "not running as a job array (pass --array to sbatch)"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
VERIFY="$SCRIPT_DIR/ko_verify.sh"
[ -x "$VERIFY" ] || die "ko_verify.sh not found or not executable next to this script: $VERIFY"

# +1 skips the header row, so --array=1 is the first sample.
ROW=$((SLURM_ARRAY_TASK_ID + 1))
LINE=$(awk -v r="$ROW" 'NR == r' "$SAMPLESHEET")
[ -n "$LINE" ] || die "samplesheet has no row $ROW (array index out of range)"

SAMPLE=$(echo "$LINE" | cut -f1)
R1=$(echo "$LINE"     | cut -f2)
R2=$(echo "$LINE"     | cut -f3)
TARGETS=$(echo "$LINE" | cut -f4)

[ -n "$SAMPLE" ]  || die "row $ROW has no sample name (is the file TAB-separated?)"
[ -n "$TARGETS" ] || die "row $ROW has no target region for sample $SAMPLE"

# Comma-separated targets become repeated --target arguments.
TARGET_ARGS=()
IFS=',' read -r -a TLIST <<< "$TARGETS"
for t in "${TLIST[@]}"; do
  [ -n "$t" ] && TARGET_ARGS+=(--target "$t")
done

echo "[array task ${SLURM_ARRAY_TASK_ID}] sample=$SAMPLE targets=$TARGETS"

exec "$VERIFY" \
  --sample  "$SAMPLE" \
  --r1      "$R1" \
  --r2      "$R2" \
  --ref     "$REF" \
  --outdir  "$OUTDIR" \
  --threads "${SLURM_CPUS_PER_TASK:-4}" \
  "${TARGET_ARGS[@]}"
