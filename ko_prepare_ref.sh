#!/bin/bash
# Build a chimeric reference (host genome + KO cassette as an extra contig) and
# index it for bwa. Run ONCE per genome+cassette pair, before any ko_verify.sh
# job. Splitting this out keeps parallel array jobs from racing to build the
# same bwa index.
#
# Usage:
#   ko_prepare_ref.sh --genome ref.fasta --cassette cassette.fasta --outdir DIR
#
# Produces DIR/chimeric.fasta (+ bwa index, .fai) and DIR/chimeric.manifest

set -e

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*" >&2; }

GENOME=""; CASSETTE=""; OUTDIR=""
while [ $# -gt 0 ]; do
  case "$1" in
    --genome)   GENOME="$2";   shift 2 ;;
    --cassette) CASSETTE="$2"; shift 2 ;;
    --outdir)   OUTDIR="$2";   shift 2 ;;
    -h|--help)  sed -n '2,12p' "$0"; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

# --- validate inputs BEFORE doing any work ---------------------------------
[ -n "$GENOME" ]   || die "--genome is required"
[ -n "$CASSETTE" ] || die "--cassette is required"
[ -n "$OUTDIR" ]   || die "--outdir is required"
[ -s "$GENOME" ]   || die "genome FASTA missing or empty: $GENOME"
[ -s "$CASSETTE" ] || die "cassette FASTA missing or empty: $CASSETTE"

# Module names are overridable: the TSL wrapper name does not reliably match
# the binary it delivers (sourcing bwa-0.7.7 has produced bwa 0.7.17), so the
# name is a hint and the version assert below is the real check.
BWA_MODULE="${BWA_MODULE:-bwa-0.7.17}"
SAMTOOLS_MODULE="${SAMTOOLS_MODULE:-samtools-1.10}"

# Activation scripts are not written for strict mode. The TSL wrapper reads
# $WRAPPER_DEBUG unset (dies under -u) and contains a harmlessly-failing
# command (dies under -e, exit 13). Turn both off across the source, back on
# after. Note "source X || true" does NOT contain this: errexit still fires
# inside the sourced file.
set +eu
source "$BWA_MODULE"
source "$SAMTOOLS_MODULE"
set -e

# Check activation by its outcome, not its exit status. sbatch defaults to
# --export=ALL, so PATH can be inherited from the submitting shell and look
# like success -- report what we actually got so that is visible in the log.
command -v bwa      >/dev/null || die "bwa not on PATH after sourcing $BWA_MODULE"
command -v samtools >/dev/null || die "samtools not on PATH after sourcing $SAMTOOLS_MODULE"
log "bwa:      $(command -v bwa)      $(bwa 2>&1 | awk '/^Version/ { print $2; exit }')"
log "samtools: $(command -v samtools) $(samtools --version | awk 'NR == 1 { print $2 }')"

mkdir -p "$OUTDIR"
REF="$OUTDIR/chimeric.fasta"
MANIFEST="$OUTDIR/chimeric.manifest"

# --- rename cassette records to a reserved namespace ------------------------
# A fixed name (KO_CASSETTE / KO_CASSETTE_n) means every downstream script can
# find the cassette contig without being told what its header said.
NCAS=$(grep -c '^>' "$CASSETTE" || true)
[ "$NCAS" -ge 1 ] || die "no FASTA records found in cassette: $CASSETTE"

if grep -q '^>KO_CASSETTE' "$GENOME"; then
  die "genome already contains a contig named KO_CASSETTE*; rename it first"
fi

CAS_RENAMED="$OUTDIR/cassette.renamed.fasta"
awk -v n="$NCAS" '
  /^>/ { i++; if (n == 1) print ">KO_CASSETTE"; else print ">KO_CASSETTE_" i; next }
  { print }
' "$CASSETTE" > "$CAS_RENAMED"

CAS_CONTIGS=$(grep '^>' "$CAS_RENAMED" | tr -d '>' | paste -sd, -)
log "cassette contigs: $CAS_CONTIGS"

# --- build and index --------------------------------------------------------
cat "$GENOME" "$CAS_RENAMED" > "$REF"
[ -s "$REF" ] || die "chimeric reference is empty: $REF"

log "indexing with samtools faidx"
samtools faidx "$REF"

log "indexing with bwa (this takes a few minutes on a fungal genome)"
bwa index "$REF"

# --- assert on the output, not on exit status -------------------------------
[ -s "${REF}.bwt" ] || die "bwa index did not produce ${REF}.bwt"
[ -s "${REF}.fai" ] || die "samtools faidx did not produce ${REF}.fai"

for c in $(echo "$CAS_CONTIGS" | tr ',' ' '); do
  awk -v c="$c" '$1 == c { found = 1 } END { exit !found }' "${REF}.fai" \
    || die "cassette contig $c is absent from ${REF}.fai"
done

NGENOME=$(grep -c '^>' "$GENOME")
{
  echo "built=$(date '+%F %T')"
  echo "genome=$GENOME"
  echo "genome_contigs=$NGENOME"
  echo "cassette=$CASSETTE"
  echo "cassette_contigs=$CAS_CONTIGS"
  echo "reference=$REF"
} > "$MANIFEST"

log "done. reference=$REF  manifest=$MANIFEST"
