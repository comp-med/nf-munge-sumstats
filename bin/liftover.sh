#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") --input VCF --from BUILD --to BUILD \\
           --src-ref FASTA --tgt-ref FASTA --chain CHAIN [options]

Required:
  --input     Input VCF/VCF.GZ (formatted sumstats)
  --from      Source build label, used in output naming (e.g. grch37)
  --to        Target build label, used in output naming (e.g. grch38)
  --src-ref   Reference FASTA for the source build
  --tgt-ref   Reference FASTA for the target build
  --chain     UCSC chain file for source -> target

Optional:
  --log       Log file                       (default: liftover_log.txt)
  -h, --help  Show this help

Outputs (in the working directory):
  formatted_sumstats_<from>.vcf.gz   sorted source-build VCF (P stored as String)
  formatted_sumstats_<to>.vcf.gz     lifted VCF, INFO stripped
  reject.vcf.bgz                     records that failed liftover
  <log>                              combined bcftools log
EOF
}

die() { echo "ERROR: $*" >&2; exit 1; }

# ARGUMENTS ----

INPUT="" FROM="" TO="" SRC_REF="" TGT_REF="" CHAIN="" LOG=liftover_log.txt

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)    INPUT=${2:?}; shift 2 ;;
        --from)     FROM=${2:?};  shift 2 ;;
        --to)       TO=${2:?};    shift 2 ;;
        --src-ref)  SRC_REF=${2:?}; shift 2 ;;
        --tgt-ref)  TGT_REF=${2:?}; shift 2 ;;
        --chain)    CHAIN=${2:?}; shift 2 ;;
        --log)      LOG=${2:?};   shift 2 ;;
        -h|--help)  usage; exit 0 ;;
        *)          usage >&2; die "Unknown argument: $1" ;;
    esac
done

for var in INPUT FROM TO SRC_REF TGT_REF CHAIN; do
    [[ -n "${!var}" ]] || { usage >&2; die "Missing required argument for ${var}"; }
done
for f in "$INPUT" "$SRC_REF" "$TGT_REF" "$CHAIN"; do
    [[ -r "$f" ]] || die "File not readable: $f"
done

SRC_VCF="formatted_sumstats_${FROM}.vcf.gz"
OUT_VCF="formatted_sumstats_${TO}.vcf.gz"

# Guard against clobbering a staged input with the same name as an output
[[ "$(basename "$INPUT")" != "$SRC_VCF" ]] || die "Input name collides with output name: $SRC_VCF"

# LOGGING ----
# Everything below (stdout + stderr, including every side of each pipe) goes to
# the log. fd 3 keeps the real stderr so failures still reach .command.err.

exec 3>&2
exec >> "$LOG" 2>&1
trap 'rc=$?; echo "ERROR: liftover failed (exit $rc), see $LOG" >&3; tail -n 20 "$LOG" >&3' ERR

echo "Liftover direction: ${FROM} -> ${TO}"
bcftools --version | head -n 1

# MAIN PIPELINE ----
# One streamed pass, no intermediate files (uncompressed BCF between stages):
#
#   reheader -> sort -> [tee: source VCF] -> norm -m+ -> liftover
#            -> sort -> norm -m- -> annotate -x INFO -> final VCF
#
# P is stored as String: tiny p-values underflow to 0 as 32-bit floats.

bcftools reheader \
    -h <(bcftools view -h "$INPUT" | sed -E '/^##FORMAT=<ID=P,/ s/Type=Float/Type=String/') \
    "$INPUT" \
  | bcftools sort -Oz \
  | tee "$SRC_VCF" \
  | bcftools norm --no-version -m+ -c wx -f "$SRC_REF" -Ou \
  | bcftools +liftover -Ou -- \
        -s "$SRC_REF" -f "$TGT_REF" -c "$CHAIN" \
        --reject reject.vcf.bgz --reject-type z \
  | bcftools sort -Ou \
  | bcftools norm --no-version -m- -Ou \
  | bcftools annotate -x INFO -Oz -o "$OUT_VCF"

# Strand/REF stats on source and lifted
bcftools +fixref "$SRC_VCF" -- -f "$SRC_REF"
bcftools +fixref "$OUT_VCF" -- -f "$TGT_REF"

# Count REF mismatches in the lifted file
total=$(bcftools view -H "$OUT_VCF" | wc -l)
ok=$(bcftools norm -N -c x -f "$TGT_REF" -Ou "$OUT_VCF" | bcftools view -H | wc -l)
mismatch=$((total - ok))
echo "REF mismatches: ${mismatch} / ${total}"

# Counts for the report
bcftools view -H "$SRC_VCF" | wc -l        # input records
bcftools view -H reject.vcf.bgz | wc -l    # rejected by liftover
bcftools stats "$OUT_VCF"

# Fail the task if the mismatch fraction is too high (0.1% here)
if awk -v m="$mismatch" -v t="$total" -v max=0.001 'BEGIN{exit !(t>0 && m/t>max)}'; then
    echo "ERROR: REF mismatch rate above threshold" >&2
    exit 1
fi

echo "Done: ${OUT_VCF}"