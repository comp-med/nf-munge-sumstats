#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") --input VCF --from BUILD --to BUILD \\
           --src-ref FASTA --tgt-ref FASTA --chain CHAIN [options]

Required:
  --input     Input VCF/VCF.GZ (formatted sumstats)
  --build     Source build label, used in output naming (e.g. grch37)
  --src-ref   Reference FASTA for the source build

Optional:
  --log       Log file                       (default: liftover_log.txt)
  -h, --help  Show this help

Outputs (in the working directory):
  formatted_sumstats_<from>.vcf.gz   sorted source-build VCF (P stored as String)
  <log>                              combined bcftools log
EOF
}

die() { echo "ERROR: $*" >&2; exit 1; }

# ARGUMENTS ----

INPUT="" FROM="" SRC_REF="" LOG=no_liftover_log.txt

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)    INPUT=${2:?}; shift 2 ;;
        --from)     FROM=${2:?};  shift 2 ;;
        --src-ref)  SRC_REF=${2:?}; shift 2 ;;
        --log)      LOG=${2:?};   shift 2 ;;
        -h|--help)  usage; exit 0 ;;
        *)          usage >&2; die "Unknown argument: $1" ;;
    esac
done

for var in INPUT FROM SRC_REF; do
    [[ -n "${!var}" ]] || { usage >&2; die "Missing required argument for ${var}"; }
done
for f in "$INPUT" "$SRC_REF"; do
    [[ -r "$f" ]] || die "File not readable: $f"
done

SRC_VCF="formatted_sumstats_${FROM}.vcf.gz"

# Guard against clobbering a staged input with the same name as an output
[[ "$(basename "$INPUT")" != "$SRC_VCF" ]] || die "Input name collides with output name: $SRC_VCF"

# LOGGING ----
# Everything below (stdout + stderr, including every side of each pipe) goes to
# the log. fd 3 keeps the real stderr so failures still reach .command.err.

exec 3>&2
exec >> "$LOG" 2>&1
trap 'rc=$?; echo "ERROR: script failed (exit $rc), see $LOG" >&3; tail -n 20 "$LOG" >&3' ERR

bcftools --version | head -n 1

bcftools reheader \
    -h <(bcftools view -h "$INPUT" | sed -E '/^##FORMAT=<ID=P,/ s/Type=Float/Type=String/') \
    "$INPUT" \
  | bcftools sort -Ou \
  | bcftools norm --no-version -m- -Ou \
  | bcftools annotate -x INFO -Oz -o "$SRC_VCF"

# Strand/REF stats on source and lifted
bcftools +fixref "$SRC_VCF" -- -f "$SRC_REF"

# Count REF mismatches in the lifted file
total=$(bcftools view -H "$SRC_VCF" | wc -l)
ok=$(bcftools norm -N -c x -f "$SRC_REF" -Ou "$SRC_VCF" | bcftools view -H | wc -l)
mismatch=$((total - ok))
echo "REF mismatches: ${mismatch} / ${total}"

# Counts for the report
bcftools view -H "$SRC_VCF" | wc -l        # input records
bcftools stats "$SRC_VCF"

# Fail the task if the mismatch fraction is too high (0.1% here)
if awk -v m="$mismatch" -v t="$total" -v max=0.001 'BEGIN{exit !(t>0 && m/t>max)}'; then
    echo "ERROR: REF mismatch rate above threshold" >&2
    exit 1
fi

echo "Done: ${SRC_VCF}"