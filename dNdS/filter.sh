#!/usr/bin/env bash
# Usage: ./filter.sh input.tsv output.tsv

set -euo pipefail

TSV_FILE="${1:-}"
OUTPUT_FILE="${2:-}"

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 input.tsv output.tsv" >&2
    exit 1
fi

# Desired columns in output (order matters)
COLUMNS="SAMPLE_NAME,CHROMOSOME,GENOME_START,GENOMIC_WT_ALLELE,GENOMIC_MUT_ALLELE"

# One-pass AWK: map headers -> indices, validate, then print selected columns
awk -v FS='\t' -v OFS='\t' -v cols="$COLUMNS" -v out="$OUTPUT_FILE" '
BEGIN {
    n = split(cols, want, ",")
}
NR == 1 {
    # Normalize potential CRLF line endings in header
    gsub(/\r$/, "", $NF)
    for (i = 1; i <= NF; i++) {
        h[$i] = i
    }
    # Check all required columns exist and record their indices
    missing = ""
    for (i = 1; i <= n; i++) {
        if (!(want[i] in h)) { missing = missing want[i] " " }
        idx[i] = h[want[i]]
    }
    if (length(missing) > 0) {
        print "Missing columns: " missing > "/dev/stderr"
        exit 1
    }
    # Write header to output (truncate/create file)
    header = want[1]
    for (i = 2; i <= n; i++) header = header OFS want[i]
    print header > out
    next
}
{
    # Normalize potential CR on each line
    gsub(/\r$/, "", $NF)
    line = $(idx[1])
    for (i = 2; i <= n; i++) line = line OFS $(idx[i])
    print line >> out
}
' "$TSV_FILE"
