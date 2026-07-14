#!/usr/bin/env bash
# Usage: ./filter.sh input.tsv gene_list.txt output.tsv

TSV_FILE="$1"
GENE_LIST="$2"
OUTPUT_FILE="$3"

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 input.tsv gene_list.txt output.tsv"
    exit 1
fi

# Detect column positions
GENE_COL=$(head -n 1 "$TSV_FILE" | tr '\t' '\n' | grep -n "^GENE_SYMBOL$" | cut -d: -f1)
MUT_COL=$(head -n 1 "$TSV_FILE" | tr '\t' '\n' | grep -n "^MUTATION_DESCRIPTION$" | cut -d: -f1)

if [[ -z "$GENE_COL" || -z "$MUT_COL" ]]; then
    echo "Could not find GENE_SYMBOL or MUTATION_DESCRIPTION column."
    exit 1
fi

# Write header
head -n 1 "$TSV_FILE" > "$OUTPUT_FILE"

# Filter
awk -v gcol="$GENE_COL" 'BEGIN{FS=OFS="\t"}
    NR==FNR {
        g=$1; sub(/\r$/,"",g);  # strip possible CR from Windows line endings
        if (g != "") genes[g]=1;
        next
    }
    NR>1 && ($gcol in genes)
' "$GENE_LIST" "$TSV_FILE" >> "$OUTPUT_FILE"