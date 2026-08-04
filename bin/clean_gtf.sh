#!/bin/bash

# Input parameters
INPUT_FILE="$1"
OUTPUT_FILE="${INPUT_FILE%.*}.cleaned.gtf"

# Clean the GTF file to keep only transcript_name or fall back to transcript id if not present
# Remove lines starting with #, lines with MT in the first column, and lines with gene in the third column
awk -F'\t' '
BEGIN { OFS="\t" }

/^#/ || $1=="MT" || $3=="gene" { next }

match($9,/transcript_name "[^"]+"/){
    id = substr($9,RSTART+17,RLENGTH-18)
    $9=id
    print
    next
}

match($9,/transcript_id "[^"]+"/){
    id = substr($9,RSTART+15,RLENGTH-16)
    $9=id
    print
}
' "$INPUT_FILE" > "$OUTPUT_FILE"
