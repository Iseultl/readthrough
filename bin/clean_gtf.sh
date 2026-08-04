#!/bin/bash

# Input parameters
INPUT_FILE="$1"
OUTPUT_FILE="${INPUT_FILE%.*}.cleaned.gtf"

# Clean the GTF file to keep only transcript_name or fall back to transcript id if not present
# Remove lines starting with #, lines with MT in the first column, and lines with gene in the third column
awk -F'\t' '
BEGIN {
    OFS="\t"
}

# Skip comments
/^#/ { next }

# Skip mitochondrial chromosome
$1 == "MT" { next }

# Skip gene features
$3 == "gene" { next }

{
    gene_id = ""
    transcript_id = ""
    gene_name = ""
    transcript_name = ""

    if (match($9, /gene_id "[^"]+"/))
        gene_id = substr($9, RSTART, RLENGTH)

    if (match($9, /transcript_id "[^"]+"/))
        transcript_id = substr($9, RSTART, RLENGTH)

    if (match($9, /gene_name "[^"]+"/))
        gene_name = substr($9, RSTART, RLENGTH)

    if (match($9, /transcript_name "[^"]+"/))
        transcript_name = substr($9, RSTART, RLENGTH)

    # Require the two essential attributes
    if (gene_id == "" || transcript_id == "")
        next

    attrs = gene_id "; " transcript_id

    if (gene_name != "")
        attrs = attrs "; " gene_name

    if (transcript_name != "")
        attrs = attrs "; " transcript_name

    attrs = attrs ";"

    $9 = attrs
    print
}
' "$INPUT_FILE" > "$OUTPUT_FILE"
