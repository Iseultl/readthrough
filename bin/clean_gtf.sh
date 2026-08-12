#!/bin/bash

set -euo pipefail

INPUT_FILE="$1"
OUTPUT_FILE="${INPUT_FILE%.*}.cleaned.gtf"

awk -F'\t' '
BEGIN {
    OFS="\t"
}

# Skip comments
/^#/ {
    next
}

# Require at least 9 GFF/GTF columns
NF < 9 {
    next
}

# Skip mitochondrial chromosome
$1 == "MT" || $1 == "chrM" || $1 == "M" {
    next
}

# Skip gene features
$3 == "gene" {
    next
}

{
    gene_id = ""
    transcript_id = ""
    gene_name = ""
    transcript_name = ""

    # ---------------------------------------------------------
    # GTF attributes
    # ---------------------------------------------------------

    if (match($9, /gene_id "[^"]+"/))
        gene_id = substr($9, RSTART + 9, RLENGTH - 10)

    if (match($9, /transcript_id "[^"]+"/))
        transcript_id = substr($9, RSTART + 15, RLENGTH - 16)

    if (match($9, /gene_name "[^"]+"/))
        gene_name = substr($9, RSTART + 10, RLENGTH - 11)

    if (match($9, /transcript_name "[^"]+"/))
        transcript_name = substr($9, RSTART + 16, RLENGTH - 17)


    # ---------------------------------------------------------
    # GFF3 attributes
    # ---------------------------------------------------------

    # ID=...
    if (gene_id == "" && match($9, /(^|;)ID=[^;]+/)) {
        value = substr($9, RSTART, RLENGTH)

        sub(/^.*ID=/, "", value)

        # A gene feature has ID=gene, while a transcript
        # can also have ID=transcript.
        if ($3 == "transcript" || \
            $3 == "mRNA" || \
            $3 == "mrna") {
            transcript_id = value
        }
        else {
            gene_id = value
        }
    }

    # Parent=...
    if (match($9, /(^|;)Parent=[^;]+/)) {
        value = substr($9, RSTART, RLENGTH)
        sub(/^.*Parent=/, "", value)

        # Parent of a transcript is generally the gene.
        if ($3 == "transcript" || \
            $3 == "mRNA" || \
            $3 == "mrna") {
            gene_id = value
        }
        # TOGA exon/CDS rows only carry the transcript parent.
        # Remember the gene for later rows that share the same transcript.
        else if (transcript_id != "") {
            if (gene_by_transcript[transcript_id] == "") {
                gene_by_transcript[transcript_id] = value
            }
        }
    }


    # ---------------------------------------------------------
    # For exon/CDS records, Parent is usually the transcript
    # ---------------------------------------------------------

    if (transcript_id == "" && \
        ($3 == "exon" || $3 == "CDS" || \
         $3 == "start_codon" || $3 == "stop_codon")) {

        if (match($9, /(^|;)Parent=[^;]+/)) {
            value = substr($9, RSTART, RLENGTH)
            sub(/^.*Parent=/, "", value)
            transcript_id = value
        }
    }

    if ((($3 == "transcript") || ($3 == "mRNA") || ($3 == "mrna")) && \
        gene_id != "" && transcript_id != "" && \
        gene_by_transcript[transcript_id] == "") {
        gene_by_transcript[transcript_id] = gene_id
    }

    if (gene_id == "" && transcript_id != "" && \
        gene_by_transcript[transcript_id] != "") {
        gene_id = gene_by_transcript[transcript_id]
    }

    if (gene_id == "" && transcript_id != "") {
        gene_id = transcript_id
    }


    # ---------------------------------------------------------
    # Require gene_id and transcript_id
    # ---------------------------------------------------------

    if (gene_id == "" || transcript_id == "")
        next


    # ---------------------------------------------------------
    # Construct clean GTF attributes
    # ---------------------------------------------------------

    attrs = "gene_id \"" gene_id "\"; "
    attrs = attrs "transcript_id \"" transcript_id "\";"

    if (gene_name != "")
        attrs = attrs " gene_name \"" gene_name "\";"

    if (transcript_name != "")
        attrs = attrs " transcript_name \"" transcript_name "\";"

    $9 = attrs

    print
}
' "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Created: $OUTPUT_FILE"