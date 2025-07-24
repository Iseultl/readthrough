#!/bin/bash

# Input parameters from Nextflow
SPLIT_DIR="$1"    # e.g. split_transcripts/split_transcripts
OUTPUT_DIR="$2"   # geneid_original_predictions
PARAM_FILE="$3"   # geneid parameter file from params.yaml

mkdir -p "$OUTPUT_DIR"

for fasta in $SPLIT_DIR/*.fa; do
    echo "Processing $(basename "$fasta")"
    output_file="$OUTPUT_DIR/$(basename "$fasta").geneid.txt"
    geneid -P "$PARAM_FILE" -s -W "$fasta" > "$output_file"
done



