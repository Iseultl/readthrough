#!/bin/bash

# Input parameters from Nextflow
INPUT_FILE="$1"    # e.g. split_transcripts/split_transcripts
OUTPUT_DIR="$2"   # geneid_original_predictions
PARAM_FILE="$3"   # geneid parameter file from params.yaml

mkdir -p "$OUTPUT_DIR"

echo "Processing $(basename "$INPUT_FILE")"
output_file="$OUTPUT_DIR/$(basename "$INPUT_FILE").geneid.txt"
geneid -P "$PARAM_FILE" -s -W "$INPUT_FILE" > "$output_file"




