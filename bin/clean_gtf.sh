#!/bin/bash

# Input parameters
INPUT_FILE="$1"
OUTPUT_FILE="${INPUT_FILE%.*}.cleaned.gtf"

# Clean the GTF file to keep only transcript_name
grep -v $'\tgene\t' ${INPUT_FILE} | awk -F'\t' '{
    split($9, fields, ";")
    for (i in fields) {
        if (fields[i] ~ /transcript_name/) {
            split(fields[i], kv, "\"")
            $9 = kv[2]
            break
        }
    }
    print $0
}' OFS='\t' > "$OUTPUT_FILE"
