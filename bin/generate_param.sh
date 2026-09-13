#!/usr/bin/env bash

set -euo pipefail

start_time=$(date +%s)

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <species_name> <genome_fasta> <cds_gtf> <output_dir>"
    exit 1
fi

SPECIES_NAME="$1"
FASTA="$2"
GFF="$3"
OUTPUT_DIR="$4"


#############################################
#### Define output names ####################
#############################################

CLEANED_FASTA="${FASTA%.fasta}.cleaned.fasta"

GFF_BASE="${GFF%.gff}"
GFF_BASE="${GFF_BASE%.gtf}"

CDS_GFF="${GFF_BASE}.CDS.gff"
SAMPLED_GFF="${GFF_BASE}.sampled.CDS.gff"


#############################################
#### Clean FASTA ############################
#############################################

if [[ "$FASTA" == *.gz ]]; then

    gunzip -c "$FASTA" |
        awk '/^>/ {print $1; next} {print}' \
        > "$CLEANED_FASTA"

else

    awk '/^>/ {print $1; next} {print}' \
        "$FASTA" > "$CLEANED_FASTA"

fi


#############################################
#### Clean GFF ##############################
#############################################

awk -F'\t' 'BEGIN {OFS="\t"} $3 == "CDS" {
    if (match($9, /transcript_id "[^"]+"/)) {
        transcript = substr($9, RSTART, RLENGTH)
        sub(/^transcript_id "/, "", transcript)
        sub(/"$/, "", transcript)
        $9 = transcript
        print
    }
}' "$GFF" > "$CDS_GFF"


#############################################
#### Sample the gff  ########################
#############################################

sample_size=2000

# Extract unique transcript IDs
awk -F'\t' '{print $9}' "$CDS_GFF" \
    | sort -u \
    | shuf -n "$sample_size" \
    > sampled_ids.txt

# Keep all CDS lines belonging to sampled transcript IDs
awk -F'\t' '
NR == FNR {
    ids[$0] = 1
    next
}
$9 in ids
' sampled_ids.txt "$CDS_GFF" \
    > "$SAMPLED_GFF"

echo "Sampled GFF file created: $SAMPLED_GFF"

#############################################
#### run the trainer ########################
#############################################

TARGET_OUT="${OUTPUT_DIR}/${SPECIES_NAME}"
mkdir -p "$TARGET_OUT"

/scripts_geneid/geneidTRAINer4docker.pl \
    -species "$SPECIES_NAME" \
    -gff "$SAMPLED_GFF" \
    -fastas "$CLEANED_FASTA" \
    -results "$TARGET_OUT" \
    -reduced no

PARAM_FILE=$(find "$TARGET_OUT" -type f -name '*.param' -print -quit)
if [[ -z "$PARAM_FILE" ]]; then
    echo "No parameter file was produced in $TARGET_OUT" >&2
    exit 1
fi
cp "$PARAM_FILE" "${OUTPUT_DIR}/${SPECIES_NAME}.param"

# Record end time
end_time=$(date +%s)
runtime=$((end_time - start_time))


echo "Script runtime: $runtime seconds"
echo "Memory usage saved in slurm_memory_usage.txt"

#############################################
#### Make Param Single Exon Mode ############
#############################################

PARAM_FILE="${OUTPUT_DIR}/${SPECIES_NAME}.param"

# Remove species lines: the trainer may insert lines with the species name or ID
# that are not valid geneid directives
sed -i "/${SPECIES_NAME}/d" "$PARAM_FILE"

# Remove any standalone "species" keyword lines (case-insensitive) not part of a known directive
sed -i '/^[Ss][Pp][Ee][Cc][Ii][Ee][Ss]$/d' "$PARAM_FILE"

# Replace the gene model section with single-exon-only rules.
# The trainer may produce a multi-exon model with acceptor/donor connections,
# but our pipeline works on transcript-relative coordinates where each "gene"
# is a single CDS, so we only need Begin -> Single -> End connections.
sed -i '/^# GENE MODEL:/,$d' "$PARAM_FILE"
cat >> "$PARAM_FILE" << 'EOF'
# GENE MODEL: Rules about gene assembling (GenAmic)
General_Gene_Model
# BEGINNING and END of prediction
Begin+                          Single+ 0:Infinity
Begin-                          Single- 0:Infinity
Single+               End+              0:Infinity
Single-               End-              0:Infinity
EOF