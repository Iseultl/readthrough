#!/bin/bash

set -euo pipefail

# The purpose of this script is to take the predicted protein file from the 
# SecORFsearch pipeline and run it through Selenoprofiles to identify the matching hits

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <predicted_protein_directory> <genome_file> <species> <output_directory>"
    exit 1
fi

# Assign input arguments to variables
PREDICTED_PROTEIN_DIR="$1"
GENOME_FILE="$2"
SPECIES="$3"
OUTPUT_DIRECTORY="$4"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIRECTORY"

GENOME_FA="$(pwd)/genome.fasta"
SELENOPROFILES_IMAGE="${SELENOPROFILES_IMAGE:-maxtico/selenoprofiles_container:v4.5.7}"
SELENOPROFILES_HOME="${SELENOPROFILES_HOME:-$(pwd)/.selenoprofiles_home}"

run_selenoprofiles() {
    if [[ "${SELENOPROFILES_CONTAINERIZED:-0}" == "1" ]]; then
        command selenoprofiles "$@"
        return
    fi

    mkdir -p "$SELENOPROFILES_HOME"

    if command -v docker >/dev/null 2>&1; then
        docker run --rm \
            -u "$(id -u):$(id -g)" \
            -e HOME=/work/home \
            -v "$PWD:$PWD" \
            -v "$SELENOPROFILES_HOME:/work/home" \
            -w "$PWD" \
            "$SELENOPROFILES_IMAGE" \
            selenoprofiles "$@"
        return
    fi

    if command -v apptainer >/dev/null 2>&1; then
        APPTAINERENV_HOME=/work/home apptainer exec --cleanenv \
            --bind "$PWD:$PWD" \
            --bind "$SELENOPROFILES_HOME:/work/home" \
            --pwd "$PWD" \
            "$SELENOPROFILES_IMAGE" \
            selenoprofiles "$@"
        return
    fi

    if command -v singularity >/dev/null 2>&1; then
        SINGULARITYENV_HOME=/work/home singularity exec --cleanenv \
            --bind "$PWD:$PWD" \
            --bind "$SELENOPROFILES_HOME:/work/home" \
            --pwd "$PWD" \
            "$SELENOPROFILES_IMAGE" \
            selenoprofiles "$@"
        return
    fi

    command selenoprofiles "$@"
}

cleanup() {
    rm -f "$GENOME_FA"
    rm -f profile_*.fasta* 2>/dev/null || true
}
trap cleanup EXIT

shopt -s nullglob

# Copy genome file to local folder and rename it to "genome.fasta" for Selenoprofiles
cp "$GENOME_FILE" "$GENOME_FA"

# Set up Selenoprofiles
run_selenoprofiles -setup

# Loop through each split FASTA file and run Selenoprofiles
for FASTA_FILE in "${PREDICTED_PROTEIN_DIR}"/*.fasta; do
    BASENAME=$(basename "$FASTA_FILE" .fasta)

    #Replace any "." in the basename with "_" to avoid issues with Selenoprofiles
    BASENAME="${BASENAME//./_}"

    # Create a profile
    run_selenoprofiles build -i "$FASTA_FILE" -o "profile_$BASENAME.fasta" -y

    # Run Selenoprofiles search
    run_selenoprofiles -o "$BASENAME" -t "$GENOME_FA" -p "profile_$BASENAME.fasta" -tblastn -s "$SPECIES"

    # Clean up stale database locks for the next iteration
    if [[ -f "$BASENAME/$SPECIES.genome/results.sqlite" ]]; then
        run_selenoprofiles database -i "$BASENAME/$SPECIES.genome/results.sqlite" -clean 2>/dev/null || true
    fi

    # Clean up & extract only the relevant output folder
    mkdir -p "$OUTPUT_DIRECTORY/$BASENAME"
    cp -r "$BASENAME"/*/output/* "$OUTPUT_DIRECTORY/$BASENAME/" || true
    rm -rf "$BASENAME"
    rm -f "profile_${BASENAME}.fasta"* 2>/dev/null || true

done