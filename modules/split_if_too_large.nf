process SPLIT_IF_TOO_LARGE {
    tag { "split_if_needed_${input_file.baseName}" }
    publishDir "${params.output_dir}/split_transcripts", mode: 'copy'
    cpus 1
    memory '2GB'
    label 'splitfasta'

    input:
    path input_file

    output:
    tuple val(input_file.baseName), path("*.fa")

    script:
    """
    # Set thresholds
    SEQ_COUNT=\$(seqkit stats ${input_file} | awk 'NR==2 {print \$4}')
    THRESHOLD=100000  # adjust this based on your 8GB RAM target

    if [ "\$SEQ_COUNT" -gt "\$THRESHOLD" ]; then
        # Split into groups of THRESHOLD sequences
        seqkit split -p \$(( (\$SEQ_COUNT + \$THRESHOLD - 1) / \$THRESHOLD )) -O . ${input_file}
        # Rename outputs for clarity
        for f in ${input_file.baseName}.split/*; do
            mv "\$f" "\$(basename \$f .fasta).fa"
        done
    else
        cp ${input_file} ${input_file.baseName}.fa
    fi
    """
}
