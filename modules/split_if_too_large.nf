process SPLIT_IF_TOO_LARGE {
    tag { "split_if_needed_${input_file.baseName}" }
    publishDir "${params.output_dir}/split_transcripts", mode: 'copy'
    cpus 1
    memory '2GB'
    label 'splitfasta'

    input:
    path input_file

    output:
    path("*.fa")

    script:
    """
    echo "Processing file: ${input_file}"
    # Use seqkit stats with tabular output and headers to robustly get sequence count
    SEQ_COUNT=\$(seqkit stats -T ${input_file} | awk -F '\t' 'NR==2 {print \$4}')
    echo "Sequence count: \$SEQ_COUNT"
    THRESHOLD=80000

    if [ -n "\$SEQ_COUNT" ] && [ "\$SEQ_COUNT" -gt "\$THRESHOLD" ]; then
        echo "Splitting file as sequence count \$SEQ_COUNT is greater than threshold \$THRESHOLD"
        # Calculate number of parts to split into
        PARTS=\$(( (\$SEQ_COUNT + \$THRESHOLD - 1) / \$THRESHOLD ))
        seqkit split -p \$PARTS -O . ${input_file}

        # Rename outputs for clarity and ensure they are in the correct output directory
        for f in ${input_file.baseName}.split/*; do
            # The output files from seqkit split will be named like 'input.part_001.fasta'
            # We want to move them to the work directory with a simpler name
            mv "\$f" "\$(basename \${f%.fasta}).fa"
        done
    else
        echo "No split needed"
        cp ${input_file} ${input_file.baseName}.fa
    fi
    """
}

