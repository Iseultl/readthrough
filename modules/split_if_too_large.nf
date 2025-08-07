process SPLIT_IF_TOO_LARGE {
    tag { "split_if_needed_${input_file.baseName}" }
    cpus 1
    memory '8GB'
    label 'splitfasta'

    input:
    path input_file

    output:
    path "*.fa", emit: split_fasta

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
    else
        echo "No split needed"
        mv ${input_file} ${input_file.baseName}.part_001.fa
    fi
    """
}

