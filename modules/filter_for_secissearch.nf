// modules/filter_for_secissearch.nf

process FILTER_FOR_SECISSEARCH {
    tag { input_fasta.baseName }
    memory '4GB' 
    label 'python'

    input:
    path(input_fasta)

    output:
    path("${input_fasta}.filtered.fa")

    script:
    """
    #!/bin/bash
    # Remove sequences who's length is greater than 4000
    filter_fasta_length.py \\
    --input ${input_fasta} \\
    --output ${input_fasta}.filtered.fa
    rm ${input_fasta}
    """
}