// modules/secissearch.nf

process SECISSEARCH {
    tag { input_fasta.baseName }
    memory '4GB' 
    label 'secissearch'

    input:
    path(input_fasta)

    output:
    path("*.gff")

    script:
    def fileName = input_fasta.baseName
    """
    #!/bin/bash
    set -euo pipefail

    mkdir -p secis_temp

    # Remove sequences who's length is greater than 4000
    filter_fasta_length.py \\
    --input ${input_fasta} \\
    --output ${input_fasta}.filtered.fa

    python /Seblastian/Seblastian.py \\
        -t ${input_fasta} \\
        -o ${fileName} \\
        -SS \\
        -temp secis_temp \\
        -no_complement \\
        -secis_energy 0.0002

    rm -r secis_temp ${input_fasta}
    """
}