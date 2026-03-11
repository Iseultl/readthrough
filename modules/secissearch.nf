// modules/secissearch.nf

process SECISSEARCH {
    tag { input_fasta.getBaseName() }

    publishDir "${params.output_dir}/secis_final_gff", mode: 'copy'

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

    python /Seblastian/Seblastian.py \\
        -t ${input_fasta} \\
        -o ${fileName} \\
        -SS \\
        -temp secis_temp \\
        -no_complement \\
        -secis_energy 0.0002
    """
}