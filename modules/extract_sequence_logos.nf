process EXTRACT_SEQUENCE_LOGOS {
    tag { "extract_sequence_logos_${result_csv}" }
    label 'python'
    publishDir "${params.output_dir}/${params.species_name}", mode: 'copy'
    cpus 1
    time '1h'
    memory '4GB'
    
    input:
    path result_csv
    path gffread_dir
    val species
     
    output:
    path "${species}_sequence_logos/*"
    
    script:
    """
    #!/bin/bash
    set -euo pipefail

    mkdir -p ${species}_sequence_logos
    
    extract_sequence_logos.py \
        --result_file ${result_csv} \
        --fasta_dir ${gffread_dir} \
        --output_dir ${species}_sequence_logos

    rm -rf ${gffread_dir} ${result_csv}
    """
}