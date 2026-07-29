process COMBINE_ORFSECIS {
    publishDir "${params.output_dir}/", mode: 'copy'
    label 'python'
    memory '4GB'

    input:
    path orfsearch_result
    path all_secis
    path filtered_secis
    val params.species_name

    output:
    path "ORFsearch_SECIS.result"

    script:
    """
    #!/bin/bash
    set -euo pipefail
    # Combine ORFsearch results with filtered SECISearch results
    add_SECIS_annotation.py --ORFsearch ${orfsearch_result} --filtered_secis ${filtered_secis} --all_secis ${all_secis} --output ${params.species_name}_ORFsearch_SECIS.result
    """
}