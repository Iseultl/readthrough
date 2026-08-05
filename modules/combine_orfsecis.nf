process COMBINE_ORFSECIS {
    publishDir "${params.output_dir}/${species_name}", mode: 'copy'
    memory '4GB'

    input:
    path orfsearch_result
    path all_secis
    path filtered_secis
    val species_name

    output:
    path "*_ORFsearch_SECIS.result"

    script:
    """
    #!/bin/bash
    set -euo pipefail
    # Combine ORFsearch results with filtered SECISearch results
    add_SECIS_annotation.py --ORFsearch ${orfsearch_result} --filtered_secis ${filtered_secis} --all_secis ${all_secis} --output ${species_name}_ORFsearch_SECIS.result
    """
}