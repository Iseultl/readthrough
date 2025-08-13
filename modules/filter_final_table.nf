process FILTER_FINAL_TABLE {
    tag { "filter_final_table_${input_file}" }
    label 'python'
    publishDir "${params.output_dir}", mode: 'copy'

    cpus 1
    memory '4GB'
    time '2h'
    
    input:
    path("ORFsearch.filter")
    
    output:
    path("ORFsearch.result")
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    filter_final_table.py \\
        --input ORFsearch.filter \\
        --output ORFsearch.result
    """
}