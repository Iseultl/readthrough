process FILTER_FINAL_TABLE {
    tag { "filter_final_table_${input_file}" }
    label 'python'
    cpus 1
    memory '4GB'
    time '2h'
    
    input:
    path("ORFsearch.filter")
    
    output:
    path("ORFsearch.temp")
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    filter_final_table.py \\
        ORFsearch.filter \\
        ORFsearch.temp
    """
}