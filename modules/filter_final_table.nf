process FILTER_FINAL_TABLE {
    tag { "filter_final_table" }
    label 'python'
    cpus 1
    time '2h'
    memory '4GB'
    
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