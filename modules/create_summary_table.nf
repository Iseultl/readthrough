process CREATE_SUMMARY_TABLE {
    tag { "create_summary_${id}" }
    label 'python'
    cpus 1
  
    input:
    tuple val(id), path(interesting_score), path(interesting_longest), path(original_score), path(original_longest)
    path relocated_gtf
    
    output:
    path("summary_results/${id}.csv")
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "Processing ID: ${id}"

    # Create output folder
    mkdir -p summary_results

    # Run script with output prefix
    create_og_re_table.py \
        --og_score ${original_score} \
        --og_length ${original_longest} \
        --re_score ${interesting_score} \
        --re_length ${interesting_longest} \
        --gff ${relocated_gtf} \
        --output summary_results/${id}
    """
}
