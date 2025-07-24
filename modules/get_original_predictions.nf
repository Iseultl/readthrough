process GET_ORIGINAL_PREDICTIONS {
    tag { "get_original_${geneid_result}" }
    label 'python'
    
    publishDir "original_predictions", mode: 'copy'
    cpus 1
    memory '2GB'
    time '1h'
    
    input:
    path geneid_output

    output:
    path "original_score.csv", emit: score
    path "original_longest.csv", emit: longest
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    get_og_predictions.py \
        --geneid ${input_dir}/combined_predictions_original.txt \
        --score original_score.csv \
        --longest original_longest.csv
    """
}
