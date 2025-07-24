process SELECT_INTERESTING {
    tag { "select_interesting_${recoded_predictions}_${relocated_gtf}" }
    
    publishDir "interesting_predictions", mode: 'copy'
    label 'python'
    cpus 1
    memory '4GB'
    time '2h'
    
    input:
    path recoded_predictions
    path relocated_gtf
    
    output:
    path "score.csv", emit: score
    path "longest.csv", emit: longest
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    select_interesting_recodings.py \
        --geneid ${recoded_predictions} \
        --gff ${relocated_gtf} \
        --score score.csv \
        --longest longest.csv
    """
}
