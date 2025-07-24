process CREATE_SUMMARY_TABLE {
    tag { "create_summary_${original_scores}_${interesting_predictions}_${relocated_gtf}_${secis_gff}" }
    label 'python'
    publishDir "summary_results", mode: 'copy'
    cpus 1
    memory '4GB'
    time '2h'
    
    input:
    path(interesting_predictions_score)
    path(interesting_predictions_longest)
    path(original_scores_score)
    path(original_scores_longest)
    path(relocated_gtf)
    
    output:
    path("summary_results/*")
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    create_og_re_table.py \
        --og_score ${original_scores_score} \
        --og_length ${original_scores_longest} \
        --re_score ${interesting_predictions_score} \
        --re_length ${interesting_predictions_longest} \
        --gff ${relocated_gtf} \
        --output summary_results
    """
}
