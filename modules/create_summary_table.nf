process CREATE_SUMMARY_TABLE {
    tag { "create_summary_${original_scores}_${interesting_predictions}_${relocated_gtf}" }
    label 'python'
    publishDir "${params.output_dir}/summary_results", mode: 'copy'

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
    
    # Extract identifier (e.g., NW_027221872.1)
    ID=\$(basename ${original_scores_score} | grep -oE '\\b[A-Z]{2}_[0-9]+\\.[0-9]+\\b')

    # Create output folder
    mkdir -p summary_results

    # Run script with output prefix
    create_og_re_table.py \\
        --og_score ${original_scores_score} \\
        --og_length ${original_scores_longest} \\
        --re_score ${interesting_predictions_score} \\
        --re_length ${interesting_predictions_longest} \\
        --gff ${relocated_gtf} \\
        --output summary_results/\${ID}
    """
}
