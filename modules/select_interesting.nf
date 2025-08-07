process SELECT_INTERESTING {
    tag { "select_interesting_${geneid_output.simpleName}" }
    
    publishDir "${params.output_dir}/interesting_predictions", mode: 'copy', pattern: "*_*.csv"
    label 'python'
    cpus 1
    memory '8GB'
    time '1h'
    
    input:
    path geneid_output
    path relocated_gtf
    
    output:
    tuple val(id), path("*_recoded_score.csv"), path("*_recoded_longest.csv"), emit: interesting_predictions
    
    script:
    def id = geneid_output.getName().replaceAll('transcripts_clean_', '').replaceAll(/_recoded.*/, '')
    """
    #!/bin/bash
    set -euo pipefail
    
    select_interesting_recodings.py \\
        --geneid ${geneid_output} \\
        --gff ${relocated_gtf} \\
        --score "${id}_score.csv" \\
        --longest "${id}_longest.csv"
    """
}
