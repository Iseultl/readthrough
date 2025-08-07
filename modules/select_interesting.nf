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
    tuple val(id), path("*_score.csv"), path("*_longest.csv"), emit: interesting_predictions
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    id=\$(echo "${geneid_output}" | grep -oE '(NC|NW|NT|AC|NG|NM|NR|NP|XM|XR|XP|YP|ZP)_[0-9]+\\.[0-9]+\\.part_[0-9]+')
    base_name=\$(basename ${geneid_output} .gff)

    select_interesting_recodings.py \\
        --geneid ${geneid_output} \\
        --gff ${relocated_gtf} \\
        --score "${id}_score.csv" \\
        --longest "${id}_longest.csv"
    """
}
