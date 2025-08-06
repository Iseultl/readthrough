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
    path "*_score.csv", emit: score
    path "*_longest.csv", emit: longest
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    base_name=\$(basename ${geneid_output} .gff)

    select_interesting_recodings.py \\
        --geneid ${geneid_output} \\
        --gff ${relocated_gtf} \\
        --score "\${base_name}_score.csv" \\
        --longest "\${base_name}_longest.csv"
    """
}
