process GET_ORIGINAL_PREDICTIONS {
    tag { "get_original_${geneid_output.simpleName}" }
    label 'python'

    publishDir "${params.output_dir}/original_predictions", mode: 'copy', pattern: "*_*.csv"

    cpus 1
    memory '8GB'
    time '1h'

    input:
    path geneid_output

    output:
    path "*_score.csv", emit: score
    path "*_longest.csv", emit: longest

    script:
    """
    set -euo pipefail

    base_name=\$(basename ${geneid_output} .gff)

    get_og_predictions.py \\
        --geneid ${geneid_output} \\
        --score "\${base_name}_original_score.csv" \\
        --longest "\${base_name}_original_longest.csv"
    """
}
