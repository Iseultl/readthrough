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
    tuple val(id), path("*_score.csv"), path("*_longest.csv"), emit: original_predictions

    script:
    """
    set -euo pipefail

    id=\$(echo "${geneid_output}" | grep -oE '(NC|NW|NT|AC|NG|NM|NR|NP|XM|XR|XP|YP|ZP)_[0-9]+\\.[0-9]+\\.part_[0-9]+')
    base_name=\$(basename ${geneid_output} .gff)

    get_og_predictions.py \\
        --geneid ${geneid_output} \\
        --score "${id}_original_score.csv" \\
        --longest "${id}_original_longest.csv"
    """
}
