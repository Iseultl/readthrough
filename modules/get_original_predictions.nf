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

    id=\$(echo "${geneid_output}" | sed -e 's/transcripts_clean_//' -e 's/_recoded.*//')
    echo "ID: \${id}"
    base_name=\$(basename "${geneid_output}" .gff)

    get_og_predictions.py \\
        --geneid ${geneid_output} \\
        --score "\${id}_original_score.csv" \\
        --longest "\${id}_original_longest.csv"
    """
}
