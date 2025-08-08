process GET_ORIGINAL_PREDICTIONS {
    tag { "get_original_${geneid_output.simpleName}" }
    label 'python'

    publishDir "${params.output_dir}/original_predictions", mode: 'copy', pattern: "*_original_*.csv"

    cpus 1
    memory '8GB'
    time '1h'

    input:
    path geneid_output

    output:
    tuple val(id), path("${id}_original_score.csv"), path("${id}_original_longest.csv"), emit: original_predictions

    script:
    id = geneid_output.getName().replaceAll('transcripts_clean_', '').replaceAll('.fa.geneid.txt_*', '')
    """
    get_og_predictions.py \
        --geneid ${geneid_output} \
        --score "${id}_original_score.csv" \
        --longest "${id}_original_longest.csv"
    """
}
