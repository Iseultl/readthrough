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
    tuple val(id), path("*_original_score.csv"), path("*_original_longest.csv"), emit: original_predictions

    def id = geneid_output.getName().replaceAll('transcripts_clean_', '').replaceAll(/_recoded.*/, '')
    
    script:
    """
    get_og_predictions.py \\
        --geneid ${geneid_output} \\
        --score "${id}_original_score.csv" \\
        --longest "${id}_original_longest.csv"
    """
}
