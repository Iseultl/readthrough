process GET_ORIGINAL_PREDICTIONS {
    tag { "get_original_${geneid_output.simpleName}" }
    label 'python'
    cpus 1
    time '1h'
    memory '4GB'
    
    input:
    path geneid_output

    output:
    tuple val(id), path("${id}_original_score.csv"), path("${id}_original_longest.csv"), emit: original_predictions

    script:
    id = geneid_output.getName().replaceAll("transcripts_clean_", "").replaceAll("_recoded", "").replaceAll(".fa.geneid.txt_original_*", "")
    """
    get_og_predictions.py \
        --geneid ${geneid_output} \
        --score "${id}_original_score.csv" \
        --longest "${id}_original_longest.csv"

    rm ${geneid_output}
    """
}
