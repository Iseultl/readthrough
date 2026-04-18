process SELECT_INTERESTING {
    tag { "select_interesting_${geneid_output.simpleName}" }
    label 'python'
    cpus 1
    memory '16GB'
    time '2h'

    input:
    path geneid_output
    path relocated_gtf

    output:
    tuple val(id), path("${id}_recoded_score.csv"), path("${id}_recoded_longest.csv"), emit: interesting_predictions

    script:
    id = geneid_output.getName().replaceAll("transcripts_clean_", "").replaceAll("_recoded.", ".").replaceAll(".fa.geneid.txt_recoded_*", "")
    """
    select_interesting_recodings.py \
        --geneid ${geneid_output} \
        --gff ${relocated_gtf} \
        --score "${id}_recoded_score.csv" \
        --longest "${id}_recoded_longest.csv"
    """
}
