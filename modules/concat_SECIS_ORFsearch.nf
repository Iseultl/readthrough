process concat_SECIS_ORFsearch {
    tag "Concatenating SECIS and ORFsearch results"
    publishDir "${params.output_dir}", mode: 'copy'
    label 'python'
    
    input:
    path merged_secis_gff
    path filtered_secis
    path orfsearch_result

    output:
    path "ORFsearch.result"

    script:
    """
    # Python script to concatenate the results 
    add_SECIS_annotation.py \\
        --orfsearch ${orfsearch_result} \\
        --secis_gff ${merged_secis_gff} \\
        --filtered_sec ${filtered_secis} \\
        --output ORFsearch.result
    """
}