process concat_SECIS_ORFsearch {
    tag "Concatenating SECIS and ORFsearch results"
    publishDir "${params.output_dir}", mode: 'copy'
    label 'python'
    memory '4GB'
     
    input:
    path merged_secis_gff
    path filtered_secis
    path orfsearch_result
    val params.species_name

    output:
    path "ORFsearch.result"

    script:
    """
    # Python script to concatenate the results 
    add_SECIS_annotation.py \\
        --ORFsearch ${orfsearch_result} \\
        --all_secis ${merged_secis_gff} \\
        --filtered_secis ${filtered_secis} \\
        --output ${params.species_name}_ORFsearch.result
    """
}