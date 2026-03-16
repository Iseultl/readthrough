process COMBINE_ORFSECIS {
    publishDir "${params.output_dir}/", mode: 'copy'
    
    input:
    path orfsearch_result
    path all_secis
    path filtered_secis

    output:
    path "ORFsearch_SECIS.result"

    script:
    """
    # Combine ORFsearch results with filtered SECISearch results
    add_SECIS_annotation.py --ORFsearch ${orfsearch_result} --filtered_secis ${filtered_secis} --all_secis ${all_secis} --output ORFsearch_SECIS.result
    """
}