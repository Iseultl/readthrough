process CONCAT_SUMMARY_RESULTS {
    publishDir "${params.output_dir}/", mode: 'copy'
    
    input:
    path summary_tables

    output:
    path "ORFsearch.filter"

    script:
    """
    cat ${summary_tables.join(' ')} > ORFsearch.filter
    """
}
