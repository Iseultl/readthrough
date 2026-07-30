process CONCAT_SUMMARY_RESULTS {
    memory '4GB'
     
    input:
    path summary_tables

    output:
    path "ORFsearch.filter"

    script:
    """
    cat ${summary_tables.join(' ')} > ORFsearch.filter
    """
}
