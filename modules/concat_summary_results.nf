process CONCAT_SUMMARY_RESULTS {
    input:
    path summary_tables

    output:
    path "ORFsearch.result"

    script:
    """
    cat ${summary_tables.join(' ')} > ORFsearch.result
    """
}
