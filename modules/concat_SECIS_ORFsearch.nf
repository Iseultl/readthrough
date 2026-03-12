#!/usr/bin/env nextflow

process concat_SECIS_ORFsearch {
    tag "Concatenating SECIS and ORFsearch results"

    input:
    path merged_secis_gff
    path filtered_secis
    path orfsearch_result

    output:
    path "concatenated_results.tsv"

    script:
    """
    # Python script to concatenate the results 
    """
}