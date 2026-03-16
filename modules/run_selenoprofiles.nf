#!/usr/bin/env nextflow

process runSelenoprofiles {
    tag "${genomeFile}"
    label 'selenoprofiles'
    publishDir "${params.output_dir}/selenoprofiles", mode: 'copy'

    input:
    path genome
    path reference_annotation
    species_name

    output:
    "selenoprofiles_output/all_predictions.gtf"
    "selenoprofiles_output/annotation_result.csv"

    script:
    """
    mkdir -p selenoprofiles_output
    selenoprofiles  -o selenoprofiles_output  -t ${genome}  -s "${species_name}"  -p eukarya -output_gtf_file selenoprofiles_output/all_predictions.gtf

    selenoprofiles assess -s selenoprofiles_output/all_predictions.gtf \\
        -e ${reference_annotation} \\
        -f ${genome} \\
        -o selenoprofiles_output/annotation_result.csv 
    """
}

