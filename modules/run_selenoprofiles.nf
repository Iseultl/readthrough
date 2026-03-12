#!/usr/bin/env nextflow

process runSelenoprofiles {
    tag "${genomeFile}"
    label 'selenoprofiles'
    publishDir "${params.output_dir}/selenoprofiles", mode: 'copy'

    input:
    path genome
    species_name

    output:
    "all_predictions.gtf"

    script:
    """
    mkdir -p selenoprofiles_output
    selenoprofiles  -o selenoprofiles_output  -t ${genome}  -s "${species_name}"  -p eukarya -output_gtf_file all_predictions.gtf
    """
}

