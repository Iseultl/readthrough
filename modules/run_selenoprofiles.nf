#!/usr/bin/env nextflow

process RUN_SELENOPROFILES {
    tag "${genome.baseName}"
    label 'selenoprofiles'
    publishDir "${params.output_dir}/selenoprofiles", mode: 'copy'


    input:
    path genome
    path reference_annotation
    val species_name

    output:
    path "selenoprofiles_output/all_predictions.gtf"
    path "selenoprofiles_output/annotation_result.csv"

    script:
    """
    mkdir -p selenoprofiles_output
    
    selenoprofiles -setup
    selenoprofiles -download

    selenoprofiles  -o selenoprofiles_output  -t ${genome}  -s "${species_name}"  -p eukarya -output_gtf_file -temp temp_folder selenoprofiles_output/all_predictions.gtf

    selenoprofiles assess -s selenoprofiles_output/all_predictions.gtf \\
        -e ${reference_annotation} \\
        -f ${genome} \\
        -o selenoprofiles_output/annotation_result.csv 
    """
}

