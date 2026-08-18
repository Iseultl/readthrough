#!/usr/bin/env nextflow

process RUN_SELENOPROFILES {
    tag "${genome.baseName}"
    label 'selenoprofiles'
    publishDir "${params.output_dir}/selenoprofiles", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 2
    memory { 8.GB * task.attempt }
    time '3h'

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

    if [ -f ~/.selenoprofiles_config.txt ]; then
        echo "Config file exists"
    else
        echo "Config file does not exist"
        selenoprofiles -setup
    fi

    if [ ! -d ~/selenoprofiles_data/selenoprotein_profiles ]; then
        printf '\\n' | selenoprofiles -download -y
    fi

    export SELENOPROFILES_DATA_DIR=~/selenoprofiles_data

    selenoprofiles  -o selenoprofiles_output  -t ${genome}  -s "${species_name}"  -p eukarya -output_gtf_file selenoprofiles_output/all_predictions.gtf -temp temp_folder

    if grep -qi 'Selenocysteine' selenoprofiles_output/all_predictions.gtf; then
        selenoprofiles assess -s selenoprofiles_output/all_predictions.gtf \\
            -e ${reference_annotation} \\
            -f ${genome} \\
            -o selenoprofiles_output/annotation_result.csv
    else
        echo "No Selenocysteine entries found in all_predictions.gtf; skipping selenoprofiles assess."
        printf 'status,details\n' > selenoprofiles_output/annotation_result.csv
        printf 'skipped,No Selenocysteine entries found in all_predictions.gtf\n' >> selenoprofiles_output/annotation_result.csv
    fi
    """
}

