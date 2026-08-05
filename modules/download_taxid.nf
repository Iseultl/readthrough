process DOWNLOAD_TAXID {
    tag "Downloading TaxID"
    memory '4GB' 
    errorStrategy 'retry'
    maxRetries 2

    input:
    val taxid
    val annotation_url
    val fasta_url

    output:
    path("genome_${taxid}.fasta"), emit: fasta
    path("genome_${taxid}.gff"), emit: gff

    script:
    """
    #!/bin/bash
    set -euo pipefail

    download_file.py --taxid ${taxid} --annotation-url "${annotation_url}" --fasta-url "${fasta_url}" --retry-log download_retry.tsv --download-delay 8 --retry-delay 15 --max-attempts 3

    annocli alias annotation.gff.gz annotation.fasta.gz --output annotation.aliasMatch.gff.gz

    gunzip -c "annotation.aliasMatch.gff.gz" > "genome_${taxid}.gff"
    gunzip -c "annotation.fasta.gz" > "genome_${taxid}.fasta"

    rm -rf annotation.*
    """
}