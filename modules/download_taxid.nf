process DOWNLOAD_TAXID {
    tag "Downloading TaxID"
    memory '4GB' 

    input:
    val taxid

    output:
    path("genome_${taxid}.fasta"), emit: fasta
    path("genome_${taxid}.gff"), emit: gff

    script:
    """
    #!/bin/bash
    set -euo pipefail

    annocli download --ref-only --taxids "${taxid}" --add-asm --fix-alias --output-dir annotation_downloads

    echo "Downloaded files for taxid ${taxid}."

    fasta_source=\$(find annotation_downloads -type f -iname "*.fna.gz" -o -iname "*.fasta.gz" | head -n 1)
    alias_output="\$(find annotation_downloads -type f -iname *.aliasMatch.g*.gz" | head -n 1)

    if [[ -z "\${alias_output}" || -z "\${fasta_source}" ]]; then
        echo "Could not locate downloaded annotation or assembly files under annotation_downloads/" >&2
        find annotation_downloads -type f | sort >&2
        exit 1
    else
        echo "Located annotation file: \${alias_output}"
        echo "Located assembly file: \${fasta_source}"
    fi

    gunzip -c "\${fasta_source}" > "genome_${taxid}.fasta"
    gunzip -c "\${alias_output}" > "genome_${taxid}.gff"

    rm -rf annotation_downloads
    """
}