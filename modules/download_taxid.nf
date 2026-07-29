process DOWNLOAD_TAXID {
    tag "Downloading TaxID"
    publishDir "${params.output_dir}", mode: 'copy' 

    input:
    val taxid

    output:
    tuple path("genome_${taxid}.fasta.gz"), path("genome_${taxid}.gff.gz")

    script:
    """
    #!/bin/bash
    set -euo pipefail

    annocli download --ref-only --taxids "${taxid}" --add-asm --mode links | head -n 1 | bash 
    echo "Downloaded files for taxid ${taxid}."

    annotation_source=\$(find annotation_downloads -type f -name "*.gff.gz" ! -name "*.aliasMatch.gff.gz" | head -n 1)
    fasta_source=\$(find annotation_downloads -type f -name "*.fna.gz" | head -n 1)

    if [[ -z "\${annotation_source}" || -z "\${fasta_source}" ]]; then
        echo "Could not locate downloaded annotation or assembly files under annotation_downloads/" >&2
        find annotation_downloads -type f | sort >&2
        exit 1
    fi

    alias_output="\${annotation_source%.gff.gz}.aliasMatch.gff.gz"

    annocli alias "\${annotation_source}" "\${fasta_source}" --output "\${alias_output}"
    rm -f "annotation_downloads/"*.aliasMappings.tsv 2>/dev/null || true

    cp "\${fasta_source}" "genome_${taxid}.fasta.gz"
    cp "\${alias_output}" "genome_${taxid}.gff.gz"

    rm -rf annotation_downloads
    """
}