process CREATE_PAIRING_FILE {
    tag "create_pairing_file"
    publishDir "chromosome_pairs", mode: 'copy'
    label 'python'

    input:
    val all_pairs

    output:
    path "chromosome_pairs.tsv"

    script:
    """
    > all_pairs.txt
    ${all_pairs.collect { chr, fasta, gff -> "echo '${chr}\t${fasta}\t${gff}'" }.join('\n')}
    """
}