// modules/splitfasta.nf

process SPLITFASTA {
    tag { "${genome_fasta.getBaseName()}" }

    label 'splitfasta'
    
    input:
    path(genome_fasta)
    
    output:
    path "split_chr/*.fa", emit: split_chr
    
    script:
    """
    mkdir -p split_chr
    seqkit split -i -O split_chr ${genome_fasta}
    """
}