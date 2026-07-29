// modules/splitfasta.nf

process SPLITFASTA {
    tag { "${genome_fasta.getBaseName()}" }
    memory '6GB'
    label 'splitfasta'
    
    input:
    path(genome_fasta)
    
    output:
    path "split_chr/*.fa", emit: split_chr
    
    script:
    """
    mkdir -p split_chr
    seqkit split -i -O split_chr ${genome_fasta}

    for f in split_chr/*; do
        mv "\$f" "\${f%.*}.fa"
    done
    """
}