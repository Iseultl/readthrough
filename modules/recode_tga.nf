process RECODE_TGA {
    tag { "recode_tga" }
    label 'python'
    cpus 1
    time '1h'
    memory '8GB'
    
    input:
    path(transcript_fasta)
    val limit
    
    output:
    path("*.fa")
    
    script:
    def base = transcript_fasta.getBaseName()
    def out_name = "${base}_recoded.fa"
    """
    recode_any_TGA.py \
        --fasta ${transcript_fasta} \
        --recodon TGC \
        --limit ${limit} \
        --output ${out_name}

    rm ${transcript_fasta}
    """
}
