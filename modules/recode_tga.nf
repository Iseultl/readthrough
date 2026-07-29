process RECODE_TGA {
    tag { "recode_tga" }
    label 'python'
    cpus 1
    time '1h'
    
    input:
    path(transcript_fasta)
    
    output:
    path("*.fa")
    
    script:
    def base = transcript_fasta.getBaseName()
    def out_name = "${base}_recoded.fa"
    """
    recode_any_TGA.py \
        --fasta ${transcript_fasta} \
        --recodon TGC \
        --output ${out_name}
    """
}
