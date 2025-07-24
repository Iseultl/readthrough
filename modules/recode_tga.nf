process RECODE_TGA {
    tag { "recode_tga" }
    
    publishDir "recoded_transcripts", mode: 'copy'
    label 'python'
    cpus 1
    memory '8GB'
    time '1h'
    
    input:
    path(transcript_fasta)
    
    output:
    path("recoded.fa")
    
    script:
    """
    recode_any_TGA.py \
        --fasta ${transcript_fasta} \
        --recodon TGC \
        --output recoded.fa
    """
}
