process RELOCATE_TRANSCRIPTS {
    tag { "relocate_${cleaned_gtf}" }
    label 'python'
    cpus 1
    time '1h'
    memory '8GB'

    input:
    path cleaned_gtf
    
    output:
    path "relocated_gtf/*", emit: relocated_gtf
    
    script:
    """
    mkdir -p relocated_gtf

    relocate_transcripts.py \
        --gff ${cleaned_gtf} \
        --output relocated_gtf/relocated.gtf

    rm ${cleaned_gtf}
    """
}
