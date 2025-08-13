process RELOCATE_TRANSCRIPTS {
    tag { "relocate_${cleaned_gtf}" }
    
    publishDir "${params.output_dir}/relocated_gtf", mode: 'copy'
    label 'python'
    cpus 1
    memory '8GB'
    time '1h'
    
    input:
    path cleaned_gtf
    
    output:
    path "relocated_gtf/*"
    
    script:
    """
    mkdir -p relocated_gtf

    relocate_transcripts.py \
        --gff ${cleaned_gtf} \
        --output relocated_gtf/relocated.gtf
    """
}
