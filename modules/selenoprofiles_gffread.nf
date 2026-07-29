// modules/selenoprofiles_gffread.nf

process SELENOPROFILES_GFFREAD {
    publishDir "${params.output_dir}/gffread_out", mode: 'copy'
    tag { "gffread_${id}" }
    label 'gffread'
    cpus 1
    memory '8GB'
    time '1h'

    input:
    tuple val(id), path(gtf), path(fasta)
    
    output:
    tuple val(genome_id), path("selenoprofiles_output/*")
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Create output directory if it doesn't exist
    mkdir -p gffread_out

    # Run gffread
    gffread -F -w transcripts.fa -g ${fasta} ${gtf} --w-add 600
    
    # Process transcripts
    awk '/^>/ {sub(/^>/, ">"); print \$1; next} {print}'  transcripts.fa > transcripts_clean_${id}.fa 
    # Move output to output directory
    mv transcripts_clean_${id}.fa selenoprofiles_output/
    
    """
}
