// modules/gffread.nf

process GFFREAD {
    tag { "gffread" }

    publishDir "gffread_out", mode: 'copy'
    label 'gffread'

    input:
    tuple val(id), path(gtf), path(fasta)
    
    output:
    path("gffread_out/*")
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Create output directory if it doesn't exist
    mkdir -p gffread_out
    
    # Run gffread
    gffread -F -w transcripts.fa -g ${fasta} ${gtf}
    
    # Process transcripts
    awk '/^>/ {sub(/^>/, ">"); print \$1; next} {print}'  transcripts.fa > transcripts_clean.fa
    
    # Move output to output directory
    mv transcripts_clean.fa gffread_out/
    
    """
}
