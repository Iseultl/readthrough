// modules/gffread.nf

process GFFREAD {
    tag { "gffread_${id}" }
    label 'gffread'
    cpus 1
    memory '4GB'
    time '1h'

    input:
    tuple path(gtf), path(fasta)
    
    output:
    path("gffread_out/*"), emit: transcripts
    path("gffread_out"), emit: gffread_dir
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Create output directory if it doesn't exist
    mkdir -p gffread_out

    # Run gffread
    gffread -F -w transcripts.fa -g ${fasta} ${gtf}
    
    # Process transcripts
    awk '/^>/ {sub(/^>/, ">"); print \$1; next} {print}'  transcripts.fa > transcripts_clean_${id}.fa 
    # Move output to output directory
    mv transcripts_clean_${id}.fa gffread_out/
    
    """
}
