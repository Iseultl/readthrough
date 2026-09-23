// modules/gffread_chr.nf

process GFFREAD_CHR {
    tag { "gffread_${id}" }
    label 'gffread'
    cpus 1
    memory '8GB'
    time '1h'

    input:
    tuple val(id), path(gtf), path(fasta)
    
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
    
    # Process transcripts. Headers are the first word of each title.
    # Seblastian (SECISSEARCH) rejects identifiers longer than 63 characters
    # (Cyanidiococcus 28725643 failed on a 64-char GenBank gnl|WGS id), so any
    # over-length id is deterministically shortened: first 55 chars + '_' +
    # 6-digit rolling hash of the full id (stays unique, stays <= 63).
    awk 'BEGIN { for (c = 1; c < 256; c++) ORD[sprintf("%c", c)] = c }
         /^>/ { sub(/^>/, ">"); id = \$1;
                if (length(id) > 63) {
                    h = 0
                    for (c = 1; c <= length(id); c++) h = (h * 31 + ORD[substr(id, c, 1)]) % 1000000
                    id = substr(id, 1, 55) "_" sprintf("%06d", h)
                }
                print id; next }
         { print }' transcripts.fa > transcripts_clean_${id}.fa
    # Move output to output directory
    mv transcripts_clean_${id}.fa gffread_out/
    
    """
}