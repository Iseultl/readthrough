// modules/gffread.nf

process GFFREAD {
    tag { "gffread" }
    label 'gffread'
    cpus 1
    memory '4GB'
    time '1h'

    input:
    path(gtf)
    path(fasta)
    
    output:
    path("gffread_out/*"), emit: transcripts
    path("gffread_out"), emit: gffread_dir
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    . /usr/local/env-activate.sh
    
    # Create output directory if it doesn't exist
    mkdir -p gffread_out

    # Run gffread
    gffread -F -w transcripts.fa -g ${fasta} ${gtf}
    
    # Process transcripts. Headers are the first word of each title.
    # Seblastian (SECISSEARCH) rejects identifiers longer than 63 characters, so
    # any over-length id is deterministically shortened: first 55 chars + '_' +
    # 6-digit rolling hash of the full id (stays unique, stays <= 63).
    awk 'BEGIN { for (c = 1; c < 256; c++) ORD[sprintf("%c", c)] = c }
         /^>/ { sub(/^>/, ">"); id = \$1;
                if (length(id) > 63) {
                    h = 0
                    for (c = 1; c <= length(id); c++) h = (h * 31 + ORD[substr(id, c, 1)]) % 1000000
                    id = substr(id, 1, 55) "_" sprintf("%06d", h)
                }
                print id; next }
         { print }' transcripts.fa > transcripts_clean.fa
    # Move output to output directory
    mv transcripts_clean.fa gffread_out/
    
    rm ${fasta} ${gtf}
    """
}
