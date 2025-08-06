process SPLIT_RECODED {
    tag { "split_${input_file}" }
    
    container "quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0"
    cpus 1
    memory '8GB'
    
    input:
    path(recoded_file)
    
    output:
    path("split_recoded_transcripts/*.fa"), emit: split_recoded_dir
    
    script:
    """
    #!/bin/bash
    set -euo pipefail 
    mkdir -p split_recoded_transcripts
    seqkit split -i -O split_recoded_transcripts ${recoded_file}
    """
}