process CREATE_PARAM_FILE {
    tag "Creating Parameter File"
    memory '4GB' 
    label 'param_file'
    errorStrategy 'retry'
    maxRetries 2

    input:
    path fasta
    path gtf
    val species_name

    output:
    path("${species_name}.param"), emit: param_file
    
    script:
    """
    #!/bin/bash
    set -euo pipefail

    bash generate_param.sh "${species_name}" "${fasta}" "${gtf}" .
    """
}