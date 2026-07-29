process CONCATENATE_GTFS {
    input:
    path gtf_files
    path gffread_out

    output:
    path "concatenated.gtf"

    script:
    """
    #!/bin/bash
    echo ${gtf_files}
    cat ${gtf_files.join(' ')} > concatenated.gtf
    """
}