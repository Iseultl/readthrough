process AGAT_GFF2GTF {
    tag "${gff_file.baseName}"
    label 'agat'

    input:
    path gff_file

    output:
    path "${gff_file.baseName}.gtf", emit: gtf_files

    script:
    """
    #!/bin/bash
    agat_convert_sp_gff2gtf.pl \\
            --gff "${gff_file}" \\
            --gtf_version 3 \\
            --output "${gff_file.baseName}.gtf"
    """
}