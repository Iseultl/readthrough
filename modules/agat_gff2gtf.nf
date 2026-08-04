process AGAT_GFF2GTF {
    tag "${gff_file.baseName}"
    label 'agat'
    memory '4GB'
    
    input:
    path gff_file

    output:
    path "${gff_file.baseName}.gtf", emit: gtf_file

    script:
    """
    #!/bin/bash
    awk -F'\t' 'BEGIN{OFS="\t"}
    {
        gene=""
        transcript=""

        if (match(\$9, /gene_id "[^"]+"/))
            gene = substr(\$9, RSTART, RLENGTH)

        if (match(\$9, /transcript_id "[^"]+"/))
            transcript = substr(\$9, RSTART, RLENGTH)

        if (gene && transcript)
            \$9 = gene_id " " gene "; " transcript_id " " transcript ";"

        print
    }' ${gff_file} > ${gff_file.baseName}.temp.gtf

    agat_convert_sp_gff2gtf.pl \\
            --gff "${gff_file.baseName}.temp.gtf" \\
            --gtf_version 3 \\
            --output "${gff_file.baseName}.gtf"

    rm ${gff_file.baseName}.temp.gtf ${gff_file}
    """
}