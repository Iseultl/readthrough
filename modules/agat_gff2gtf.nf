process AGAT_GFF2GTF {
    tag "${gff_file.baseName}"
    label 'agat'
    memory '6GB'
    
    input:
    path gff_file

    output:
    path "${gff_file.baseName}.gtf", emit: gtf_file

    script:
    """
    #!/bin/bash
    awk -F'\\t' 'BEGIN{OFS="\\t"}
    {
        n = split(\$9, attrs, ";")

        id = ""
        parent = ""

        for (i = 1; i <= n; i++) {
            gsub(/^[[:space:]]+|[[:space:]]+\$/, "", attrs[i])

            if (attrs[i] ~ /^ID=/)
                id = attrs[i]
            else if (attrs[i] ~ /^Parent=/)
                parent = attrs[i]
        }

        out = ""
        if (id != "")
            out = id
        if (parent != "")
            out = (out == "" ? parent : out ";" parent)

        \$9 = out
        print
    }' "${gff_file}" > "${gff_file.baseName}.temp.gff" 

    agat_convert_sp_gff2gtf.pl \\
            --gff "${gff_file.baseName}.temp.gff" \\
            --gtf_version 3 \\
            --output "${gff_file.baseName}.gtf"

    if [ ! -s "${gff_file.baseName}.gtf" ]; then
        echo "Error: GTF file is empty or not created."
        mv "${gff_file.baseName}.temp.gff" "${gff_file.baseName}.gtf"
        rm "${gff_file}"
    else
        echo "GTF file created successfully."
        rm "${gff_file.baseName}.temp.gff" "${gff_file}"
    fi
    
    """
}