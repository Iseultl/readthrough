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
    awk -F'\\t' 'BEGIN { OFS="\\t" }

    # Leave comments untouched
    /^#/ {
        print
        next
    }

    {
        original_attrs = \$9

        gene_id = ""
        transcript_id = ""
        id = ""
        parent = ""

        n = split(\$9, attrs, ";")

        for (i = 1; i <= n; i++) {

            # Remove leading/trailing whitespace
            attr = attrs[i]
            gsub(/^[[:space:]]+|[[:space:]]+\$/, "", attr)

            # ----------------------------
            # GTF-style attributes
            # ----------------------------
            if (attr ~ /^gene_id[[:space:]]+/)
                gene_id = attr

            else if (attr ~ /^transcript_id[[:space:]]+/)
                transcript_id = attr

            # ----------------------------
            # GFF3-style attributes
            # ----------------------------
            else if (attr ~ /^ID=/)
                id = attr

            else if (attr ~ /^Parent=/)
                parent = attr
        }

        out = ""

        # ----------------------------------------
        # GTF input
        # Keep gene_id and transcript_id
        # ----------------------------------------
        if (gene_id != "" || transcript_id != "") {

            if (gene_id != "")
                out = gene_id

            if (transcript_id != "")
                out = (out == "" ? transcript_id : out "; " transcript_id)

            out = out ";"
        }

        # ----------------------------------------
        # GFF3 input
        # Keep ID and Parent
        # ----------------------------------------
        else if (id != "" || parent != "") {

            if (id != "")
                out = id

            if (parent != "")
                out = (out == "" ? parent : out ";" parent)
        }

        # ----------------------------------------
        # Unknown/unrecognised format
        # Do NOT destroy attributes
        # ----------------------------------------
        else {
            out = original_attrs
        }

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