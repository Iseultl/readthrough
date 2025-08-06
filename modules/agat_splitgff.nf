process AGAT_SPLITGFF {
    tag "${genome_gtf.baseName}"

    label 'agat'

    input:
    path genome_gtf

    output:
    path "agat_gff2gtf/*.gff", emit: gff_files

    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Create output directory
    mkdir -p agat_gff2gtf
    
    # Split GFF by chromosome/sequence
    awk -F '\\t' -v outdir="agat_gff2gtf" '
        \$0 !~ /^#/ && NF >= 9 {
            print > outdir "/" \$1 ".gff"
        }
    ' "${genome_gtf}"
    
    # Remove single-line GFF files (header-only or empty files)
    find agat_gff2gtf -type f -name "*.gff" | while read -r file; do
        if [ \$(wc -l < \"\$file\") -le 1 ]; then
            rm -f "\$file"
        fi
    done
    
    """
}