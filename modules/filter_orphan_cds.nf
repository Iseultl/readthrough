process FILTER_ORPHAN_CDS {
    tag { "filter_orphan_cds" }
    label 'python'
    cpus 1
    memory '4GB'
    time '30m'

    input:
    path gff

    output:
    path("*.orphan_cds_removed.gff")

    script:
    def out_name = "${gff.baseName}.orphan_cds_removed.gff"
    """
    filter_orphan_cds.py \
        --gff ${gff} \
        --output ${out_name}

    rm ${gff}
    """
}
