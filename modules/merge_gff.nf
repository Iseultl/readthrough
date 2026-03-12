process MERGE_GFF {
    publishDir "${params.output_dir}/secis_final_gff", mode: 'copy'
    container 'ubuntu:22.04'

    input:
    path gff_files

    output:
    path "all_secis_combined.gff"

    script:
    """
    cat *.gff > all_secis_combined.gff
    """
}