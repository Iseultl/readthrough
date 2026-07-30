process MERGE_GFF {
    container 'ubuntu:22.04'
    memory '4GB'
    
    input:
    path gff_files

    output:
    path "all_secis_combined.gff"

    script:
    """
    cat *.gff > all_secis_combined.gff
    """
}