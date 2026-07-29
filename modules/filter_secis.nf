process FILTER_SECIS {
    publishDir "${params.output_dir}/secis_final_gff", mode: 'copy'
    container 'ubuntu:22.04'
    memory '4GB'
    input:
    path "all_secis_combined.gff"

    output:
    path "filtered_secis.gff"

    script:
    """
    awk '{
        infernal=0; covels=0;
        for (i=1; i<=NF; i++) {
            if (\$i ~ /^Infernal_score=/) {
                split(\$i,a,"=");
                infernal=a[2]+0;
            }
            if (\$i ~ /^Covels_score=/) {
                split(\$i,b,"=");
                covels=b[2]+0;
            }
        }
        if (infernal >= 15 && covels >= 15) print
    }' all_secis_combined.gff > filtered_secis.gff
    """
}
