process RELOCATE_TRANSCRIPTS {
    tag { "relocate_${cleaned_gtf}" }
    
    publishDir "${params.output_dir}/relocated_gtf", mode: 'copy'
    label 'python'
    cpus 1
    memory '8GB'
    time '1h'
    
    input:
    path cleaned_gtf
    
    output:
    path "relocated_gtf/*"
    
    script:
    """
    mkdir -p relocated_gtf

    awk 'BEGIN{OFS="\\t"} 
     /^#/ {print; next} 
     {
       tid = \$9;
       split(tid,a,"-");
       gid = a[1];
       \$9 = "gene_id \\"" gid "\\"; transcript_id \\"" tid "\\";";
       print
     }' ${cleaned_gtf} > ${cleaned_gtf}.fixed

    relocate_transcripts.py \
        --gff ${cleaned_gtf}.fixed \
        --output relocated_gtf/relocated.gtf
    """
}
