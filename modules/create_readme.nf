process CreateReadme {
    tag "Creating README"
    publishDir "${params.output_dir}", mode: 'copy' 
    label 'samtools'

    input:
    path genome
    path annotation_gtf

    output:
    path "README.txt"

    script:
    """
    # Calculate genome length using samtools
    samtools faidx ${genome} | awk '{sum+=\$2} END {print "Genome Length: " sum}' > README.txt

    # Extract number of annotated genes from GTF file
    cut -f3 ${annotation_gtf} | grep -c 'gene'  | awk '{print "Number of Annotated Genes: " \$1}' >> README.txt
    """
}