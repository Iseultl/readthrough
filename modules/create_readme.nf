process CREATE_README {
    tag "Creating README"
    publishDir "${params.output_dir}", mode: 'copy' 
    label 'samtools'

    input:
    val species_name
    path genome
    path annotation_gtf

    output:
    path "README.txt"

    script:
    """
    # Create README file with genome length and number of annotated genes
    printf "${species_name}\n" > README.txt

    # Create FASTA index
    samtools faidx ${genome}

    # Sum sequence lengths (column 2 of the .fai file)
    awk '{sum+=\$2} END {print "Genome Length: " sum}' ${genome}.fai >> README.txt


    # Extract number of annotated genes from GTF file
    cut -f3 ${annotation_gtf} | grep -c 'gene'  | awk '{print "Number of Annotated Genes: " \$1}' >> README.txt
    """
}