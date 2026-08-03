// modules/handle_zipped_input.nf

process UNZIP_IF_NEEDED {
    tag { "unzip" }
    memory '4GB'
    
    input:
    tuple path(fasta), path(gff)
    val ${params.species_taxid}
    
    output:
    tuple path("*.fasta"), path("*.gff")
    
    script:
    """
    #!/bin/bash
    gunzip -c ${fasta} > genome_${params.species_taxid}.fasta
    gunzip -c ${gff} > genome_${params.species_taxid}.gff
    """
}
