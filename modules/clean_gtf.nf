process CLEAN_GTF {
    tag { "clean_gtf_${genome_gtf}" }

    cpus 1

    input:
    path(genome_gtf)

    output:
    path "*.cleaned.gtf", emit: cleaned_gtf

    script:
    """
    #!/bin/bash
    # Run the GTF cleaning script
    bash clean_gtf.sh ${genome_gtf}
    """
}
