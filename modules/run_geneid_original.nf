process RUN_GENEID_ORIGINAL {
    tag { "run_geneid_original_${input_file}" }
    
    publishDir "${params.output_dir}/geneid_original_predictions", mode: 'copy'
    cpus 1
    memory '8GB'
    
    input:
    path(input_file)
    path(param_file)
    
    output:
    path("geneid_original_predictions/*.txt")
    
    script:
    """
    # Run geneid on the transcripts
    run_geneid_original.sh \
        ${input_file} \
        geneid_original_predictions \
        ${param_file}
    """
}
