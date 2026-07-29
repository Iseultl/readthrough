process SELENOPROFILES_ANALYSIS {
    tag "${genome_id}"
    label 'selenoprofiles'
    publishDir "${params.output_dir}/${genome_id}/selenoprofiles", mode: 'copy'
    
    input:
    tuple val(genome_id), path(genome_fasta), path(gff)
    path(params.candidates)
    val(params.species_name)
    
    output:
    tuple val(genome_id), path("SecORFsearch_profiles_output")
    
    script:
    """
    mkdir -p SecORFsearch_profiles_output

    if [ -f ~/.selenoprofiles_config.txt ]; then
        echo "Config file exists"
    else
        echo "Config file does not exist"
        selenoprofiles -setup
    fi
    
    for f in ${params.candidates}/*.fasta; do
        if [ ! -f "\$f" ]; then
            echo "Candidate proteins file \$f does not exist." >&2
            exit 1
        fi
        BASENAME=\$(basename "\$f" .fasta)

        #Replace any "." in the basename with "_" to avoid issues with Selenoprofiles
        BASENAME="\${BASENAME//./_}"

        # Create a profile
        selenoprofiles build -i "\$f" -o "profile_\$BASENAME.fasta" -y

        # Run Selenoprofiles search
        selenoprofiles -o "\$BASENAME" -t ${genome_fasta} -p "profile_\$BASENAME.fasta" -tblastn -s ${params.species_name} -output_gtf_file SecORFsearch_profiles_output/\$BASENAME.gtf -temp temp_folder

        # Clean up stale database locks for the next iteration
        if [[ -f "\$BASENAME/${params.species_name}.genome/results.sqlite" ]]; then
            run_selenoprofiles database -i "\$BASENAME/${params.species_name}.genome/results.sqlite" -clean 2>/dev/null || true
        fi

        # Clean up & extract only the relevant output folder
        mkdir -p "SecORFsearch_profiles_output/\$BASENAME"
        cp -r "\$BASENAME"/*/output/* "SecORFsearch_profiles_output/\$BASENAME/" || true
        rm -rf "\$BASENAME"
        rm -f "profile_\${BASENAME}.fasta"* 2>/dev/null || true
    done
    
    """
}