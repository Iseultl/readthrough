// modules/handle_zipped_input.nf

process UNZIP_IF_NEEDED {
    tag { "${file_type}" }
    memory '4GB'
    
    input:
    tuple val(file_type), path(input_file)
    
    output:
    tuple val(file_type), path("${output_name}"), emit: unzipped_file
    
    script:
    is_gzipped = input_file.name.endsWith('.gz')
    output_name = is_gzipped ? input_file.name - '.gz' : input_file.name
    
    if (is_gzipped) {
        """
        gunzip -c ${input_file} > ${output_name}
        """
    } else {
        """
        cp ${input_file} ${output_name}
        """
    }
}
