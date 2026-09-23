process SPLIT_IF_TOO_LARGE {
    tag { "split_if_needed_${input_file.baseName}" }
    cpus 1
    memory '4GB'
    label 'splitfasta'

    input:
    path input_file

    output:
    path "*.fa.split/*.fa", emit: split_fasta

    script:
    """
    seqkit split2 --by-size 40000 ${input_file}
    """
}
