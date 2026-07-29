process SPLIT_IF_TOO_LARGE {
    tag { "split_if_needed_${input_file.baseName}" }
    cpus 1
    memory '4GB'
    maxForks 1
    label 'splitfasta'

    input:
    path input_file

    output:
    path "*.fa", emit: split_fasta

    script:
    """
    echo "Processing file: ${input_file}"
    echo "[\$(date '+%F %T')] split_if_too_large: starting size check for ${input_file}"
    # Use seqkit stats with tabular output and headers to robustly get sequence count
    echo "[\$(date '+%F %T')] split_if_too_large: running seqkit stats"
    SEQ_COUNT=\$(seqkit stats -T ${input_file} | awk 'NR==2 {print \$4}')
    echo "Sequence count: \$SEQ_COUNT"
    THRESHOLD=40000
    echo "[\$(date '+%F %T')] split_if_too_large: threshold=\$THRESHOLD"

    if [ -n "\$SEQ_COUNT" ] && [ "\$SEQ_COUNT" -gt "\$THRESHOLD" ]; then
        echo "Splitting file as sequence count \$SEQ_COUNT is greater than threshold \$THRESHOLD"
        echo "[\$(date '+%F %T')] split_if_too_large: entering split branch"
        awk -v threshold="\$THRESHOLD" -v base="${input_file.baseName}" '
            BEGIN {
                part = 1
                seqs_in_part = 0
                out = sprintf("%s.part_%03d.fa", base, part)
                print "split_if_too_large: opening " out > "/dev/stderr"
            }
            /^>/ {
                if (seqs_in_part >= threshold) {
                    print "split_if_too_large: closing " out " after " seqs_in_part " sequences" > "/dev/stderr"
                    close(out)
                    part++
                    seqs_in_part = 0
                    out = sprintf("%s.part_%03d.fa", base, part)
                    print "split_if_too_large: opening " out > "/dev/stderr"
                }
                seqs_in_part++
            }
            {
                print > out
            }
            END {
                print "split_if_too_large: final close " out " after " seqs_in_part " sequences" > "/dev/stderr"
                close(out)
            }
        ' ${input_file}
        echo "[\$(date '+%F %T')] split_if_too_large: split branch completed"
    else
        echo "No split needed"
        echo "[\$(date '+%F %T')] split_if_too_large: entering no-split branch"
        cp ${input_file} ${input_file.baseName}.part_001.fa
        echo "[\$(date '+%F %T')] split_if_too_large: no-split branch completed"
    fi
    echo "[\$(date '+%F %T')] split_if_too_large: task finished"
    """
}

