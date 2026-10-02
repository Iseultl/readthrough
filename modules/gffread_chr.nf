// modules/gffread_chr.nf
//
// Batched (main_protists.nf only, since 2026-09-28): one job per hash bucket
// of scaffolds (see batchByScaffold in main_protists.nf), looping over the
// paired (gtf, fasta) scaffolds inside the job.
//
// The gtf and fasta files of the bucket arrive as ONE merged list (see
// modules/scaffold_batches.nf for the channel rules), so the pair for each
// scaffold is matched BY FILENAME inside the job: every gtf is
// "<scaffold>.gtf" and its genome part is "<stem>.part_<scaffold>.fa".
//
// Emits the same per-scaffold gffread_out/transcripts_clean_<scaffold>.fa
// files as the original one-scaffold-per-job version, so everything
// downstream (recode, split, geneid, secissearch, logos) is unaffected.

process GFFREAD_CHR {
    tag { "gffread_batch_${bucket}" }
    label 'gffread'
    cpus 1
    memory '8GB'
    // No explicit time: a bucket may hold hundreds of scaffolds (Conticribra
    // 7,719 / 24 buckets ~ 322/job); the cluster profile default is 8h.

    input:
    tuple val(bucket), path(files)

    output:
    path("gffread_out/*"), emit: transcripts
    path("gffread_out"), emit: gffread_dir

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Create output directory if it doesn't exist
    mkdir -p gffread_out

    for gtf in *.gtf; do
        scaffold="\$(basename "\$gtf" .gtf)"

        # Find the matching genome part for this scaffold (flat root staging).
        fasta=""
        for fa in *.fa; do
            case "\$fa" in
                *.part_"\$scaffold".fa) fasta="\$fa"; break ;;
            esac
        done
        if [ -z "\$fasta" ]; then
            echo "ERROR: no genome part matching scaffold '\$scaffold'" >&2
            exit 1
        fi

        # Run gffread
        gffread -F -w transcripts.fa -g "\$fasta" "\$gtf"

        # Process transcripts. Headers are the first word of each title.
        # Seblastian (SECISSEARCH) rejects identifiers longer than 63 characters
        # (Cyanidiococcus 28725643 failed on a 64-char GenBank gnl|WGS id), so any
        # over-length id is deterministically shortened: first 55 chars + '_' +
        # 6-digit rolling hash of the full id (stays unique, stays <= 63).
        awk 'BEGIN { for (c = 1; c < 256; c++) ORD[sprintf("%c", c)] = c }
             /^>/ { sub(/^>/, ">"); id = \$1;
                    if (length(id) > 63) {
                        h = 0
                        for (c = 1; c <= length(id); c++) h = (h * 31 + ORD[substr(id, c, 1)]) % 1000000
                        id = substr(id, 1, 55) "_" sprintf("%06d", h)
                    }
                    print id; next }
             { print }' transcripts.fa > "gffread_out/transcripts_clean_\${scaffold}.fa"

        rm transcripts.fa
    done
    """
}
