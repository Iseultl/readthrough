// modules/scaffold_batches.nf
//
// Batched variants of the per-scaffold stages, used by main_protists.nf for
// genomes with many small scaffolds (e.g. Cylindrotheca_closterium: 2,534
// scaffolds -> ~15,000 short SLURM jobs, dominated by submission/container
// overhead; Conticribra: 7,719 scaffolds -> ~46,000).
//
// Each process receives a hash bucket (0..scaffold_batches-1) plus the list
// of per-scaffold files assigned to it, and loops over the files inside ONE
// job. This trades many short jobs for a few longer ones WITHOUT increasing
// peak memory (the loop is sequential: one scaffold at a time).
//
// Every process still emits the same per-scaffold output files (same names,
// same layout) as the original one-job-per-scaffold modules, so all
// downstream name/id derivation (geneid part names, SELECT_INTERESTING /
// GET_ORIGINAL_PREDICTIONS replaceAll chains, ...) is unaffected.
//
// CHANNEL/STAGING RULES (verified empirically on this cluster, NF 26.04.6):
//  - a SINGLE list `path` input WITHOUT stageAs is staged flat into the work
//    dir root, base names preserved (same pattern as the production
//    EXTRACT_SEQUENCE_LOGOS collect() input) -> the loops below use root globs.
//  - NEVER add stageAs to a list input: 'stageAs: "in/*.gtf"' stages the
//    files as in/1.gtf, in/2.gtf, ... (numbered) and the scaffold names are
//    lost (this silently broke the gtf<->fasta combine once).
//  - groupTuple() is COLUMN-WISE: for (a,b,c) tuples it emits
//    (a, [b...], [c...]) — extract item[1]/item[2] as separate file lists
//    (the "input file name collision: nfs, scratch01" errors seen in
//    2026-09-28 toy runs were NOT a staging bug: indexing it[0]/it[1] on a
//    File/Path object returns its path COMPONENTS, which poisoned the
//    channels with literal 'nfs'/'scratch01' paths).
//  - GFFREAD_CHR still receives its gtf+fasta files as ONE merged list
//    (simplest shape; the job re-pairs by filename).
//
// The mammal/model pipelines (main_mammal.nf, main_model.nf) keep using the
// original single-file modules (agat_gff2gtf.nf, clean_gtf.nf, recode_tga.nf,
// split_if_too_large.nf, secissearch.nf) — do not touch those here.

process AGAT_GFF2GTF_BATCH {
    tag { "agat_gff2gtf_batch_${bucket}" }
    label 'agat'
    memory '6GB'

    input:
    tuple val(bucket), path(gff_files)

    output:
    path "*.gtf", emit: gtf_file

    script:
    """
    #!/bin/bash
    . /usr/local/env-activate.sh

    for gff in *.gff; do
        base="\$(basename "\$gff" .gff)"

        awk -F'\\t' 'BEGIN { OFS="\\t" }

    # Leave comments untouched
    /^#/ {
        print
        next
    }

    {
        original_attrs = \$9

        gene_id = ""
        transcript_id = ""
        id = ""
        parent = ""

        n = split(\$9, attrs, ";")

        for (i = 1; i <= n; i++) {

            # Remove leading/trailing whitespace
            attr = attrs[i]
            gsub(/^[[:space:]]+|[[:space:]]+\$/, "", attr)

            # ----------------------------
            # GTF-style attributes
            # ----------------------------
            if (attr ~ /^gene_id[[:space:]]+/)
                gene_id = attr

            else if (attr ~ /^transcript_id[[:space:]]+/)
                transcript_id = attr

            # ----------------------------
            # GFF3-style attributes
            # ----------------------------
            else if (attr ~ /^ID=/)
                id = attr

            else if (attr ~ /^Parent=/)
                parent = attr
        }

        out = ""

        # ----------------------------------------
        # GTF input
        # Keep gene_id and transcript_id
        # ----------------------------------------
        if (gene_id != "" || transcript_id != "") {

            if (gene_id != "")
                out = gene_id

            if (transcript_id != "")
                out = (out == "" ? transcript_id : out "; " transcript_id)

            out = out ";"
        }

        # ----------------------------------------
        # GFF3 input
        # Keep ID and Parent
        # ----------------------------------------
        else if (id != "" || parent != "") {

            if (id != "")
                out = id

            if (parent != "")
                out = (out == "" ? parent : out ";" parent)
        }

        # ----------------------------------------
        # Unknown/unrecognised format
        # DO NOT destroy attributes
        # ----------------------------------------
        else {
            out = original_attrs
        }

        \$9 = out
        print
    }' "\$gff" > "\${base}.temp.gff"

        agat_convert_sp_gff2gtf.pl \\
                --gff "\${base}.temp.gff" \\
                --gtf_version 3 \\
                --output "\${base}.gtf"

        if [ ! -s "\${base}.gtf" ]; then
            echo "Error: GTF file is empty or not created."
            mv "\${base}.temp.gff" "\${base}.gtf"
        else
            echo "GTF file created successfully."
            rm "\${base}.temp.gff"
        fi
        rm "\$gff"
    done
    """
}

process CLEAN_GTF_BATCH {
    tag { "clean_gtf_batch_${bucket}" }
    memory '4GB'
    cpus 1

    input:
    tuple val(bucket), path(gtf_files)

    output:
    path "*.cleaned.gtf", emit: cleaned_gtf

    script:
    """
    #!/bin/bash
    for gtf in *.gtf; do
        bash clean_gtf.sh "\$gtf"
        rm "\$gtf"
    done
    """
}

process RECODE_TGA_BATCH {
    tag { "recode_tga_batch_${bucket}" }
    label 'python'
    cpus 1
    memory '8GB'

    input:
    tuple val(bucket), path(transcript_fastas)
    val limit

    output:
    path "*_recoded.fa"

    script:
    // getBaseName() (used by the original RECODE_TGA module) strips the file
    // extension, so the recoded output is "<name>_recoded.fa" (no inner .fa).
    // bash basename with the .fa suffix argument matches that exactly.
    """
    #!/bin/bash
    for f in *.fa; do
        base="\$(basename "\$f" .fa)"
        out_name="\${base}_recoded.fa"

        recode_any_TGA.py \\
            --fasta "\$f" \\
            --recodon TGC \\
            --limit ${limit} \\
            --output "\${out_name}"

        rm "\$f"
    done
    """
}

process SPLIT_IF_TOO_LARGE_BATCH {
    tag { "split_if_needed_batch_${bucket}" }
    cpus 1
    memory '4GB'
    label 'splitfasta'

    input:
    tuple val(bucket), path(input_files)

    // optional: a bucket whose input files are ALL empty (transcriptless
    // scaffolds -> 0-byte recoded fasta) makes seqkit split2 emit NOTHING
    // (empty input -> empty .split folder, exit 0, verified), and a required
    // glob would then fail the task with "Missing output file(s)". An empty
    // bucket must simply contribute no part files downstream.
    output:
    path "*.fa.split/*.fa", emit: split_fasta, optional: true

    script:
    """
    #!/bin/bash
    for f in *.fa; do
        seqkit split2 --by-size 40000 "\$f"
    done
    """
}

process SECISSEARCH_BATCH {
    tag { "secissearch_batch_${bucket}" }
    memory '4GB'
    label 'secissearch'

    input:
    tuple val(bucket), path(input_fastas)

    // optional: if EVERY file in a bucket is empty (transcriptless scaffolds),
    // Seblastian is skipped for all of them and NO *.gff is produced; a
    // required glob then fails the task with "Missing output file(s) *.gff"
    // (Giardia_duodenalis job 29112739, 2026-10-02: 23 transcriptless
    // scaffolds, 14 of 24 buckets fully transcriptless). An empty bucket must
    // simply contribute no secis .gff downstream (collectFile skips it).
    output:
    path "*.gff", optional: true

    script:
    // Per-file output prefix computed with the exact same Groovy expression
    // as the original SECISSEARCH module (input_fasta.baseName - '.filtered'),
    // so the emitted *.gff names are identical to the unbatched pipeline.
    // Empty transcript files are SKIPPED, not fed to Seblastian: scaffolds
    // with no transcripts in the GFF (e.g. region/mobile_genetic_element-only
    // scaffolds) yield 0-byte gffread output, and cmsearch dies on them with
    // "Couldn't determine format of sequence file", failing the whole batch
    // job (Trichomonas_vaginalis 28990805, 2026-09-30: 6 such scaffolds).
    def cmds = input_fastas.collect { input_fasta ->
        def fileName = input_fasta.baseName - '.filtered'
        """
        if [ -s ${input_fasta} ]; then
            mkdir -p secis_temp
            python /Seblastian/Seblastian.py \\
                -t ${input_fasta} \\
                -o ${fileName} \\
                -SS \\
                -temp secis_temp \\
                -no_complement \\
                -secis_energy 0.0002
            rm -r secis_temp ${input_fasta}
        else
            echo "WARNING: skipping empty transcript file ${input_fasta} (transcriptless scaffold)" >&2
            rm -f ${input_fasta}
        fi
        """
    }
    """
    #!/bin/bash
    set -euo pipefail
    ${cmds.join('\n')}
    """
}
