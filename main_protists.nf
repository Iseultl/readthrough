#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl = 2

// Input files
params.genome_gtf = params.genome_gtf ?: '/no_backup/rg/ileahy/Mouse_Analysis/gencode.vM37.annotation.gtf'
params.genome_fasta = params.genome_fasta ?: '/no_backup/rg/ileahy/Mouse_Analysis/GRCm39.primary_assembly.genome.fa'
params.lyric_gtf = params.lyric_gtf ?: '/no_backup/rg/ileahy/Mouse_Analysis/lyric_output/lyric_predictions.gtf'
params.species_name = params.species_name ?: 'Mus musculus'
params.output_dir = params.output_dir ?: '/no_backup/rg/ileahy/Mouse_Analysis/secis_independent_output'
params.geneid_param = params.geneid_param ?: '/Users/iseult/Desktop/Geneid_Recoding/testing_false_positives/human3iso.param'
params.help = params.help ?: false
// Number of hash buckets the per-scaffold stages are grouped into
// (modules/scaffold_batches.nf). Each per-scaffold stage runs once per bucket
// and loops over its scaffolds inside the job. For large genomes with few big
// scaffolds most buckets stay empty and behaviour is ~one job per scaffold,
// as before.
params.scaffold_batches = params.scaffold_batches ?: 24

// Print help message if no parameters are provided
def printHelp() {
    log.info """
    SECIS Independent Pipeline
    =============
    
    Usage:
    nextflow run main.nf --genome_gtf <file.gff> --genome_fasta <file.fa> --lyric_gtf <file.gtf> [options]
    
    Mandatory arguments:
      --genome_gtf <file>      Input GFF/GTF file
      --genome_fasta <file>    Input FASTA file
      --lyric_gtf <file>       Input LyRic GTF file

    Optional arguments:
      --output_dir <dir>       Output directory (default: ./results)
      --max_cpus <int>         Maximum number of CPUs (default: 4)
      --max_memory <mem>       Maximum memory (default: 8GB)
      --geneid_param <file>    Path to geneid parameter file
      --debug                  Enable debug mode (default: false)
    """
    exit 0
}

// Print pipeline header
def printHeader() {
    log.info """
    ========================================
    ORFsearch Pipeline - Nextflow Pipeline
    ========================================
    Input GFF:  ${params.genome_gtf}
    Input FASTA: ${params.genome_fasta}
    Output dir:  ${params.output_dir}
    CPU's:       ${params.max_cpus}
    Memory:      ${params.max_memory}
    ========================================
    """
}

// Workflow completion message
def workflowCompletionMessage() {
    log.info """
    ========================================
    Pipeline completed successfully!
    Results are in: ${params.output_dir}
    ========================================
    """
}

// Load modules
include { UNZIP_IF_NEEDED } from './modules/handle_zipped_input'
include { FILTER_ORPHAN_CDS } from './modules/filter_orphan_cds'
include { AGAT_SPLITGFF } from './modules/agat_splitgff'
include { RELOCATE_TRANSCRIPTS } from './modules/relocate_transcripts'
include { SPLITFASTA } from './modules/splitfasta'
// Batched per-scaffold stages (see modules/scaffold_batches.nf): one job per
// hash bucket of scaffolds instead of one job per scaffold, to avoid tens of
// thousands of short SLURM jobs on small-scaffold genomes. The mammal/model
// pipelines keep the original single-file modules.
include { AGAT_GFF2GTF_BATCH } from './modules/scaffold_batches'
include { CLEAN_GTF_BATCH } from './modules/scaffold_batches'
include { RECODE_TGA_BATCH } from './modules/scaffold_batches'
include { SPLIT_IF_TOO_LARGE_BATCH } from './modules/scaffold_batches'
include { SECISSEARCH_BATCH } from './modules/scaffold_batches'
include { GFFREAD_CHR } from './modules/gffread_chr'
include { RUN_GENEID_ORIGINAL } from './modules/run_geneid_original'
include { CONCAT_SUMMARY_RESULTS } from './modules/concat_summary_results'
include { SELECT_INTERESTING } from './modules/select_interesting'
include { GET_ORIGINAL_PREDICTIONS } from './modules/get_original_predictions'
include { CREATE_SUMMARY_TABLE } from './modules/create_summary_table'
include { FILTER_FINAL_TABLE } from './modules/filter_final_table'
include { EXTRACT_SEQUENCE_LOGOS } from './modules/extract_sequence_logos'
include { FILTER_SECIS } from './modules/filter_secis'
include { COMBINE_ORFSECIS } from './modules/combine_orfsecis'
include { CREATE_README } from './modules/create_readme'
include { RUN_SELENOPROFILES } from './modules/run_selenoprofiles'
include { RUN_SECMARKER } from './modules/run_secmarker'

// Helper function to extract chromosome name (without extension)
def get_chr_name(file) {
    return file.getBaseName().replaceFirst(/\.fa$|\.gtf$/, '')
}

// Number of hash buckets the per-scaffold stages are grouped into
// (modules/scaffold_batches.nf). Each per-scaffold stage runs once per bucket
// and loops over its scaffolds inside the job (sequential loop, so peak
// memory is unchanged). For large genomes with few big scaffolds most buckets
// stay empty and behaviour is ~one job per scaffold, as before.
def scaffoldBatchCount() {
    return (params.scaffold_batches ?: 24) as int
}

// Group a channel of per-scaffold files into hash buckets (see above).
// Bucket key = hash of the file base name (per-scaffold stages name every
// file after its scaffold). Stateless, so no ordering assumption is made on
// the lazy channels; the key only needs to partition the channel — stages do
// not need to share buckets. Emits: tuple(val bucket, list of files), one
// item per non-empty bucket.
def batchByScaffold(ch) {
    def K = scaffoldBatchCount()
    return ch
        .map { f -> tuple(((f.name as String).hashCode() % K).abs(), f) }
        .groupTuple()
        .map { bucket, items -> tuple(bucket, items.toList()) }
}

workflow {
    if (params.help) {
        printHelp()
        return
    }
    printHeader()

    // Scaffold batching: the per-scaffold stages below run via the batched
    // modules (modules/scaffold_batches.nf, GFFREAD_CHR) grouped by
    // batchByScaffold() — see the helper definitions above. This avoids tens
    // of thousands of short SLURM jobs on small-scaffold genomes
    // (Cylindrotheca: 2,534 scaffolds -> ~15k jobs before the per-part geneid
    // stages). Every batched process still emits the same per-scaffold output
    // files, so all downstream name/id derivation is unaffected.

    // Step -1: Handle zipped input files
    input_files = Channel.of(
        tuple('genome_fasta', file(params.genome_fasta)),
        tuple('genome_gtf',   file(params.genome_gtf)),
        tuple('lyric_gtf',    file(params.lyric_gtf))
    )

    unzipped_files = UNZIP_IF_NEEDED(input_files)

    // Separate outputs by file type
    genome_fasta_unzipped = unzipped_files
        .filter { type, f -> type == 'genome_fasta' }
        .map { type, f -> f }

    genome_gtf_unzipped = unzipped_files
        .filter { type, f -> type == 'genome_gtf' }
        .map { type, f -> f }

    lyric_gtf_unzipped = unzipped_files
        .filter { type, f -> type == 'lyric_gtf' }
        .map { type, f -> f }

    // Step 0: Create readme
    CREATE_README(params.species_name, genome_fasta_unzipped, genome_gtf_unzipped)

    // Step 0.1: Run Selenoprofiles
    // selenoprofiles_results = RUN_SELENOPROFILES(
    //    genome_fasta_unzipped,
    //    genome_gtf_unzipped,
    //    params.species_name
    // )

    // Step 0.2: Run Secmarker
    secmarker_results = RUN_SECMARKER(genome_fasta_unzipped)

    // Step 0.3: Remove orphan CDS records from the LyRic GFF (transdecoder
    // ORF scans with .pN parents and no transcript line; gffread would
    // otherwise materialize them as spurious CDS-only transcripts)
    lyric_gtf_clean = FILTER_ORPHAN_CDS(lyric_gtf_unzipped)

    // Step 1: Split GTF/GFF
    split_results = AGAT_SPLITGFF(lyric_gtf_clean)

    // Step 2: Collect all .gff files from the output directory
    gff_files_ch = split_results.gff_files
    gff_files_ch = gff_files_ch.flatten()

    // Step 3: Run AGAT_GFF2GTF in parallel to standardise each GFF file
    // (batched: scaffold_batches jobs, each looping over its bucket; the glob
    // emit carries one List per job, flatten back to per-scaffold files)
    gtf_files_ch = AGAT_GFF2GTF_BATCH(batchByScaffold(gff_files_ch)).gtf_file.flatten()

    // Create a channel of cleaned GTF files for downstream processing
    // (batched, flattened as above)
    split_gff_dir_ch = CLEAN_GTF_BATCH(batchByScaffold(gtf_files_ch)).cleaned_gtf.flatten()
    
    // Create paired channels for GTF and FASTA files
    // First, create a channel for GTF files with chromosome names
    base_names_gtf = gtf_files_ch.map { file ->
        def chr = file.name.replaceFirst(/\.gtf$/, '')
        tuple(chr, file)
    }
    
    // Step 4: FASTA Processing Pipeline
    split_fasta_dir_ch = SPLITFASTA(genome_fasta_unzipped)
    base_names_fasta = split_fasta_dir_ch.flatten().map { file ->
        def fname = file.name  // e.g. horse_genome.part_NW_027222397.1.fa
        def chr_match = fname =~ /part_(.+)\.fa/
        def chr = chr_match ? chr_match[0][1] : null
        tuple(chr, file)
    }
    
    // Step 5: Create paired gtf & fasta channel
    paired_ch = base_names_gtf.combine(base_names_fasta, by: 0)  
     
    // Step 6: Create relocated to transcript gff 
    concatenated_gtf = split_gff_dir_ch.collectFile(name: 'concatenated.gtf')
    relocated_gtf = RELOCATE_TRANSCRIPTS(concatenated_gtf).relocated_gtf
    relocated_gtf_val = relocated_gtf.first()

    // Step 7: GFFREAD — pass batches of paired (gtf, fasta) scaffolds.
    // The pair for each scaffold stays in the same bucket (grouped AFTER the
    // combine, so both files of a matched pair travel together).
    // The bucket's gtf and fasta files are merged into ONE list input —
    // two list inputs in one process collide at staging time (see
    // modules/gffread_chr.nf); the job re-pairs them by filename.
    // Note: 1-parameter map closures with explicit list indexing — this Nextflow
    // version passes a combine's 3-element tuple item to the map closure as a
    // single list argument, which multi-parameter closures cannot be dispatched
    // against (MissingMethodException at MapOp).
    gffread_batches = paired_ch
        .map { item ->
            tuple(((item[0] as String).hashCode() % scaffoldBatchCount()).abs(), item[1], item[2])
        }
        .groupTuple()
        .map { item ->
            // NOTE: groupTuple() is COLUMN-WISE — for 3-tuples it emits
            // (key, List<value1>, List<value2>), NOT (key, List<(v1,v2)>).
            // So item[1] is the bucket's gtf files, item[2] the fasta files.
            def gtfs = item[1].toList()
            def fastas = item[2].toList()
            // GFFREAD_CHR takes ONE list input (two list inputs in one process
            // are avoided; see modules/gffread_chr.nf); the job re-pairs
            // gtf<->fasta by filename.
            tuple(item[0], gtfs + fastas)
        }
    gffread_out = GFFREAD_CHR(gffread_batches)
    // Glob emit -> one List of per-scaffold files per batch job; flatten back
    // to a per-scaffold channel (consumed by recode, secissearch and logos).
    transcripts_ch = gffread_out.transcripts.flatten()

    // Step 10. Recode all transcripts (batched)
    recoded_transcripts = RECODE_TGA_BATCH(batchByScaffold(transcripts_ch), 100000).flatten()

    // Step 11. Split recoded transcripts if too large (batched)
    split_transcripts_ch = SPLIT_IF_TOO_LARGE_BATCH(batchByScaffold(recoded_transcripts)).split_fasta.flatten()
     
    // Step 12: Run geneid on all the transcript sequences 
    geneid_results_ch = RUN_GENEID_ORIGINAL(split_transcripts_ch, params.geneid_param).geneid_results 
    
    // Step 14 & 15: Process original and recoded predictions
    interesting_predictions = SELECT_INTERESTING(geneid_results_ch).interesting_predictions
    
    original_predictions = GET_ORIGINAL_PREDICTIONS(geneid_results_ch).original_predictions
    
    // Combine the channels based on the scaffold ID
    combined_predictions = interesting_predictions.combine(original_predictions, by: 0)
    
    // Step 16: Final Output - Pass the combined channel to the summary table process
    summary_tables = CREATE_SUMMARY_TABLE(combined_predictions, relocated_gtf_val)
    
    result = summary_tables.summary_table.collectFile(name: 'ORFsearch.filter')
    result.view { item -> "Concatenated summary results: ${item}" }

    // Step 17: Filter final table to handle the duplicates from split_if_too_large
    ORFsearch_result = FILTER_FINAL_TABLE(result)
    ORFsearch_result.view { item -> "Filtered ORFsearch results: ${item}" }
    // Step 18: Extract sequences for logos
    // gffread transcripts are emitted per scaffold; collect() makes them a single
    // list input so EXTRACT_SEQUENCE_LOGOS runs ONCE with all sequences.
    // (Passing the per-scaffold channel directly would run the process once per
    // scaffold and every instance overwrites the same four output files;
    // collectFile() does not work on dirs: EISDIR.)
    extracted_sequences = EXTRACT_SEQUENCE_LOGOS(ORFsearch_result, transcripts_ch.collect(), params.species_name)

    // Step 19: Run SECISearch on the transcripts (batched)
    secissearch_results = SECISSEARCH_BATCH(batchByScaffold(transcripts_ch)).flatten()
    merged_secis_gff = secissearch_results.collectFile(name: 'all_secis_combined.gff')

    // Step 20: Filter SECISearch results
    filtered_secis = FILTER_SECIS(merged_secis_gff)

    // Step 21: Combine ORFsearch results with filtered SECISearch results
    ORFsearch_result = COMBINE_ORFSECIS(ORFsearch_result, merged_secis_gff, filtered_secis, params.species_name)
    ORFsearch_result.view { item -> "Combined ORFsearch and SECISearch results: ${item}" }
    workflowCompletionMessage()
}
