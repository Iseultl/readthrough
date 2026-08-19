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
include { AGAT_GFF2GTF } from './modules/agat_gff2gtf'
include { AGAT_SPLITGFF } from './modules/agat_splitgff'
include { CLEAN_GTF } from './modules/clean_gtf'
include { RELOCATE_TRANSCRIPTS } from './modules/relocate_transcripts'
include { SPLITFASTA } from './modules/splitfasta'
include { GFFREAD_CHR } from './modules/gffread_chr'
include { RECODE_TGA } from './modules/recode_tga'
include { SPLIT_IF_TOO_LARGE } from './modules/split_if_too_large'
include { RUN_GENEID_ORIGINAL } from './modules/run_geneid_original'
include { CONCAT_SUMMARY_RESULTS } from './modules/concat_summary_results'
include { SELECT_INTERESTING } from './modules/select_interesting'
include { GET_ORIGINAL_PREDICTIONS } from './modules/get_original_predictions'
include { CREATE_SUMMARY_TABLE } from './modules/create_summary_table'
include { FILTER_FINAL_TABLE } from './modules/filter_final_table'
include { EXTRACT_SEQUENCE_LOGOS } from './modules/extract_sequence_logos'
include { SECISSEARCH } from './modules/secissearch'
include { FILTER_SECIS } from './modules/filter_secis'
include { COMBINE_ORFSECIS } from './modules/combine_orfsecis'
include { CREATE_README } from './modules/create_readme'
include { RUN_SELENOPROFILES } from './modules/run_selenoprofiles'
include { RUN_SECMARKER } from './modules/run_secmarker'

// Helper function to extract chromosome name (without extension)
def get_chr_name(file) {
    return file.getBaseName().replaceFirst(/\.fa$|\.gtf$/, '')
}

workflow {
    if (params.help) {
        printHelp()
        return
    }
    printHeader()
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
    selenoprofiles_results = RUN_SELENOPROFILES(
        genome_fasta_unzipped,
        genome_gtf_unzipped,
        params.species_name
    )

    // Step 0.2: Run Secmarker
    secmarker_results = RUN_SECMARKER(genome_fasta_unzipped)

    // Step 1: Split GTF/GFF
    split_results = AGAT_SPLITGFF(lyric_gtf_unzipped)

    // Step 2: Collect all .gff files from the output directory
    gff_files_ch = split_results.gff_files
    gff_files_ch = gff_files_ch.flatten()

    // Step 3: Run AGAT_GFF2GTF in parallel to standardise each GFF file
    gtf_files_ch = AGAT_GFF2GTF(gff_files_ch)
    
    // Create a channel of cleaned GTF files for downstream processing
    split_gff_dir_ch = CLEAN_GTF(gtf_files_ch)
    
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
    relocated_gtf = RELOCATE_TRANSCRIPTS(concatenated_gtf)
    
    // Step 7: Now pass the paired channel to GFFREAD
    gffread_out = GFFREAD_CHR(paired_ch)
     
    // Step 10. Recode all transcripts 
    recoded_transcripts = RECODE_TGA(gffread_out.transcripts)

    // Step 11. Split recoded transcripts if too large
    split_transcripts_ch = SPLIT_IF_TOO_LARGE(recoded_transcripts).split_fasta.flatten()
     
    // Step 12: Run geneid on all the transcript sequences 
    geneid_results_ch = RUN_GENEID_ORIGINAL(split_transcripts_ch, params.geneid_param).out.geneid_results 
    geneid_results_ch.view { item -> "Geneid results: ${item}"}
    // Step 14 & 15: Process original and recoded predictions
    interesting_predictions = SELECT_INTERESTING(geneid_results_ch, relocated_gtf).interesting_predictions
    interesting_predictions.view { item -> "Interesting predictions: ${item}" }
    original_predictions = GET_ORIGINAL_PREDICTIONS(geneid_results_ch).original_predictions
    original_predictions.view { item -> "Original predictions: ${item}" }
    // Combine the channels based on the scaffold ID
    combined_predictions = interesting_predictions.combine(original_predictions, by: 0)
    combined_predictions.view { item -> "Combined predictions: ${item}" }
    // Step 16: Final Output - Pass the combined channel to the summary table process
    summary_tables = CREATE_SUMMARY_TABLE(combined_predictions, relocated_gtf)
    summary_tables.view { item -> "Summary tables: ${item}" }
    result = CONCAT_SUMMARY_RESULTS(summary_tables.collect())
    result.view { item -> "Concatenated summary results: ${item}" }
    // Step 17: Filter final table to handle the duplicates from split_if_too_large
    ORFsearch_result = FILTER_FINAL_TABLE(result)
    ORFsearch_result.view { item -> "Filtered ORFsearch results: ${item}" }
    // Step 18: Extract sequences for logos
    extracted_sequences = EXTRACT_SEQUENCE_LOGOS(ORFsearch_result, gffread_out.gffread_dir, params.species_name)

    // Step 19: Run SECISearch on the transcripts
    secissearch_results = SECISSEARCH(gffread_out.transcripts)
    merged_secis_gff = secissearch_results.collectFile(name: 'all_secis_combined.gff')

    // Step 20: Filter SECISearch results
    filtered_secis = FILTER_SECIS(merged_secis_gff)

    // Step 21: Combine ORFsearch results with filtered SECISearch results
    ORFsearch_result = COMBINE_ORFSECIS(ORFsearch_result, merged_secis_gff, filtered_secis, params.species_name)
    ORFsearch_result.view { item -> "Combined ORFsearch and SECISearch results: ${item}" }
    workflowCompletionMessage()
}
