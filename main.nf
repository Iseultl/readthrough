#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl = 2

// Input files
params.genome_gtf = params.genome_gtf ?: '/no_backup/rg/ileahy/Mouse_Analysis/gencode.vM37.annotation.gtf'
params.genome_fasta = params.genome_fasta ?: '/no_backup/rg/ileahy/Mouse_Analysis/GRCm39.primary_assembly.genome.fa'
params.transcripts_clean = params.transcripts_clean // Optional: provide a transcripts_clean file
params.output_dir = params.output_dir ?: '/no_backup/rg/ileahy/Mouse_Analysis/secis_independent_output'
params.geneid_param = params.geneid_param ?: '/Users/iseult/Desktop/Geneid_Recoding/testing_false_positives/human3iso.param'

// Print help message if no parameters are provided
if (params.help) {
    log.info """
    SECIS Independent Pipeline
    =============
    
    Usage:
    nextflow run main.nf --genome_gtf <file.gff> --genome_fasta <file.fa> [options]
    
    Mandatory arguments:
      --genome_gtf <file>      Input GFF/GTF file
      --genome_fasta <file>    Input FASTA file
      
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

// Load modules
include { AGAT_GFF2GTF } from './modules/agat_gff2gtf'
include { AGAT_SPLITGFF } from './modules/agat_splitgff'
include { CLEAN_GTF } from './modules/clean_gtf'
include { CONCATENATE_GTFS } from './modules/concatenated_gff'
include { RELOCATE_TRANSCRIPTS } from './modules/relocate_transcripts'
include { SPLITFASTA } from './modules/splitfasta'
include { GFFREAD } from './modules/gffread'
include { RECODE_TGA } from './modules/recode_tga'
include { SPLIT_RECODED } from './modules/split_recoded'
include { SPLIT_IF_TOO_LARGE } from './modules/split_if_too_large'
include { RUN_GENEID_ORIGINAL } from './modules/run_geneid_original'
include { CONCAT_SUMMARY_RESULTS } from './modules/concat_summary_results'
include { SELECT_INTERESTING } from './modules/select_interesting'
include { GET_ORIGINAL_PREDICTIONS } from './modules/get_original_predictions'
include { CREATE_SUMMARY_TABLE } from './modules/create_summary_table'
include { FILTER_FINAL_TABLE } from './modules/filter_final_table'

// Helper function to extract chromosome name (without extension)
def get_chr_name(file) {
    return file.getBaseName().replaceFirst(/\.fa$|\.gtf$/, '')
}

workflow {
    // Step 1: Split GFF by chromosome
    split_results = AGAT_SPLITGFF(params.genome_gtf)

    // Step 2: Collect all .gff files from the output directory
    gff_files_ch = split_results.gff_files
    gff_files_ch = gff_files_ch.flatten()

    // Step 3: Run AGAT_GFF2GTF in parallel to standardise each GFF file
    gtf_files_ch = AGAT_GFF2GTF(gff_files_ch)
    
    split_gff_dir_ch = CLEAN_GTF(gtf_files_ch)
     
    base_names_gff = split_gff_dir_ch.map { file ->
        def chr = file.name.replaceFirst(/\.cleaned\.gtf$/, '')
        tuple(chr, file)
    }
    
    // Step 4: FASTA Processing Pipeline
    split_fasta_dir_ch = SPLITFASTA(params.genome_fasta)
    base_names_fasta = split_fasta_dir_ch.flatten().map { file ->
        def fname = file.name  // e.g. horse_genome.part_NW_027222397.1.fa
        def chr_match = fname =~ /part_(.+)\.fa/
        def chr = chr_match ? chr_match[0][1] : null
        tuple(chr, file)
    }
    
    // Step 5: Create paired gtf & fasta channel
    paired_ch = base_names_gff.combine(base_names_fasta, by: 0)  
     
    // Step 6: Create relocated to transcript gff 
    CONCATENATE_GTFS(split_gff_dir_ch.collect())
    RELOCATE_TRANSCRIPTS(CONCATENATE_GTFS.out)
    
    // Step 7: Now pass the paired channel to GFFREAD
    gffread_out = GFFREAD(paired_ch)
     
    // Step 10. Recode all transcripts 
    recoded_transcripts = RECODE_TGA(GFFREAD.out)

    // Step 11. Split recoded transcripts if too large
    split_transcripts_ch = SPLIT_IF_TOO_LARGE(recoded_transcripts).split_fasta.flatten()
     
    // Step 12: Run geneid on all the transcript sequences 
    geneid_results_ch = RUN_GENEID_ORIGINAL(split_transcripts_ch, params.geneid_param)  
 
    // Step 14 & 15: Process original and recoded predictions
    interesting_predictions = SELECT_INTERESTING(geneid_results_ch, RELOCATE_TRANSCRIPTS.out).interesting_predictions
    original_predictions = GET_ORIGINAL_PREDICTIONS(geneid_results_ch).original_predictions

    // Combine the channels based on the scaffold ID
    combined_predictions = interesting_predictions.combine(original_predictions, by: 0)
    combined_predictions.view()
    return
    // Step 16: Final Output - Pass the combined channel to the summary table process
    summary_tables = CREATE_SUMMARY_TABLE(combined_predictions, RELOCATE_TRANSCRIPTS.out)

    result = CONCAT_SUMMARY_RESULTS(summary_tables)
    
    // Step 17: Filter final table to handle the duplicates from split_if_too_large
    ORFsearch_result = FILTER_FINAL_TABLE(result)
}

// Workflow completion message
workflow.onComplete {
    log.info """
    ========================================
    Pipeline completed successfully!
    Results are in: ${params.output_dir}
    ========================================
    """
}
