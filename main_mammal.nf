#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl = 2

// Input value
params.species_taxid = params.species_taxid ?: '10090'
params.species_name = params.species_name ?: 'Mus musculus'
params.output_dir = params.output_dir ?: '/results'
params.geneid_param = params.geneid_param ?: 'human3iso.param'
params.help = params.help ?: false

// Print help message if no parameters are provided
if (params.help) {
    log.info """
    Selenoprotein Hunter Pipeline
    =============
    
    Usage:
    nextflow run main.nf --taxid <taxid> --species <name> [options]
    
    Mandatory arguments:
      --taxid <taxid>          Taxonomic ID of the species
      --species <name>         Species name (for labeling) 

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
 Selenoprotein Hunter Pipeline - Nextflow Pipeline
========================================
Input TaxID:  ${params.species_taxid}
Input Species:  ${params.species_name}
Output dir:  ${params.output_dir}
CPU's:       ${params.max_cpus}
Memory:      ${params.max_memory}
========================================
"""

// Load modules
include { DOWNLOAD_TAXID } from './modules/download_taxid'
include { UNZIP_IF_NEEDED } from './modules/handle_zipped_input'
include { AGAT_GFF2GTF } from './modules/agat_gff2gtf'
include { AGAT_SPLITGFF } from './modules/agat_splitgff'
include { CLEAN_GTF } from './modules/clean_gtf'
include { CONCATENATE_GTFS } from './modules/concatenated_gff'
include { RELOCATE_TRANSCRIPTS } from './modules/relocate_transcripts'
include { SPLITFASTA } from './modules/splitfasta'
include { GFFREAD } from './modules/gffread'
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
include { MERGE_GFF } from './modules/merge_gff'
include { FILTER_SECIS } from './modules/filter_secis'
include { COMBINE_ORFSECIS } from './modules/combine_orfsecis'
include { CREATE_README } from './modules/create_readme'
include { SELENOPROFILES_ANALYSIS } from './modules/secorfsearch_profiles'
include { SELENOPROFILES_GFFREAD } from './modules/selenoprofiles_gffread'
include { RUN_SELENOPROFILES } from './modules/run_selenoprofiles'


// Helper function to extract chromosome name (without extension)
def get_chr_name(file) {
    return file.getBaseName().replaceFirst(/\.fa$|\.gtf$/, '')
}

workflow {
    // Step 0: Download genome files based on taxid
    downloaded_files = DOWNLOAD_TAXID(params.species_taxid)
    
    // Print the downloaded files for debugging
    downloaded_files.map { fasta, gff ->
        println "$fasta $gff"
    }
    
    input_files = downloaded_files.flatMap { fasta_file, gff_file ->
        [
            tuple('genome_${params.species_taxid}.fasta.gz', fasta_file),
            tuple('genome_${params.species_taxid}.gff.gz', gff_file)
        ]
    }
    
    unzipped_files = UNZIP_IF_NEEDED(input_files)

    // Separate outputs by file type
    genome_fasta_unzipped = unzipped_files
        .filter { type, file -> type.endsWith(".fasta.gz") }
        .map { type, file -> file }
    genome_fasta_unzipped
    .multiMap { fasta ->
        readme: fasta
        split: fasta
        selenoprofiles: fasta
        selenoprofiles_base: fasta
    }
    .set { fasta }
    
    genome_gff_unzipped = unzipped_files
        .filter { type, file -> type.endsWith(".gff.gz") }
        .map { type, file -> file }
    genome_gff_unzipped
        .multiMap { gff ->
            readme: gff
            split: gff
            selenoprofiles: gff
            selenoprofiles_base: gff
        }
        .set { gff }
    
    // Step 1: Create readme
    readme = CREATE_README(params.species_name, fasta.readme, gff.readme)
    
    genome_for_selenoprofiles =
        fasta.selenoprofiles
            .combine(gff.selenoprofiles)
            .map { fasta, gff ->
                tuple(params.species_taxid, fasta, gff)
            }
    
    // Step 3: Split GTF/GFF
    split_results = AGAT_SPLITGFF(gff.split)
    
    // Step 4: Collect all .gff files from the output directory
    gff_files_ch = split_results.gff_files
    gff_files_ch = gff_files_ch.flatten()

    // Step 5: Run AGAT_GFF2GTF in parallel to standardise each GFF file 
    gtf_files_ch = AGAT_GFF2GTF(gff_files_ch)
    
    // Create a channel of cleaned GTF files for downstream processing
    split_gff_dir_ch = CLEAN_GTF(gtf_files_ch)
    
    // Create paired channels for GTF and FASTA files
    // First, create a channel for GTF files with chromosome names
    base_names_gtf = gtf_files_ch.map { file ->
        def chr = file.name.replaceFirst(/\.gtf$/, '')
        tuple(chr, file)
    }

    // Step 6: FASTA Processing Pipeline
    split_fasta_dir_ch = SPLITFASTA(fasta.split).split_chr
    base_names_fasta = split_fasta_dir_ch.flatten().map { file ->
        def fname = file.name  // e.g. horse_genome.part_NW_027222397.1.fa
        def chr_match = fname =~ /part_(.+)\.fa/
        def chr = chr_match ? chr_match[0][1] : null
        tuple(chr, file)
    }
    
    // Step 7: Create paired gtf & fasta channel
    paired_ch = base_names_gtf.combine(base_names_fasta, by: 0)  
    
    // Step 8: Now pass the paired channel to GFFREAD
    gffread_out = GFFREAD(paired_ch)

    // Step 9: Create relocated to transcript gff 
    CONCATENATE_GTFS(split_gff_dir_ch.collect(), gffread_out.collect())
    relocated_gtf = RELOCATE_TRANSCRIPTS(CONCATENATE_GTFS.out)
     
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
    
    // Step 16: Final Output - Pass the combined channel to the summary table process
    summary_tables = CREATE_SUMMARY_TABLE(combined_predictions, RELOCATE_TRANSCRIPTS.out)

    result = CONCAT_SUMMARY_RESULTS(summary_tables.collect())
    
    // Step 17: Filter final table to handle the duplicates from split_if_too_large
    ORFsearch_result = FILTER_FINAL_TABLE(result)

    // Step 18: Extract sequences for logos
    extracted_sequences = EXTRACT_SEQUENCE_LOGOS(ORFsearch_result, gffread_out.collect())

    // Step 19: Run SECISearch on the transcripts
    secissearch_results = SECISSEARCH(split_transcripts_ch)
    secissearch_results.collect().set { all_gffs }
    merged_secis_gff = MERGE_GFF(all_gffs)

    // Step 20: Filter SECISearch results
    filtered_secis = FILTER_SECIS(merged_secis_gff)

    // Step 21: Combine ORFsearch results with filtered SECISearch results
    ORFsearch_result = COMBINE_ORFSECIS(ORFsearch_result, merged_secis_gff, filtered_secis)

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
