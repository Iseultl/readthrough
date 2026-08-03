# Project Overview
The purpose of the nextflow.yml is to download all the available mammalian genomes and annotations file and test the open reading frames ending in TGA codons for readthrough capability. With these candidates we will be able to identify the candidate genes. We will also extract a 600 nucleotide window around the TGA site which will allow us to align our hits and observe if the protein sequence is maintained after the TGA site for all species. 

# Nextflow workflow guidelines

## Verification & Environmental guardrails
- This pipeline requires environmental set up that is detailed in the build-tools.yml using the Dockerfile in the docker/ folder. 
- Check that environmental set up does not take all of the available memory in each VM as this will prevent the rest of the pipeline from running as it also requires a lot of memory. 
- The Dockerfile for environmental setup is new - please ensure it is inline with all of the nextflow steps - in particular the processes which call on python files as we should not be calling from the Nextflow config file anymore. 
- Check that all the inputs and outputs for each step of the nextflow pipeline transfer their inputs and outputs to each other smoothly. 

## Key Development Commands 
- The first step is to build the matrix which reads the mammals_to_search.tsv file and builds a matrix for input into the nextflow pipeline. 
- The next step is the nextflow pipeline which uses the matrix as input, and runs the main_mammal.nf script usig the species taxid, species name, and geneid parameter file. 
- The output of the nextflow pipeline is the *SECIS.result file and the folder containing the sequence_logos. This will allow us to identify the genes with positive recoding scores and also align the genes across species using the outputted sequence files. 

## Verification requirements
- Ensure that the environmental setups are not taking up all the available memory. 
- Workflow is complete if the output *SECIS.result file is extracted and the seqence logos are extracted with sequences present in the output files. 
