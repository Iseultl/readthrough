#!/usr/bin/env python3

# The purpose of this script is to relocate he gtf file coordinates to the 
# position along the transcript rather than their genomic positions

import argparse
import pandas as pd
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time

# Load GFF file into dataframe
def load_gff(gff_file):
    df = pd.read_csv(gff_file, sep='\t', comment='#', header=None)
    df.columns = ['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attributes']
    df['Start'] -= 1  # Convert to 0-based indexing
    
    # Extract gene names from Attributes column
    extracted_gene_id = df['Attributes'].str.extract(r'gene_id "(.*?)"')[0]
    extracted_gene_id = extracted_gene_id.str.replace("gene-", "", regex=False)

    extracted_transcript_id = df['Attributes'].str.extract(r'transcript_id "(.*?)"')[0]
    extracted_transcript_id = extracted_transcript_id.str.replace("rna-", "", regex=False)

    df['Attributes'] = extracted_gene_id + "_" + extracted_transcript_id
    print(df.head())
    df = df[df['Attributes'] != '.'].copy()
    print(df.head()) 
    return df

def group_attributes(gff):
    sorted_gff = gff.sort_values(by=["Attributes", "Start"]).copy()
    grouped_by_attr = dict(tuple(sorted_gff.groupby("Attributes")))
    return grouped_by_attr 

# First get the new total transcript lengths
def relocate_trancripts(grouped):
    relocated_transcripts = []
    
    for attr, transcript_group in grouped.items():
        transcript_dct = {}
        
        # Extract rows
        transcript_row = transcript_group[transcript_group["Feature"] == "transcript"]
        exon_rows = transcript_group[transcript_group["Feature"] == "exon"].copy()
        strand = transcript_row.iloc[0]["Strand"]
        chr = transcript_row.iloc[0]['Chromosome']
        
        total_transcript_length = 0
        
        for _, row in exon_rows.iterrows():
            exon_length = row['End'] - row['Start']
            total_transcript_length = total_transcript_length + exon_length 
            
        transcript_dct = {
                "Chromosome": chr,
                "Source": '.',
                "Feature": "transcript",
                "Start": 0,
                "End": total_transcript_length,
                "Score": '.',
                "Strand": strand,
                "Frame": '.',
                "Attributes": transcript_row.iloc[0]['Attributes'] 
            }
        relocated_transcripts.append(transcript_dct) 
    return relocated_transcripts

# Second get the CDS locations 
def relocate_CDS(grouped):

    relocated_CDS = []
    relocated_exons = []
    for _, transcript_group in grouped.items():
        exon_dct = {}
        CDS_dct = {}
 
        # Extract rows
        exon_rows = transcript_group[transcript_group["Feature"] == "exon"]
        CDS_rows = transcript_group[transcript_group["Feature"] == "CDS"]
        transcript_row = transcript_group[transcript_group["Feature"] == "transcript"]
        if exon_rows.empty or CDS_rows.empty or transcript_row.empty:
            continue 
        else:
            transcript_start = transcript_row.iloc[0]['Start']
            transcript_end = transcript_row.iloc[0]['End']
            strand = exon_rows.iloc[0]["Strand"]
            chr = exon_rows.iloc[0]['Chromosome']
            
            # First exon starts at 0
            exon_position = 0
            
            for i, (_, row) in enumerate(exon_rows.iterrows()):
                exon_len = row['End'] - row['Start']
                exon_dct = {
                    "Chromosome": chr,
                    "Source": '.',
                    "Feature": "exon",
                    "Start": exon_position + 1,
                    "End": exon_position + exon_len,
                    "Score": '.',
                    "Strand": strand,
                    "Frame": '.',
                    "Attributes": row['Attributes'] 
                }
                
                
                exon_position = exon_position + exon_len
                relocated_exons.append(exon_dct)
            
            # First CDS starts at end of UTR
            # Except in cases where the first exon and CDS start at the same position
            # Need to handle cases of more than one exon before first CDS 
            UTR_5 = 0
            
            for _, row in exon_rows.iterrows():
                # Capture the exon(s) that occur before the CDS
                if row['End'] < CDS_rows['Start'].min():
                    exon_before_CDS = row['End'] - row['Start']
                    UTR_5 = UTR_5 + exon_before_CDS
                if (row['Start'] < CDS_rows['Start'].min()) and (CDS_rows['End'].min() == row['End']):
                    len_before_CDS = CDS_rows['Start'].min() - row['Start']
                    UTR_5 = UTR_5 + len_before_CDS
                if (row['Start'] < CDS_rows['Start'].min()) and (CDS_rows['End'].min() < row['End']):
                    len_before_CDS = CDS_rows['Start'].min() - row['Start']
                    UTR_5 = UTR_5 + len_before_CDS 

            CDS_position = UTR_5

            for _, row in CDS_rows.iterrows():
                CDS_len = row['End'] - row['Start']
                CDS_dct = {
                    "Chromosome": chr,
                    "Source": '.',
                    "Feature": "CDS",
                    "Start": CDS_position + 1,
                    "End": CDS_position + CDS_len,
                    "Score": '.',
                    "Strand": strand,
                    "Frame": '.',
                    "Attributes": row['Attributes'] 
                }
                CDS_position = CDS_position + CDS_len
                relocated_CDS.append(CDS_dct)
            
    return relocated_CDS, relocated_exons
        
# Relocate the Sec
def relocate_sec(grouped):
    
    relocated_Sec = []
    for _, transcript_group in grouped.items():
        sec_dct = {}
        
        # Extract rows
        exon_rows = transcript_group[transcript_group["Feature"] == "exon"]
        sec_rows = transcript_group[transcript_group["Feature"] == "Selenocysteine"]
        if exon_rows.empty or sec_rows.empty:
            continue 
        else:
            strand = exon_rows.iloc[0]["Strand"]
            chr = exon_rows.iloc[0]['Chromosome']
            
            sec_pos = 0
            for _, row in exon_rows.iterrows():
                if (row['Start'] < sec_rows['Start'].min()) and (sec_rows['Start'].min() < row['End']):
                    exon_pos = sec_rows['Start'].min() - row['Start']
                    break
                else:
                    exon_len = row['End'] - row['Start']
                    sec_pos = sec_pos + exon_len

            sec_pos = int(sec_pos)+int(exon_pos)              

            Sec_dct = {
                    "Chromosome": chr,
                    "Source": '.',
                    "Feature": "Selenocysteine",
                    "Start": sec_pos + 1,
                    "End": sec_pos + 3,
                    "Score": '.',
                    "Strand": strand,
                    "Frame": '.',
                    "Attributes": exon_rows.iloc[0]['Attributes']  
                }
            relocated_Sec.append(Sec_dct)
    
    return relocated_Sec
       
# Relocate the secis
def relocate_secis(grouped):
    relocated_secis = []
    
    for _, transcript_group in grouped.items():
        transcript_dct = {}
        
        # Extract rows
        transcript_row = transcript_group[transcript_group["Feature"] == "transcript"]
        exon_rows = transcript_group[transcript_group["Feature"] == "exon"]
        secis_rows = transcript_group[transcript_group["Feature"] == "secis"] 
        if transcript_row.empty or exon_rows.empty or secis_rows.empty:
            continue
        else:
            strand = transcript_row.iloc[0]["Strand"]
            chr = transcript_row.iloc[0]['Chromosome']
            
            total_transcript_length = 0
            for _, row in exon_rows.iterrows():
                exon_length = row['End'] - row['Start']
                total_transcript_length = total_transcript_length + exon_length
            
            for _, row in secis_rows.iterrows():
                if strand == "+":
                    start_dist = transcript_row['End'].max() - row['Start']
                    end_dist = transcript_row['End'].max() - row['End']
                    new_start = total_transcript_length - start_dist
                    new_end = total_transcript_length - end_dist 
                    secis_dct = {
                            "Chromosome": chr,
                            "Source": '.',
                            "Feature": "secis",
                            "Start": new_start,
                            "End": new_end,
                            "Score": '.',
                            "Strand": strand,
                            "Frame": '.',
                            "Attributes": row['Attributes']  
                        }
                else:
                    start_dist = row['Start'] - transcript_row['Start'].min()
                    end_dist = row['End'] - transcript_row['Start'].min()
                    new_start = total_transcript_length - end_dist
                    new_end = total_transcript_length - start_dist
                    
                    secis_dct = {
                            "Chromosome": chr,
                            "Source": '.',
                            "Feature": "secis",
                            "Start": start_dist,
                            "End": end_dist,
                            "Score": '.',
                            "Strand": "-",
                            "Frame": '.',
                            "Attributes": row['Attributes']  
                        } 
                relocated_secis.append(secis_dct)
    
    return relocated_secis



def handle_negs(df):
    new_df = df.copy()

    # Precompute a mask for all negative strand transcripts
    neg_transcripts = df[(df['Strand'] == "-") & (df['Feature'] == 'transcript')]

    # Map Attributes → End for all negative strand transcripts
    transcript_end_map = dict(zip(neg_transcripts['Attributes'], neg_transcripts['End']))

    # Mask for all rows on the negative strand that need adjusting
    features_to_adjust = ['secis', 'Selenocysteine', 'CDS', 'exon']
    mask = (new_df['Strand'] == "-") & (new_df['Feature'].isin(features_to_adjust))

    # Filter only relevant rows for speed
    affected = new_df[mask].copy()

    # Map each row to its transcript end
    affected['transcript_len'] = affected['Attributes'].map(transcript_end_map)

    # Drop rows without transcript info (e.g., missing map)
    affected = affected.dropna(subset=['transcript_len'])

    # Adjust Start and End positions
    start_old = affected['Start'].copy()
    end_old = affected['End'].copy()

    affected['Start'] = affected['transcript_len'] - end_old
    affected['End'] = affected['transcript_len'] - start_old

    # Update original dataframe
    new_df.update(affected[['Start', 'End']])

    # Set strand to '+' for the affected transcripts and features
    to_flip = (new_df['Strand'] == "-") & (
        new_df['Feature'].isin(features_to_adjust + ['transcript'])
    ) & new_df['Attributes'].isin(transcript_end_map.keys())

    new_df.loc[to_flip, 'Strand'] = "+"

    return new_df
                 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Recode TGA codons at selenocysteine positions in a genome FASTA file")
    parser.add_argument('--gff', type=str, required=True, help="Path to GFF file seleno_secis_profiles_transcripts.gff")
    parser.add_argument('--output', type=str, required=True, help="Path to output gtf file - positioned on the spliced transcript")
    args = parser.parse_args()
    
    start = time.time()
    gff = load_gff(args.gff)
    groups = group_attributes(gff)
    transcript_lst = relocate_trancripts(groups)
    print("Relocate transcripts took", time.time() - start, "seconds")
    relocated_CDS, relocated_exons = relocate_CDS(groups) 
    print("Relocate CDS took", time.time() - start, "seconds")
    relocated_sec = relocate_sec(groups)
    print("Relocate sec took", time.time() - start, "seconds")
    relocated_secis = relocate_secis(groups) 
    print("Relocate secis took", time.time() - start, "seconds")

    columns = ['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attributes']
    df = pd.DataFrame.from_dict(transcript_lst+relocated_CDS+relocated_exons+relocated_sec+relocated_secis)
 
    df = handle_negs(df)
    print("Relocate negatives took", time.time() - start, "seconds")
    df = df.sort_values(by=['Attributes', 'Start']).reset_index(drop=True)
    df.to_csv(args.output, index=False, sep='\t', header=False)
    
    
# Command for running the script
"""
singularity exec singularities/python.sif python -m cProfile ~/git/gitlab/Recoding/relocate_transcripts.py --gff /no_backup/rg/ileahy/recoding/gff/transcripts_for_recoding.gff --output /no_backup/rg/ileahy/recoding/gff/relocated_to_transcript.gff
singularity exec singularities/python.sif python -m cProfile ~/git/gitlab/Recoding/relocate_transcripts.py --gff /no_backup/rg/ileahy/recoding/gff/gencode.v47.annotation_clean.gtf --output /no_backup/rg/ileahy/recoding/gff/relocated_to_transcript.gff
"""