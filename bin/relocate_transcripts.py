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


def extract_feature_ids(attributes, feature):
    gene_id = pd.NA
    transcript_id = pd.NA

    if pd.isna(attributes):
        return gene_id, transcript_id

    gene_match = pd.Series([attributes]).str.extract(r'gene_id "([^"]+)"')[0].iloc[0]
    if pd.notna(gene_match):
        gene_id = gene_match
    else:
        gene_match = pd.Series([attributes]).str.extract(r'(?:^|;)gene=([^;]+)')[0].iloc[0]
        if pd.notna(gene_match):
            gene_id = gene_match

    transcript_match = pd.Series([attributes]).str.extract(r'transcript_id "([^"]+)"')[0].iloc[0]
    if pd.notna(transcript_match):
        transcript_id = transcript_match
    else:
        transcript_match = pd.Series([attributes]).str.extract(r'(?:^|;)transcript=([^;]+)')[0].iloc[0]
        if pd.notna(transcript_match):
            transcript_id = transcript_match
        else:
            id_match = pd.Series([attributes]).str.extract(r'(?:^|;)ID=([^;]+)')[0].iloc[0]
            if pd.notna(id_match):
                if feature in {'gene'}:
                    gene_id = id_match
                elif feature in {'transcript', 'mRNA', 'mrna'}:
                    transcript_id = id_match
                else:
                    transcript_id = id_match

    parent_match = pd.Series([attributes]).str.extract(r'(?:^|;)Parent=([^;]+)')[0].iloc[0]
    if pd.notna(parent_match):
        if feature in {'transcript', 'mRNA', 'mrna'} and pd.isna(gene_id):
            gene_id = parent_match
        elif feature in {'exon', 'CDS', 'start_codon', 'stop_codon'} and pd.isna(transcript_id):
            transcript_id = parent_match

    if pd.notna(gene_id):
        gene_id = str(gene_id).replace('gene-', '').replace('rna-', '')
    if pd.notna(transcript_id):
        transcript_id = str(transcript_id).replace('transcript-', '').replace('rna-', '').replace('exon-', '').replace('cds-', '')

    return gene_id, transcript_id

# Load GFF file into dataframe
def load_gff(gff_file):
    df = pd.read_csv(gff_file, sep='\t', comment='#', header=None)
    df.columns = ['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attributes']
    df['Start'] -= 1  # Convert to 0-based indexing

    ids = df.apply(lambda row: extract_feature_ids(row['Attributes'], row['Feature']), axis=1, result_type='expand')
    ids.columns = ['gene_id', 'transcript_id']

    df['gene_id'] = ids['gene_id']
    df['transcript_id'] = ids['transcript_id']

    transcript_mask = df['Feature'].isin(['transcript', 'mRNA', 'mrna'])
    transcript_gene_map = (
        df.loc[transcript_mask & df['gene_id'].notna() & df['transcript_id'].notna(), ['transcript_id', 'gene_id']]
        .drop_duplicates(subset=['transcript_id'])
        .set_index('transcript_id')['gene_id']
    )
    missing_gene_mask = df['gene_id'].isna() & df['transcript_id'].notna()
    df.loc[missing_gene_mask, 'gene_id'] = df.loc[missing_gene_mask, 'transcript_id'].map(transcript_gene_map)

    # Keep only rows where we can map both levels; downstream code groups on this key.
    df = df.dropna(subset=['gene_id', 'transcript_id']).copy()

    df['Attributes'] = df['gene_id'].astype(str) + " ; " + df['transcript_id'].astype(str)
 
    return df

def load_secis(secis_file):
    df = pd.read_csv(secis_file, sep='\t', comment='#', header=None)
    df.columns = ['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attributes']
    df['Start'] -= 1  # Convert to 0-based indexing
    return df


def group_attributes(gff):
    sorted_gff = gff.sort_values(by=["Attributes", "Start"]).copy()
    grouped_by_attr = dict(tuple(sorted_gff.groupby("Attributes")))
    return grouped_by_attr 

# First get the new total transcript lengths
def relocate_trancripts(gff):
    exon_df = gff[gff['Feature'] == 'exon'].copy()
    if exon_df.empty:
        return []
    exon_df['length'] = exon_df['End'] - exon_df['Start']
    transcript_info = exon_df.groupby('Attributes').agg({'length': 'sum', 'Strand': 'first', 'Chromosome': 'first'})
    relocated_transcripts = []
    for attr, row in transcript_info.iterrows():
        transcript_dct = {
            "Chromosome": row['Chromosome'],
            "Source": '.',
            "Feature": "transcript",
            "Start": 0,
            "End": row['length'],
            "Score": '.',
            "Strand": row['Strand'],
            "Frame": '.',
            "Attributes": attr
        }
        relocated_transcripts.append(transcript_dct)
    return relocated_transcripts

# Second get the CDS locations
def relocate_CDS(gff):
    exon_df = gff[gff['Feature'] == 'exon'].copy()
    CDS_df = gff[gff['Feature'] == 'CDS'].copy()
    if exon_df.empty or CDS_df.empty:
        return [], []

    exon_df['length'] = exon_df['End'] - exon_df['Start']
    CDS_df['length'] = CDS_df['End'] - CDS_df['Start']

    # Get attrs with both exons and CDS
    attrs_with_exons = set(exon_df['Attributes'])
    attrs_with_cds = set(CDS_df['Attributes'])
    valid_attrs = attrs_with_exons & attrs_with_cds

    exon_df = exon_df[exon_df['Attributes'].isin(valid_attrs)]
    CDS_df = CDS_df[CDS_df['Attributes'].isin(valid_attrs)]

    # Sort exons and compute cumsum for exon positions
    exon_df = exon_df.sort_values(['Attributes', 'Start'])
    exon_df['cumsum'] = exon_df.groupby('Attributes')['length'].cumsum()
    exon_df['Start_new'] = exon_df['cumsum'] - exon_df['length'] + 1
    exon_df['End_new'] = exon_df['cumsum']

    # Get min CDS start per attr
    min_cds_start = CDS_df.groupby('Attributes')['Start'].min()

    # Compute UTR5
    exon_df = exon_df.merge(min_cds_start.rename('min_cds_start'), left_on='Attributes', right_index=True, how='left')
    exon_df['utr_contrib'] = 0
    mask_full = exon_df['End'] <= exon_df['min_cds_start']
    exon_df.loc[mask_full, 'utr_contrib'] = exon_df.loc[mask_full, 'length']
    mask_partial = (exon_df['Start'] < exon_df['min_cds_start']) & (exon_df['min_cds_start'] < exon_df['End'])
    exon_df.loc[mask_partial, 'utr_contrib'] = exon_df.loc[mask_partial, 'min_cds_start'] - exon_df.loc[mask_partial, 'Start']
    utr5 = exon_df.groupby('Attributes')['utr_contrib'].sum()

    # Relocate exons
    relocated_exons = []
    for _, row in exon_df.iterrows():
        exon_dct = {
            "Chromosome": row['Chromosome'],
            "Source": '.',
            "Feature": "exon",
            "Start": row['Start_new'],
            "End": row['End_new'],
            "Score": '.',
            "Strand": row['Strand'],
            "Frame": '.',
            "Attributes": row['Attributes']
        }
        relocated_exons.append(exon_dct)

    # Relocate CDS
    CDS_df = CDS_df.sort_values(['Attributes', 'Start'])
    CDS_df['cumsum'] = CDS_df.groupby('Attributes')['length'].cumsum()
    CDS_df = CDS_df.merge(utr5.rename('utr5'), left_on='Attributes', right_index=True, how='left')
    CDS_df['Start_new'] = CDS_df['cumsum'] - CDS_df['length'] + CDS_df['utr5'] + 1
    CDS_df['End_new'] = CDS_df['cumsum'] + CDS_df['utr5']

    relocated_CDS = []
    for _, row in CDS_df.iterrows():
        CDS_dct = {
            "Chromosome": row['Chromosome'],
            "Source": '.',
            "Feature": "CDS",
            "Start": row['Start_new'],
            "End": row['End_new'],
            "Score": '.',
            "Strand": row['Strand'],
            "Frame": '.',
            "Attributes": row['Attributes']
        }
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
            for _, sec_row in sec_rows.iterrows():
                sec_pos = 0
                exon_pos = 0
                for _, exon_row in exon_rows.iterrows():
                    if (exon_row['Start'] < sec_row['Start']) and (sec_row['Start'] < exon_row['End']):
                        exon_pos = sec_row['Start'] - exon_row['Start']
                        break
                    else:
                        exon_len = exon_row['End'] - exon_row['Start']
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

    for attr, transcript_group in grouped.items():
        # Extract rows
        secis_rows = transcript_group[transcript_group["Feature"] == "secis"]
        if secis_rows.empty:
            continue

        strand = secis_rows.iloc[0]["Strand"]
        # Since SECIS are already on transcript coordinates, copy directly
        for _, row in secis_rows.iterrows():
            start = int(row['Start'])
            end = int(row['End'])

            secis_dct = {
                "Chromosome": row['Chromosome'],
                "Source": '.',
                "Feature": "secis",
                "Start": start,
                "End": end,
                "Score": '.',
                "Strand": strand,
                "Frame": '.',
                "Attributes": attr  # Use the group key as attributes
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
    features_to_adjust = ['Selenocysteine', 'CDS', 'exon']
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

def edit_secis(df):
    new_df = df.copy()
    new_df['transcript_id'] = new_df['Chromosome']
    new_df['Chromosome'] = None
    return new_df
    
def concat_secis_gtf(secis, gtf):
    # Copy inputs
    gtf_copy = gtf.copy()
    secis_copy = secis.copy()
    # Extract transcript IDs and gene IDs from Attributes
    if 'gene_id' not in gtf_copy.columns or 'transcript_id' not in gtf_copy.columns:
        gtf_copy[['gene_id', 'transcript_id']] = gtf_copy['Attributes'].str.split(' ; ', n=1, expand=True)
    secis_copy = secis_copy.drop(columns=['Chromosome'])
    # Merge SECIS with GTF gene info based on transcript ID
    merged = secis_copy.merge(
        gtf_copy[['gene_id', 'transcript_id', 'Chromosome']],
        on='transcript_id',
        how='left'
    )
    print(merged.head()) 
    # If gene_id is missing (no match found), keep original Attributes
    merged['gene_id'] = merged['gene_id'].ffill()
    merged['Chromosome'] = merged['Chromosome'].ffill()
    merged['Attributes'] = merged['gene_id'].astype(str) + ' ; ' + merged['transcript_id'].astype(str)
    
    # Clean up
    merged = merged.drop(columns=['gene_id', 'transcript_id'])
    merged = merged.drop_duplicates()
    return merged.to_dict('records')
                 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Recode TGA codons at selenocysteine positions in a genome FASTA file")
    parser.add_argument('--gff', type=str, required=True, help="Path to transcript GFF file")
    parser.add_argument('--secis', type=str, required=False, help="Path to secis GFF file")
    parser.add_argument('--output', type=str, required=True, help="Path to output gtf file - positioned on the spliced transcript")
    args = parser.parse_args()
    
    start = time.time()
    gff = load_gff(args.gff)
    groups = group_attributes(gff)
    transcript_lst = relocate_trancripts(gff)
    print("Relocate transcripts took", time.time() - start, "seconds")
    relocated_CDS, relocated_exons = relocate_CDS(gff)
    print("Relocate CDS took", time.time() - start, "seconds")
    relocated_sec = relocate_sec(groups)
    print("Relocate sec took", time.time() - start, "seconds")
    if args.secis:
        secis_gff = load_gff(args.secis)
        secis_gff = edit_secis(secis_gff)
        relocated_secis = concat_secis_gtf(secis_gff, gff)
    else:
        relocated_secis = []
    print("Relocate secis took", time.time() - start, "seconds")

    all_dcts = transcript_lst+relocated_CDS+relocated_exons+relocated_sec+relocated_secis
    if not all_dcts:
        print("No rows to write")
        df = pd.DataFrame()
    else:
        df = pd.DataFrame(all_dcts)
 
    df = handle_negs(df)
    print("Relocate negatives took", time.time() - start, "seconds")
    print(df.head())
    df = df.sort_values(by=['Attributes', 'Start']).reset_index(drop=True)
    print("Sort took", time.time() - start, "seconds")
    df.to_csv(args.output, index=False, sep='\t', header=False)
    print("Write took", time.time() - start, "seconds")
    
    
# Command for running the script
"""
singularity exec ~/singularities/python.sif python -m cProfile ~/git/gitlab/Recoding/relocate_transcripts.py --gff /no_backup/rg/ileahy/recoding/gff/transcripts_for_recoding.gff --output /no_backup/rg/ileahy/recoding/gff/relocated_to_transcript.gff
singularity exec ~/singularities/python.sif python -m cProfile ~/git/gitlab/Recoding/relocate_transcripts.py --gff /no_backup/rg/ileahy/recoding/gff/gencode.v47.annotation_clean.gtf --output /no_backup/rg/ileahy/recoding/gff/relocated_to_transcript.gff
"""