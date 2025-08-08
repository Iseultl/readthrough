#!/usr/bin/env python3

# The purpose of this script is to select the in-frame recodings 
# extract their score and length and then output all the best 
# recodings to the score/longest table

import pandas as pd
import re
import argparse
import numpy as np



def extract_tga_from_header(header):
    match = re.search(r'_TGA_(\d+)_', header)
    return int(match.group(1)) if match else None

def is_in_frame(tga_pos, start, end):
    if tga_pos is not None:
        return start <= tga_pos <= end and ((tga_pos + 1) - start) % 3 == 0
    else:
        return False

def process_recoding_predictions(geneid_txt):
    best_by_score = {}
    best_by_length = {}
    total_predictions = 0
    out_of_frame = 0
    
    def get_transcript_name(seq_name):
        if "_TGA" not in seq_name:
            return None
        else:
            return re.search(r'^(.*?)_TGA', seq_name).group(1)
    
    def get_gene_name(seq_name):
        if "_TGA" not in seq_name:
            return None
        else:
            return re.search(r'^(.*?)_TGA', seq_name).group(1).split('-')[0]

    with open(geneid_txt) as f:
        current_seq = None
        current_tga_pos = None

        for line in f:
            if line.startswith("# Sequence"):
                prediction_count = 0
                if current_seq and "original" not in current_seq:
                    print(current_seq)
                    transcript_name = get_transcript_name(current_seq)
                    gene_name = get_gene_name(transcript_name)
                    if transcript_name not in best_by_score:
                        # No valid in-frame predictions found, add placeholder
                        placeholder = {
                            'seqnames': current_seq + "_" + str(prediction_count), 'start': np.nan, 'end': np.nan, 'strand': '+',
                            'type': np.nan, 'score': np.nan, 'gene_name': gene_name,
                            'transcript_name': transcript_name, 'transcript_biotype': 'protein_coding',
                            'length': np.nan
                        }
                        best_by_score[transcript_name] = placeholder
                        best_by_length[transcript_name] = placeholder

                current_seq = line.split()[2]
                if "original" not in current_seq:
                    current_tga_pos = extract_tga_from_header(current_seq)
                    
                else:
                    current_tga_pos = None

            elif line.startswith("  Single"):
                
                total_predictions += 1
                fields = line.strip().split()
                start = int(fields[1])
                end = int(fields[2]) - 3

                if is_in_frame(current_tga_pos, start, end):
                    prediction_count += 1
                    transcript_name = get_transcript_name(current_seq)
                    gene_name = get_gene_name(transcript_name)
                    prediction = {
                        'seqnames': current_seq + "_" + str(prediction_count),
                        'start': start,
                        'end': end,
                        'strand': '+',
                        'type': fields[0],
                        'score': float(fields[3]),
                        'gene_name': gene_name,
                        'transcript_name': transcript_name,
                        'transcript_biotype': 'protein_coding',
                        'length': abs(end - start) + 1
                    }
                     
                    # Update best by score
                    if transcript_name not in best_by_score or pd.isna(best_by_score[transcript_name]['score']) or prediction['score'] > best_by_score[transcript_name]['score']:
                        best_by_score[transcript_name] = prediction

                    # Update best by length
                    current_best_len = best_by_length.get(transcript_name, {}).get('length')
                    if current_best_len is None or pd.isna(current_best_len) or prediction['length'] > current_best_len:
                        best_by_length[transcript_name] = prediction
                else:
                    out_of_frame += 1

        # Handle the very last sequence if it has no predictions 
        if current_seq and "original" not in current_seq:
            transcript_name = get_transcript_name(current_seq)
            gene_name = get_gene_name(transcript_name)
            if transcript_name not in best_by_score:
                placeholder = {
                    'seqnames': current_seq + "_" + str(prediction_count), 'start': np.nan, 'end': np.nan, 'strand': '+',
                    'type': np.nan, 'score': np.nan, 'gene_name': gene_name,
                    'transcript_name': transcript_name, 'transcript_biotype': 'protein_coding',
                    'length': np.nan
                }
                best_by_score[transcript_name] = placeholder
                best_by_length[transcript_name] = placeholder

    score_df = pd.DataFrame(list(best_by_score.values()))
    longest_df = pd.DataFrame(list(best_by_length.values()))

    print("Total Geneid Predictions (non-original):", total_predictions)
    print("Filtered Out of Frame:", out_of_frame)
    print("Total Transcripts in Score Output:", len(score_df))
    print("Total Transcripts in Longest Output:", len(longest_df))

    return score_df, longest_df



    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Select interesting transcript sequences")
    # Add arguments for input and output directories
    parser.add_argument('--geneid', type=str, required=True, help="Path to geneid scores of all ORFs")
    parser.add_argument('--score', type=str, required=True, help="Path to output file")
    parser.add_argument('--longest', type=str, required=True, help="Path to output file")
    parser.add_argument('--gff', type=str, required=True, help="Path to relocated to transcript gff")
    args = parser.parse_args()
    
    score_df, longest_df = process_recoding_predictions(args.geneid)

    score_df.to_csv(args.score, index=False)
    longest_df.to_csv(args.longest, index=False)
    
# Command for running script
"""
python select_interesting_recodings.py --geneid_scores /Users/iseult/Desktop/Geneid_Recoding/testing_false_positives/geneid_results.txt --score temp_test --longest temp_test --gff /Users/iseult/Desktop/Geneid_Recoding/testing_false_positives/relocated_to_transcript.gff
python select_interesting_recodings.py --geneid /no_backup/rg/ileahy/Horse_Analysis/secis_independent_output/geneid_original_predictions/geneid_original_predictions/transcripts_clean_NC_091684.1_recoded.part_001.fa.geneid.txt --score temp_test --longest temp_test --gff ~/git/gitlab/secis_independent/relocated_gtf/relocated_gtf/relocated.gtf 
""" 