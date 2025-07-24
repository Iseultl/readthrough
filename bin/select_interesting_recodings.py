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

def select_target_predictions(geneid_txt):
    data = []
    total_predictions = 0
    out_of_frame = 0
    added_empty = set()

    with open(geneid_txt) as f:
        current_seq = None
        current_tga_pos = None
        found_prediction = False
        prediction_count = 0

        for line in f:
            if line.startswith("# Sequence"):
                # Add empty row for previous sequence if needed
                if current_seq and not found_prediction and "original" not in current_seq and current_seq not in added_empty:
                    data.append([
                        '1', np.nan, np.nan, '+', np.nan, np.nan,
                        current_seq, current_seq, 'protein_coding'
                    ])
                    added_empty.add(current_seq)

                # Start new sequence
                current_seq = line.split()[2]
                current_tga_pos = extract_tga_from_header(current_seq)
                found_prediction = False
                prediction_count = 0

            elif line.startswith("  Single"):
                total_predictions += 1
                fields = line.strip().split()
                start = int(fields[1])
                end = int(fields[2]) - 3

                if "original" in current_seq:
                    continue  # skip "original" sequences

                if is_in_frame(current_tga_pos, start, end):
                    found_prediction = True
                    structure = fields[0]
                    score = float(fields[3])
                    strand = '+'
                    seqname = current_seq + '+' + str(prediction_count)
                    gene_name = current_seq

                    data.append([
                        '1', start, end, strand, structure, score,
                        gene_name, seqname, 'protein_coding'
                    ])
                    prediction_count += 1
                else:
                    out_of_frame += 1

        # Handle last sequence if it had no predictions
        if current_seq and not found_prediction and "original" not in current_seq and current_seq not in added_empty:
            data.append([
                '1', np.nan, np.nan, '+', np.nan, np.nan,
                current_seq, current_seq, 'protein_coding'
            ])
            added_empty.add(current_seq)

    df = pd.DataFrame(data, columns=[
        'seqnames', 'start', 'end', 'strand', 'type', 'score',
        'gene_name', 'transcript_name', 'transcript_biotype'
    ])
    
    # Optionally add a length column if needed
    df["length"] = df["end"] - df["start"]

    print("Total Geneid Predictions:", total_predictions)
    print("Filtered Out of Frame:", out_of_frame)
    print("Total Transcripts in Output:", len(df['transcript_name'].unique()))
    
    return df


def get_length(df):
    df['length'] = (df['end'] - df['start']) + 1
    return df

def transcript_name(df):
    df['transcript_name'] = df['gene_name'].str.split('_').str[0]
    return df

def best_score(df):
    # Separate rows with valid scores
    df_valid = df.dropna(subset=["score"])
    # Compute best scores from valid rows
    score_idx = df_valid.groupby("transcript_name")["score"].idxmax()
    score_df = df_valid.loc[score_idx]

    # Identify transcript_names with no predictions
    all_transcripts = set(df["transcript_name"])
    scored_transcripts = set(score_df["transcript_name"])
    missing_transcripts = all_transcripts - scored_transcripts

    # Include original rows for missing predictions
    missing_rows = df[df["transcript_name"].isin(missing_transcripts)]
    missing_rows_unique = missing_rows.drop_duplicates(subset=["transcript_name"], keep="first")
     
    # Combine and return
    combined_df = pd.concat([score_df, missing_rows_unique], ignore_index=True)
    return combined_df.reset_index(drop=True)

def best_length(df):
    # Separate rows with valid lengths
    df_valid = df.dropna(subset=["length"])
    # Compute longest predictions
    longest_idx = df_valid.groupby("transcript_name")["length"].idxmax()
    longest_df = df_valid.loc[longest_idx]

    # Identify transcript_names with no predictions
    all_transcripts = set(df["transcript_name"])
    length_transcripts = set(longest_df["transcript_name"])
    missing_transcripts = all_transcripts - length_transcripts

    # Include only one row per missing transcript (take the first)
    missing_rows = df[df["transcript_name"].isin(missing_transcripts)]
    missing_rows_unique = missing_rows.drop_duplicates(subset=["transcript_name"], keep="first")

    # Combine and return
    combined_df = pd.concat([longest_df, missing_rows_unique], ignore_index=True)
    return combined_df.reset_index(drop=True)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Select interesting transcript sequences")
    # Add arguments for input and output directories
    parser.add_argument('--geneid_scores', type=str, required=True, help="Path to geneid scores of all ORFs")
    parser.add_argument('--score', type=str, required=True, help="Path to output file")
    parser.add_argument('--longest', type=str, required=True, help="Path to output file")
    parser.add_argument('--gff', type=str, required=True, help="Path to relocated to transcript gff")
    args = parser.parse_args()
    
    df = select_target_predictions(args.geneid_scores)
    
    df = get_length(df)
    df = transcript_name(df)
    
    print(df[df["transcript_name"] == 'C1orf54-203'])

    score = best_score(df)
    print(score)
    longest = best_length(df)
    print(longest)
    score.to_csv(args.score + '_score.csv', index=False)
    longest.to_csv(args.longest + '_longest.csv', index=False)
    
# Command for running script
"""
python select_interesting_recodings.py --geneid_scores /Users/iseult/Desktop/Geneid_Recoding/testing_false_positives/geneid_results.txt --score temp_test --longest temp_test --gff /Users/iseult/Desktop/Geneid_Recoding/testing_false_positives/relocated_to_transcript.gff
""" 