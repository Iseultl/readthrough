#!/usr/bin/env python3

# The purpose of this script is to select the 
# best predictions on the original transcript sequences

import pandas as pd
import argparse
import re
import numpy as np

def extract_predictions(geneid_txt):
    predictions = []
    current_seq = None
    found_prediction = False  # Track if any prediction was found for the current sequence

    with open(geneid_txt) as f:
        for line in f:
            if line.startswith("# Sequence"):
                # Before moving to the next sequence, check if the previous one had no predictions
                if current_seq and not found_prediction:
                    predictions.append({
                        'seq': current_seq,
                        'start': np.nan,
                        'end': np.nan,
                        'strand': np.nan,
                        'type': np.nan,
                        'length': np.nan,
                        'score': np.nan,
                        'gene_name': current_seq.split('-')[0]
                    })

                current_seq = line.split()[2]
                found_prediction = False  # Reset for new sequence

            elif current_seq and "Single " in line and current_seq not in line:
                # Parse prediction
                fields = line.strip().split()
                structure = fields[0]
                start = int(fields[1])
                end = int(fields[2]) - 3
                score = float(fields[3])
                aa_seq = fields[-1]
                predictions.append({
                    'seq': current_seq,
                    'start': start,
                    'end': end,
                    'strand': '+',
                    'type': structure,
                    'length': abs(end - start) + 1,
                    'score': score,
                    'gene_name': '-'.join(current_seq.split('-')[:2])
                })
                found_prediction = True

        # After the loop ends, make sure to check the last sequence
        if current_seq and "original" in current_seq and not found_prediction:
            predictions.append({
                'seq': current_seq,
                'start': np.nan,
                'end': np.nan,
                'strand': np.nan,
                'type': np.nan,
                'length': np.nan,
                'score': np.nan,
                'gene_name': current_seq.split('-')[0]
            })

    df = pd.DataFrame.from_dict(predictions)
    return df


def get_transcript_id(df):
    df['transcript_name'] = df['seq'].str.split('_').str[0]
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
    
    # Combine and return
    combined_df = pd.concat([score_df, missing_rows], ignore_index=True)
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

    # Include original rows for missing predictions
    missing_rows = df[df["transcript_name"].isin(missing_transcripts)]
    
    # Combine and return
    combined_df = pd.concat([longest_df, missing_rows], ignore_index=True)
    return combined_df.reset_index(drop=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract the coding regions of each transcript and input to coding_potential script")
    # Add arguments for input and output directories
    parser.add_argument('--score', type=str, required=True, help="Path to output file - csv of predictions")
    parser.add_argument('--longest', type=str, required=True, help="Path to output file - csv of predictions")
    parser.add_argument('--geneid', type=str, required=True, help="Path to geneid txt")
     
    args = parser.parse_args()
    
    geneid = extract_predictions(args.geneid)
    print(geneid)
    geneid = get_transcript_id(geneid)
    
    score = best_score(geneid)
    longest = best_length(geneid)
    longest.to_csv(args.longest + '_longest.csv', index=False)
    score.to_csv(args.score + '_score.csv', index=False)
    
       
    
    

# Code for running script
"""
python get_og_predictions.py --geneid /Users/iseult/Desktop/Geneid_Recoding/testing_false_positives/geneid_results_originals.txt --output temp_test_og
"""