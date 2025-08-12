#!/usr/bin/env python3

import sys
import pandas as pd
from collections import defaultdict
from typing import List, Dict, Tuple

def remove_duplicate_headers(input_file: str, output_file: str) -> str:
    """
    Remove duplicate header lines from the input file, keeping only the first instance.
    
    Args:
        input_file: Path to the input file
        output_file: Path to save the output file
        
    Returns:
        Path to the output file with deduplicated headers
    """
    header_seen = False
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('transcript_name,') or line.startswith('transcript_name\t'):
                if not header_seen:
                    outfile.write(line)
                    header_seen = True
            else:
                outfile.write(line)
    return output_file

def find_duplicate_transcripts(input_file: str) -> Dict[str, List[int]]:
    """
    Find all transcript names that appear multiple times in the file.
    
    Args:
        input_file: Path to the input file with transcript data
        
    Returns:
        Dictionary mapping transcript names to lists of line numbers (0-based) where they appear
    """
    transcript_lines = defaultdict(list)
    with open(input_file, 'r') as f:
        # Skip header
        next(f)
        for i, line in enumerate(f, 1):  # Start from 1 to match line numbers
            if line.strip():
                transcript = line.split(',', 1)[0].strip('"\'')
                if transcript:  # Skip empty transcript names
                    transcript_lines[transcript].append(i)
    
    # Return only transcripts that appear more than once
    return {k: v for k, v in transcript_lines.items() if len(v) > 1}

def select_best_predictions(input_file: str, output_file: str) -> None:
    """
    Process the input file to keep only the best predictions for each transcript.
    For each transcript, selects the row with the highest score for both original and recoded predictions.
    
    Args:
        input_file: Path to the input file with transcript data
        output_file: Path to save the filtered output
    """
    # Read the data
    df = pd.read_csv(input_file)
    
    # Group by transcript_name and select the row with the best scores
    def get_best_row(group):
        # Calculate a combined score for ranking
        group['combined_score'] = (
            group['og_score_score'].fillna(0) + 
            group['re_score_score'].fillna(0)
        )
        # Get the row with the maximum combined score
        best_idx = group['combined_score'].idxmax()
        return group.loc[best_idx].drop('combined_score')
    
    # Apply the function to each group and reset index
    result_df = df.groupby('transcript_name').apply(get_best_row).reset_index(drop=True)
    
    # Save the result
    result_df.to_csv(output_file, index=False)

def main():
    if len(sys.argv) != 3:
        print("Usage: python filter_final_table.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Step 1: Remove duplicate headers
    print("Removing duplicate headers...")
    temp_file = "temp_dedup_headers.csv"
    remove_duplicate_headers(input_file, temp_file)
    
    # Step 2: Find and process duplicates
    print("Finding duplicate transcripts...")
    duplicates = find_duplicate_transcripts(temp_file)
    
    if duplicates:
        print(f"Found {len(duplicates)} transcripts with multiple entries")
        print("Selecting best predictions for each transcript...")
        select_best_predictions(temp_file, output_file)
        print(f"Filtered results saved to {output_file}")
    else:
        print("No duplicate transcripts found. Copying file as is...")
        import shutil
        shutil.copy2(temp_file, output_file)
    
    # Clean up
    import os
    if os.path.exists(temp_file):
        os.remove(temp_file)
    
    print("Done!")

if __name__ == "__main__":
    main()