#!/usr/bin/env python3
# The purpose of this script is to extract the sequence around the recoded TGA site

import argparse
import os
import sys
import pandas as pd

def read_fasta_dir(fasta_dir, target_ids):
    seqs = {}
    for filename in os.listdir(fasta_dir):
        if filename.endswith('.fa'):
            filepath = os.path.join(fasta_dir, filename)
            with open(filepath) as f:
                current_id = None
                current_seq = []
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_id and current_id in target_ids:
                            seqs[current_id] = ''.join(current_seq)
                        current_id = line[1:].split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line)
                if current_id and current_id in target_ids:
                    seqs[current_id] = ''.join(current_seq)
    return seqs

def main():
    parser = argparse.ArgumentParser(description='Extract sequences around TGA sites for positive and negative transcripts')
    parser.add_argument('--result_file', required=True, help='SecORFsearch result file (CSV)')
    parser.add_argument('--fasta_dir', required=True, help='Directory with FASTA files from gffread')
    parser.add_argument('--output_dir', required=True, help='Output directory for FASTA files')
    args = parser.parse_args()

    # Read target list
    df = pd.read_csv(args.result_file)
    print(f"Loaded {len(df)} rows from result file", file=sys.stderr)

    # Get target transcript IDs (from df and target list)
    target_ids = set(df['transcript_name'].unique())
    print(f"Target IDs: {len(target_ids)}", file=sys.stderr)

    # Read sequences only for target IDs
    seqs = read_fasta_dir(args.fasta_dir, target_ids)
    print(f"Loaded {len(seqs)} sequences", file=sys.stderr)

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Open output files
    files = {
        'positives_603': open(os.path.join(args.output_dir, 'positives_603.fa'), 'w'),
        'negatives_603': open(os.path.join(args.output_dir, 'negatives_603.fa'), 'w'),
        'positives_stop_603': open(os.path.join(args.output_dir, 'positives_stop_603.fa'), 'w'),
        'negatives_stop_603': open(os.path.join(args.output_dir, 'negatives_stop_603.fa'), 'w'),
    }

    for idx, row in df.iterrows():
        tid = row['transcript_name']
        print(tid, file=sys.stderr)
        try:
            pos = int(row['TGA_site_score'])
            stop = int(row['re_score_end'])
            re_score = float(row['re_score_score'])
            og_score = float(row['og_score_score'])
            print(f"Parsed scores for {tid}: pos={pos}, stop={stop}, re={re_score}, og={og_score}", file=sys.stderr)
        except (ValueError, KeyError) as e:
            print(f"Warning: Failed to parse row for {tid}: {e}", file=sys.stderr)
            continue
        if pd.isna(re_score) or pd.isna(og_score):
            print(f"Warning: Skipping {tid} due to NA in scores", file=sys.stderr)
            continue
        if tid not in seqs:
            print(f"Warning: {tid} not found in FASTA", file=sys.stderr)
            continue
        print(f"Found sequence for {tid}, length: {len(seqs[tid])}", file=sys.stderr)
        seq = seqs[tid]
        
        # Check for TGA within +/-1 of annotated pos
        actual_pos = None
        if pos-1 >= 0 and pos+2 <= len(seq) and seq[pos-1:pos+2] == 'TGA':
            actual_pos = pos
        elif pos-2 >= 0 and pos+1 <= len(seq) and seq[pos-2:pos+1] == 'TGA':
            actual_pos = pos - 1
        elif pos >= 0 and pos+3 <= len(seq) and seq[pos:pos+3] == 'TGA':
            actual_pos = pos + 1
        if actual_pos is None:
            print(f"Warning: {tid} TGA not found within +/-1 at pos {pos}")
            continue  # Skip if TGA not found within +/-1
        print(f"Found TGA for {tid} at actual_pos {actual_pos}", file=sys.stderr)

        # Check for stop codon within +/-1 of annotated stop
        actual_stop = None
        if stop-1 >= 0 and stop+2 <= len(seq) and (seq[stop-1:stop+2] == 'TGA' or seq[stop-1:stop+2] == 'TAG' or seq[stop-1:stop+2] == 'TAA'):
            actual_stop = stop
        elif stop-2 >= 0 and stop+1 <= len(seq) and (seq[stop-2:stop+1] == 'TGA' or seq[stop-2:stop+1] == 'TAG' or seq[stop-2:stop+1] == 'TAA'):
            actual_stop = stop - 1
        elif stop >= 0 and stop+3 <= len(seq) and (seq[stop:stop+3] == 'TGA' or seq[stop:stop+3] == 'TAG' or seq[stop:stop+3] == 'TAA'):
            actual_stop = stop + 1
        if actual_stop is None:
            print(f"Warning: {tid} stop codon not found within +/-1 at stop {stop}")
            continue  # Skip if stop codon not found within +/-1
        print(f"Found stop for {tid} at actual_stop {actual_stop}", file=sys.stderr)

        # Extract -300 and +300 around full TGA codon (223bp total), padded with N.
        left_300 = actual_pos - 300
        right_300 = actual_pos + 3 + 300
        left_pad_300 = max(0, -left_300)
        right_pad_300 = max(0, right_300 - len(seq))
        start_300 = max(0, left_300)
        end_300 = min(len(seq), right_300)
        extracted_300 = 'N' * left_pad_300 + seq[start_300:end_300] + 'N' * right_pad_300
        assert extracted_300[300:303] == 'TGA'

        # Extract -300 and +300 around full stop codon (223bp total), padded with N.
        left_stop_300 = actual_stop - 300
        right_stop_300 = actual_stop + 3 + 300
        left_pad_stop_300 = max(0, -left_stop_300)
        right_pad_stop_300 = max(0, right_stop_300 - len(seq))
        start_stop_300 = max(0, left_stop_300)
        end_stop_300 = min(len(seq), right_stop_300)
        extracted_stop_300 = 'N' * left_pad_stop_300 + seq[start_stop_300:end_stop_300] + 'N' * right_pad_stop_300

        header = f">{tid}_{actual_pos}"
        header_stop = f">{tid}_{actual_stop}"
        is_positive = (re_score - og_score) >= -1.8 and re_score >= -1.5
        print(f"{tid} is_positive: {is_positive} (re={re_score}, og={og_score})", file=sys.stderr)
        if is_positive:
            files['positives_603'].write(f"{header}\n{extracted_300}\n")
            files['positives_stop_603'].write(f"{header_stop}\n{extracted_stop_300}\n")
            print(f"Wrote positive for {tid}", file=sys.stderr)
        else:
            files['negatives_603'].write(f"{header}\n{extracted_300}\n")
            files['negatives_stop_603'].write(f"{header_stop}\n{extracted_stop_300}\n")
            print(f"Wrote negative for {tid}", file=sys.stderr)

    # Close all files
    for f in files.values():
        f.close()

if __name__ == '__main__':
    main()
    
# Code for running the script
"""
python bin/extract_sequence_logos.py --result_file /Users/iseult/Desktop/Human_Analysis/SecORFsearch_results.csv --fasta_dir /no_backup/rg/ileahy/Human_Analysis/gffread/gffread --output_dir /no_backup/rg/ileahy/Human_Analysis/secorfsearch_out/sequence_logos
"""