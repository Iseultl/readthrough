#!/usr/bin/env python3
# The purpose of this script is to extract the sequence around the recoded TGA site

import argparse
import os
import sys
import pandas as pd

STOP_CODONS = {"TAA", "TAG", "TGA"}

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

######
# Check positions are stop or sec codons
######

def _codon_at(seq, start_1based):
    start_0 = start_1based - 1
    if start_0 < 0 or start_0 + 3 > len(seq):
        return None
    return seq[start_0 : start_0 + 3]


def infer_sec_start(seq,
                    sec_raw_1based):

    refs = [sec_raw_1based]

    offsets = [0, 1, -1]
    candidates = []

    for ref in refs:
        for off in offsets:
            pos = ref + off
            codon = _codon_at(seq, pos)

            if codon == "TGA":
                candidates.append(
                    (
                        abs(off),
                        abs(pos - sec_raw_1based),
                        pos
                    )
                )

    if candidates:
        candidates.sort()
        return candidates[0][2], "TGA"

    return None, None


def infer_stop_from_cds_end(seq, cds_end_1based, sec_start_1based):
    """Infer stop-codon start using CDS-end guidance with robust offset handling."""
    refs = [cds_end_1based]

    stop_positions = [
        i + 1 for i in range(len(seq) - 2) if seq[i : i + 3] in STOP_CODONS
    ]
    inframe_after_sec = [
        p for p in stop_positions if p > sec_start_1based and (p - sec_start_1based) % 3 == 0
    ]

    offsets = [1, -1, 0]
    direct_candidates = []
    for ref in refs:
        for off in offsets:
            pos = ref + off
            codon = _codon_at(seq, pos)
            if codon in STOP_CODONS and pos > sec_start_1based:
                inframe_penalty = 0 if (pos - sec_start_1based) % 3 == 0 else 1
                direct_candidates.append((inframe_penalty, abs(off), abs(pos - ref), pos, codon))

    if direct_candidates:
        direct_candidates.sort()
        _, _, _, pos, codon = direct_candidates[0]
        return pos, codon

    if inframe_after_sec:
        target_ref = refs[-1]
        pos = min(inframe_after_sec, key=lambda p: abs(p - target_ref))
        return pos, _codon_at(seq, pos)

    downstream_stops = [p for p in stop_positions if p > sec_start_1based]
    if downstream_stops:
        target_ref = refs[-1]
        pos = min(downstream_stops, key=lambda p: abs(p - target_ref))
        return pos, _codon_at(seq, pos)

    return None, None

def extract_window(seq, codon_start_1based, upstream, downstream):
    """Extract window around codon start with N padding.

    Window length is upstream + 3 + downstream.
    """
    codon_start_0 = codon_start_1based - 1
    left = codon_start_0 - upstream
    right = codon_start_0 + 3 + downstream

    left_pad = max(0, -left)
    right_pad = max(0, right - len(seq))

    start = max(0, left)
    end = min(len(seq), right)
    return ("N" * left_pad) + seq[start:end] + ("N" * right_pad)


def write_fasta(path, records):
    with open(path, "w", encoding="utf-8") as handle:
        for header, seq in records:
            handle.write(f">{header}\n{seq}\n")

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


    sec_pos_records = []
    sec_neg_records = []
    stop_pos_records = []
    stop_neg_records = []

    for row in df.itertuples(index=False):
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
        actual_pos, codon = infer_sec_start(seq, pos)
        if actual_pos is None:
            print(f"Warning: {tid} TGA not found within +/-1 at pos {pos}")
            continue  # Skip if TGA not found within +/-1
        print(f"Found TGA for {tid} at actual_pos {actual_pos}", file=sys.stderr)

        # Check for stop codon within +/-1 of annotated stop
        actual_stop, stop_codon = infer_stop_from_cds_end(seq, stop, actual_pos)
        if actual_stop is None:
            print(f"Warning: {tid} stop codon not found within +/-1 at stop {stop}")
            continue  # Skip if stop codon not found within +/-1
        print(f"Found stop for {tid} at actual_stop {actual_stop}", file=sys.stderr)

        # Extract sequences
        extracted_300 = extract_window(seq, actual_pos, 300, 300)
        extracted_stop_300 = extract_window(seq, actual_stop, 300, 300)
        
        
        
        is_positive = (re_score - og_score) >= -1.8 and re_score >= -1.5
        print(f"{tid} is_positive: {is_positive} (re={re_score}, og={og_score})", file=sys.stderr)
                
        if is_positive:
            sec_pos_records.append((f"{tid}|sec_start={actual_pos}|codon={codon}", extracted_300))
            stop_pos_records.append((f"{tid}|sec_start={actual_stop}|codon={stop_codon}", extracted_stop_300))       
        else:
            sec_neg_records.append((f"{tid}|sec_start={actual_pos}|codon={codon}", extracted_300))
            stop_neg_records.append((f"{tid}|sec_start={actual_stop}|codon={stop_codon}", extracted_stop_300))

    # Write to files
    write_fasta(os.path.join(args.output_dir, 'positives_603.fa'), sec_pos_records)
    write_fasta(os.path.join(args.output_dir, 'negatives_603.fa'), sec_neg_records)
    write_fasta(os.path.join(args.output_dir, 'positives_stop_603.fa'), stop_pos_records)
    write_fasta(os.path.join(args.output_dir, 'negatives_stop_603.fa'), stop_neg_records)
    


if __name__ == '__main__':
    main()
    
# Code for running the script
"""
python bin/extract_sequence_logos.py --result_file /Users/iseult/Desktop/Human_Analysis/SecORFsearch_results.csv --fasta_dir /no_backup/rg/ileahy/Human_Analysis/gffread/gffread --output_dir /no_backup/rg/ileahy/Human_Analysis/secorfsearch_out/sequence_logos
"""