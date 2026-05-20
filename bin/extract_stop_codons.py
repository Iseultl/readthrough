#!/usr/bin/env python3
"""
Extract and analyze stop codons from coding transcripts.

This script:
1. Reads a GFF file with relocated CDS coordinates (from relocate_transcripts.py)
2. Reads transcript sequences (from gffread)
3. Extracts the stop codon from the last CDS position of each transcript
4. Prints distribution of stop codons across the genome
"""

import argparse
from collections import Counter
from operator import pos
from Bio import SeqIO
import pandas as pd
import os
from pathlib import Path


def load_gff(gff_file):
    """Load GFF file into a dataframe."""
    df = pd.read_csv(gff_file, sep='\t', comment='#', header=None)
    df.columns = [
        'Chromosome', 'Source', 'Feature', 'Start', 'End', 
        'Score', 'Strand', 'Frame', 'Attributes'
    ]
    return df


def load_sequences_from_directory(gffread_dir):
    """
    Load transcript sequences from chromosome-specific FASTA files in a directory.
    Returns a function that lazily loads sequences from the correct chromosome file.
    """
    cache = {}
    
    def get_sequence(transcript_id, chromosome):
        """Get sequence for a transcript from the appropriate chromosome file."""
        # Try to load from cache first
        if chromosome not in cache:
            # Build filename for this chromosome
            fasta_file = os.path.join(gffread_dir, f'transcripts_clean_{chromosome}.fa')
            
            if not os.path.exists(fasta_file):
                return None
            
            # Load sequences for this chromosome
            sequences = {}
            try:
                for record in SeqIO.parse(fasta_file, 'fasta'):
                    sequences[record.id] = str(record.seq)
                cache[chromosome] = sequences
            except Exception as e:
                print(f"Warning: Could not load {fasta_file}: {e}")
                cache[chromosome] = {}
                return None
        
        # Return sequence if found in cache
        return cache[chromosome].get(transcript_id)
    
    return get_sequence


def extract_transcript_id(attributes):
    """Extract transcript ID from GFF attributes column."""
    # Handle different attribute formats
    if 'transcript_id' in attributes:
        # GTF format
        parts = attributes.split('transcript_id "')
        if len(parts) > 1:
            tid = parts[1].split('"')[0]
            return tid
    elif 'ID=' in attributes:
        # GFF3 format
        parts = attributes.split('ID=')
        if len(parts) > 1:
            tid = parts[1].split(';')[0]
            return tid
    elif ';' in attributes:
        # Simple format (as used in relocate_transcripts.py)
        parts = attributes.split(';')
        if len(parts) >= 2:
            tid = parts[1].strip()
            return tid
    
    return attributes.strip()


def extract_stop_codons(gff_file, gffread_dir, output_file=None):
    """
    Extract stop codons from transcripts.
    
    Args:
        gff_file: Path to GFF file with relocated CDS coordinates
        gffread_dir: Path to directory containing chromosome-specific FASTA files
        output_file: Optional file to write detailed results
    
    Returns:
        Dictionary with stop codon statistics
    """
    # Load data
    gff_df = load_gff(gff_file)
    get_sequence = load_sequences_from_directory(gffread_dir)
    
    # Filter for CDS features
    cds_df = gff_df[gff_df['Feature'] == 'CDS'].copy()
    
    # Group by transcript and find the last CDS for each
    stop_codons = []
    stop_codon_counts = Counter()
    missing_transcripts = set()
    
    for transcript_id, group in cds_df.groupby('Attributes'):
        # Sort by End position to find the last CDS
        group_sorted = group.sort_values('End')
        last_cds = group_sorted.iloc[-1]
        
        # Extract transcript ID for FASTA lookup
        tid = extract_transcript_id(str(transcript_id))
        chromosome = last_cds['Chromosome']
        
        # Get sequence from appropriate chromosome file
        seq= get_sequence(tid, chromosome)
        
        # Check if sequence exists
        if seq is None:
            if tid not in missing_transcripts:
                missing_transcripts.add(tid)
                print(f"Warning: Transcript {tid} not found in {chromosome} FASTA file")
            continue
        
        # Get the end position of the last CDS (1-based in GFF, but we need 0-based for Python)
        # The stop codon starts right after the last CDS ends
        stop = int(last_cds['End'])  # This is 1-based, pointing to last nucleotide of CDS
        
        # Extract stop codon (positions cds_end to cds_end+2, 0-based indexing)
        try:
            # Check for stop codon within +/-1 of annotated stop
            actual_stop = None
            actual_stop_sequence = None
            
            # Check at position stop-1
            if stop-1 >= 0 and stop+2 <= len(seq):
                candidate = seq[stop-1:stop+2]
                if candidate.upper() in ('TGA', 'TAG', 'TAA'):
                    actual_stop = stop - 1
                    actual_stop_sequence = candidate
            
            # Check at position stop (default)
            if actual_stop is None and stop >= 0 and stop+3 <= len(seq):
                candidate = seq[stop:stop+3]
                if candidate.upper() in ('TGA', 'TAG', 'TAA'):
                    actual_stop = stop
                    actual_stop_sequence = candidate
            
            # Check at position stop+1
            if actual_stop is None and stop+1 >= 0 and stop+4 <= len(seq):
                candidate = seq[stop+1:stop+4]
                if candidate.upper() in ('TGA', 'TAG', 'TAA'):
                    actual_stop = stop + 1
                    actual_stop_sequence = candidate
            
            if actual_stop is None or actual_stop_sequence is None:
                print(f"Warning: {tid} stop codon not found within +/-1 at stop {stop}")
                continue  # Skip if stop codon not found within +/-1
            
            # Validate that we got 3 nucleotides
            if len(actual_stop_sequence) < 3:
                print(f"Warning: Could not extract full stop codon for {tid} (got {len(actual_stop_sequence)} bp)")
                continue
            
            stop_codon_upper = actual_stop_sequence.upper()
            stop_codons.append({
                'transcript_id': tid,
                'stop_codon': stop_codon_upper,
                'chromosome': chromosome,
                'strand': last_cds['Strand'],
                'cds_end_position': stop
            })
            stop_codon_counts[stop_codon_upper] += 1
            
        except IndexError:
            print(f"Warning: Sequence too short for {tid} (length: {len(seq)}, trying to access around position {stop})")
            continue
    
    # Print results
    print("\n" + "="*60)
    print("STOP CODON DISTRIBUTION")
    print("="*60)
    print(f"Total transcripts analyzed: {len(stop_codons)}")
    print(f"Unique stop codons found: {len(stop_codon_counts)}")
    print()
    
    # Sort by frequency
    sorted_codons = sorted(stop_codon_counts.items(), key=lambda x: x[1], reverse=True)
    
    print("Stop Codon Frequencies:")
    print("-" * 60)
    total_count = sum(count for _, count in sorted_codons)
    for codon, count in sorted_codons:
        percentage = (count / total_count) * 100
        print(f"{codon}: {count:>6} ({percentage:>6.2f}%)")
    
    print("="*60)
    
    # Write detailed output if requested
    if output_file:
        with open(output_file, 'w') as f:
            f.write("transcript_id\tstop_codon\tchromosome\tstrand\tcds_end_position\n")
            for record in stop_codons:
                f.write(f"{record['transcript_id']}\t{record['stop_codon']}\t"
                       f"{record['chromosome']}\t{record['strand']}\t{record['cds_end_position']}\n")
        print(f"\nDetailed results written to: {output_file}")
    
    return {
        'stop_codons': stop_codons,
        'counts': stop_codon_counts,
        'total': len(stop_codons)
    }


def main():
    parser = argparse.ArgumentParser(
        description='Extract and analyze stop codons from transcripts'
    )
    parser.add_argument(
        'gff_file',
        help='GFF file with relocated CDS coordinates (from relocate_transcripts.py)'
    )
    parser.add_argument(
        'gffread_dir',
        help='Directory containing chromosome-specific FASTA files (from gffread)'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output file for detailed stop codon information (optional)',
        default=None
    )
    
    args = parser.parse_args()
    
    extract_stop_codons(args.gff_file, args.gffread_dir, args.output)


if __name__ == '__main__':
    main()
