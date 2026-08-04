#!/usr/bin/env python3

# The purpose of this script is to take rna 
# transcripts as input 
# and recode any instance of a TGA to a TGC, 
# producing a new transcript for each recoding instance

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Load FASTA sequences into a dictionary
def load_fasta(fasta_file):
    return {record.id: record.seq for record in SeqIO.parse(fasta_file, 'fasta')}  

# Identify positions of TGA codons in each transcript sequence
def find_tga_occurrences(sequences):
    tga_positions = {}
    
    for transcript_name, seq in sequences.items():
        seq_str = str(seq)  # Convert Seq object to string
        positions = [pos for pos in range(len(seq_str) - 2) if seq_str[pos:pos+3] == "TGA"]
        
        if positions:
            tga_positions[transcript_name] = positions  # Map transcript name to TGA positions
            # print("Found TGAs for ", transcript_name)
        else:
            continue
    
    return tga_positions
          
# Perform recoding
def recode(sequences, recodon, tga_positions, output_fasta):
    records = []
    
    for transcript_name, seq in sequences.items():
        
        seq_str = str(seq)  # Convert Seq object to string
        positions = tga_positions.get(transcript_name, [])  # Get TGA positions or empty list
        if len(seq_str) > 4000:
            # Skip sequences longer than 4000 bases
            continue
        # Store the original sequence first
        original_record = SeqRecord(Seq(seq_str), id=f"{transcript_name}_original", description="")
        records.append(original_record)
        
        for pos in positions:
            # Create a new sequence with TGA replaced by the recodon at position `pos`
            new_seq = seq_str[:pos] + recodon + seq_str[pos+3:]  
            # Create a new record for FASTA output
            record_id = f"{transcript_name}_TGA_{pos}_recode"
            record = SeqRecord(Seq(new_seq), id=record_id, description="")
            records.append(record)

    # Write all modified sequences to output FASTA
    with open(output_fasta, "w") as out_fasta:
        SeqIO.write(records, out_fasta, "fasta")

    print(f"Rewritten sequences saved to {output_fasta}")
        
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Recode TGA codons at selenocysteine positions in a genome FASTA file")
    parser.add_argument('--fasta', type=str, required=True, help="Path to transcript FASTA file") 
    parser.add_argument('--recodon', type=str, required=True, choices=['TGT', 'TGC'], help="Codon to replace TGA (TGT or TGC)")
    parser.add_argument('--output', type=str, required=True, help="Output FASTA file of recoded sequences")
     
    args = parser.parse_args()
    
    fasta = load_fasta(args.fasta)  # Load sequences
    
    TGAs = find_tga_occurrences(fasta)  # Find TGA positions
    recode(fasta, args.recodon, TGAs, args.output)  # Recode sequences and save



    
# Command for running script 
'''
python recode_any_TGA.py --fasta gffread_out/transcripts_clean.fa --gff subset_mouse.all_secis.gff --recodon TGC --output recoded_mouse_test.fa
'''  