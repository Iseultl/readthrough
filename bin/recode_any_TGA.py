#!/usr/bin/env python3

# The purpose of this script is to take rna 
# transcripts as input 
# and recode any instance of a TGA to a TGC, 
# producing a new transcript for each recoding instance

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Identify positions of TGA codons in each transcript sequence
def find_tga_occurrences(seq_str):
    positions = []
    start = 0

    while True:
        pos = seq_str.find("TGA", start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1

    return positions


# Perform recoding without holding all records in memory
def recode(fasta_file, recodon, output_fasta):
    with open(output_fasta, "w") as out_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_str = str(record.seq)

            if len(seq_str) > 8000:
                # Skip sequences longer than 8000 bases
                continue

            original_record = SeqRecord(Seq(seq_str), id=f"{record.id}_original", description="")
            SeqIO.write(original_record, out_fasta, "fasta")

            for pos in find_tga_occurrences(seq_str):
                new_seq = seq_str[:pos] + recodon + seq_str[pos + 3:]
                record_id = f"{record.id}_TGA_{pos}_recode"
                recoded_record = SeqRecord(Seq(new_seq), id=record_id, description="")
                SeqIO.write(recoded_record, out_fasta, "fasta")

    print(f"Rewritten sequences saved to {output_fasta}")
        
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Recode TGA codons at selenocysteine positions in a genome FASTA file")
    parser.add_argument('--fasta', type=str, required=True, help="Path to transcript FASTA file") 
    parser.add_argument('--recodon', type=str, required=True, choices=['TGT', 'TGC'], help="Codon to replace TGA (TGT or TGC)")
    parser.add_argument('--output', type=str, required=True, help="Output FASTA file of recoded sequences")
     
    args = parser.parse_args()
    
    recode(args.fasta, args.recodon, args.output)  # Recode sequences and save



    
# Command for running script 
'''
python recode_any_TGA.py --fasta gffread_out/transcripts_clean.fa --gff subset_mouse.all_secis.gff --recodon TGC --output recoded_mouse_test.fa
'''  