#!/usr/bin/env python3

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="Remove FASTA sequences longer than a specified length."
)

parser.add_argument("-i", "--input", required=True,
                    help="Input FASTA file")
parser.add_argument("-o", "--output", required=True,
                    help="Output FASTA file")
parser.add_argument("-m", "--max_length", type=int, default=4000,
                    help="Maximum allowed sequence length (default: 4000)")

args = parser.parse_args()

records = SeqIO.parse(args.input, "fasta")

kept_records = (
    record for record in records
    if len(record.seq) <= args.max_length
)

kept = SeqIO.write(kept_records, args.output, "fasta")

# Determine how many were removed
total = sum(1 for _ in SeqIO.parse(args.input, "fasta"))
removed = total - kept

print(f"Kept {kept} sequences.")
print(f"Removed {removed} sequences (> {args.max_length} nt).")