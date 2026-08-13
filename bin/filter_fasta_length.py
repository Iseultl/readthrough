#!/usr/bin/env python3

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="Remove FASTA sequences longer than a specified length and duplicate FASTA identifiers."
)

parser.add_argument("-i", "--input", required=True,
                    help="Input FASTA file")
parser.add_argument("-o", "--output", required=True,
                    help="Output FASTA file")
parser.add_argument("-m", "--max_length", type=int, default=4000,
                    help="Maximum allowed sequence length (default: 4000)")

args = parser.parse_args()

seen_ids = set()

total = 0
removed_length = 0
removed_duplicates = 0
kept_records = []

for record in SeqIO.parse(args.input, "fasta"):
    total += 1

    # Remove sequences longer than the maximum length
    if len(record.seq) > args.max_length:
        removed_length += 1
        continue

    # Remove duplicate FASTA identifiers
    if record.id in seen_ids:
        removed_duplicates += 1
        continue

    seen_ids.add(record.id)
    kept_records.append(record)

kept = SeqIO.write(kept_records, args.output, "fasta")

print(f"Total sequences: {total}")
print(f"Kept {kept} sequences.")
print(f"Removed {removed_length} sequences (> {args.max_length} nt).")
print(f"Removed {removed_duplicates} sequences with duplicate FASTA identifiers.")