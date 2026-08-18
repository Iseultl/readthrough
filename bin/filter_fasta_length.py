#!/usr/bin/env python3

import argparse
from pathlib import Path
from Bio import SeqIO


parser = argparse.ArgumentParser(
    description=(
        "Remove FASTA sequences longer than a specified length and "
        "duplicate FASTA identifiers. Optionally split the output into "
        "multiple FASTA files."
    )
)

parser.add_argument(
    "-i", "--input",
    required=True,
    help="Input FASTA file"
)

parser.add_argument(
    "-o", "--output",
    required=True,
    help="Output FASTA file or base name for split output files"
)

parser.add_argument(
    "-m", "--max_length",
    type=int,
    default=4000,
    help="Maximum allowed sequence length (default: 4000)"
)

parser.add_argument(
    "-s", "--split_size",
    type=int,
    default=1000,
    help="Maximum number of sequences per output FASTA file (default: 1000)"
)

args = parser.parse_args()


seen_ids = set()

total = 0
removed_length = 0
removed_duplicates = 0
kept_records = []


# Filter sequences
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


# Split output into chunks
output_path = Path(args.output)

num_files = 0

for start in range(0, len(kept_records), args.split_size):

    chunk = kept_records[start:start + args.split_size]

    part_number = (start // args.split_size) + 1

    output_file = (
        output_path.parent /
        f"{output_path.stem}_part_{part_number:03d}{output_path.suffix}"
    )

    SeqIO.write(chunk, output_file, "fasta")

    num_files += 1

    print(
        f"Wrote {len(chunk)} sequences to: {output_file}"
    )


# Summary
print("\nSummary")
print(f"Total sequences: {total}")
print(f"Kept {len(kept_records)} sequences.")
print(f"Removed {removed_length} sequences (> {args.max_length} nt).")
print(f"Removed {removed_duplicates} sequences with duplicate FASTA identifiers.")
print(f"Created {num_files} output FASTA file(s).")