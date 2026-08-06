#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(
    description="Remove FASTA sequences longer than a specified length."
)
parser.add_argument(
    "--input",
    "-i",
    required=True,
    help="Input FASTA file"
)
parser.add_argument(
    "--output",
    "-o",
    required=True,
    help="Output FASTA file"
)
parser.add_argument(
    "--max_length",
    "-m",
    type=int,
    default=4000,
    help="Maximum allowed sequence length (default: 4000)"
)

args = parser.parse_args()


def write_record(header, seq, out_handle):
    """Write a FASTA record if it meets the length threshold."""
    if header is not None and len(seq) <= args.max_length:
        out_handle.write(header)
        for i in range(0, len(seq), 60):
            out_handle.write(seq[i:i+60] + "\n")


kept = 0
removed = 0

with open(args.input) as infile, open(args.output, "w") as outfile:

    header = None
    sequence = []

    for line in infile:
        if line.startswith(">"):
            if header is not None:
                seq = "".join(sequence)
                if len(seq) <= args.max_length:
                    kept += 1
                    write_record(header, seq, outfile)
                else:
                    removed += 1

            header = line
            sequence = []

        else:
            sequence.append(line.strip())

    # Write final record
    if header is not None:
        seq = "".join(sequence)
        if len(seq) <= args.max_length:
            kept += 1
            write_record(header, seq, outfile)
        else:
            removed += 1

print(f"Kept {kept} sequences.")
print(f"Removed {removed} sequences (> {args.max_length} nt).")