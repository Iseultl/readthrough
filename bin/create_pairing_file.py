#!/usr/bin/env python3

from pathlib import Path
import argparse
import sys

parser = argparse.ArgumentParser(description="Create a pairing file for chromosome FASTA and GTF files")
parser.add_argument('--fasta', type=str, required=True, help="Path to directory containing FASTA files")
parser.add_argument('--gtf', type=str, required=True, help="Path to directory containing GTF files")
parser.add_argument('--output', type=str, required=True, help="Path to output file")
args = parser.parse_args()

fasta_files = {f.name: f for f in Path(args.fasta).glob('*.fa')}
gtf_files = {f.name: f for f in Path(args.gtf).glob('*.gtf')}
print(gtf_files)

paired_count = 0
unmatched_fasta = []
unmatched_gtf = set(gtf_files.keys())

with open(args.output, 'w') as out:
    for chrom in gtf_files:
        chrom = chrom.replace(".gtf", "")
        fasta_name = "horse_genome.part_" + chrom + ".fa"
        if fasta_name in fasta_files:
            out.write(f"{chrom}\t{gtf_files[chrom+'.gtf']}\t{fasta_files[fasta_name]}\n")
            paired_count += 1
            unmatched_gtf.discard(chrom)
        else:
            unmatched_fasta.append(chrom)

# Report any issues
fasta_count = 0
if unmatched_fasta:
    print(f"[WARNING] No matching GTF found for the following FASTA files:", file=sys.stderr)
    fasta_count = fasta_count + 1

gtf_count = 0
if unmatched_gtf:
    print(f"[WARNING] No matching FASTA found for the following GTF files:", file=sys.stderr)
    gtf_count = gtf_count + 1


if paired_count == 0:
    print("[ERROR] No chromosome pairs were created. Please check your FASTA and GTF directories.", file=sys.stderr)
    sys.exit(1)

print(f"[INFO] Successfully created {paired_count} chromosome pair(s).")
