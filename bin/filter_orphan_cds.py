#!/usr/bin/env python3

# Remove GFF3 CDS lines whose Parent does not correspond to a transcript
# (RNA / mRNA) record in the file. Such orphan CDS records (e.g. transdecoder
# ORF scans whose parents carry .pN suffixes but have no transcript line,
# present in the merged gidRef GFFs) make gffread materialize spurious
# CDS-only "transcripts". Chlorella_sorokiniana: 19,907 spurious sequences
# (46,367 gffread vs 26,460 real transcripts).
# Two-pass, streaming; comment lines and all non-CDS features pass through.

import argparse


def attr(field9, key):
    prefix = key + "="
    for a in field9.split(';'):
        if a.startswith(prefix):
            return a[len(prefix):]
    return None


def main():
    parser = argparse.ArgumentParser(description="Remove orphan CDS lines from a GFF3 file")
    parser.add_argument('--gff', required=True, help="Input GFF3 file")
    parser.add_argument('--output', required=True, help="Output GFF3 file")
    args = parser.parse_args()

    # Pass 1: collect transcript-level IDs
    transcript_ids = set()
    with open(args.gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 9 and parts[2] in ('RNA', 'mRNA'):
                tid = attr(parts[8], 'ID')
                if tid:
                    transcript_ids.add(tid)

    # Pass 2: drop CDS lines with an unknown Parent
    removed = kept = 0
    with open(args.gff) as f, open(args.output, 'w') as out:
        for line in f:
            if not line.startswith('#'):
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 9 and parts[2] == 'CDS':
                    parent = attr(parts[8], 'Parent')
                    if parent is not None and parent not in transcript_ids:
                        removed += 1
                        continue
            out.write(line)
            kept += 1

    print(f"filter_orphan_cds: {len(transcript_ids)} transcript IDs, "
          f"kept {kept} lines, removed {removed} orphan CDS lines -> {args.output}")


if __name__ == '__main__':
    main()
