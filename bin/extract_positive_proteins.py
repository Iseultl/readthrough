#!/usr/bin/env python3
"""
Extract positive, SECIS-containing protein sequences for one species.

Selection (identical criteria to bin/extract_sequence_logos.py):
    positive = (re_score_score - og_score_score) >= -1.8
                AND re_score_score >= -1.5
plus the additional filter:
    the transcript has at least one identified SECIS element
    (all_secis_start is not NA in the final ORFsearch+SECIS result).

For each selected transcript the recoded CDS
(re_score_start..re_score_end, 1-based inclusive, on the gffread cDNA
transcript; re_score_end already excludes the stop codon) is translated
with the standard code and the in-frame TGA at TGA_site_score (0-based,
as stored in the result table) is replaced by selenocysteine (U).

Outputs (in --outdir):
    <species>_positive_proteins.fa      one record per selected transcript
    <species>_extraction_summary.tsv    per-protein details
    <species>_extraction_log.tsv        skipped transcripts + reasons
"""

import argparse
import os
import sys

import pandas as pd
from Bio.Seq import Seq


def shorten_id(tid):
    """Replicate the GFFREAD_CHR header shortening (commit 760ab22):
    ids > 63 chars become first 55 chars + '_' + 6-digit rolling hash
    h = (h * 31 + ord(c)) % 1000000 over the FULL id."""
    if len(tid) <= 63:
        return tid
    h = 0
    for c in tid:
        h = (h * 31 + ord(c)) % 1000000
    return tid[:55] + "_" + "%06d" % h


def read_fasta(path):
    seqs = {}
    with open(path) as f:
        cur = None
        parts = []
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(parts)
                cur = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if cur is not None:
            seqs[cur] = "".join(parts)
    return seqs


def find_tga(seq, tga_0based, re_start, re_end):
    """Locate the TGA codon, 1-based, allowing the same +/-1 nt tolerance
    as infer_sec_start() in extract_sequence_logos.py. Must be in the CDS
    and in frame with re_score_start."""
    cands = []
    for off in (0, 1, -1):
        p = tga_0based + 1 + off
        if re_start <= p <= re_end - 2 and (p - re_start) % 3 == 0 \
                and seq[p - 1:p + 2] == "TGA":
            cands.append((abs(off), p))
    if not cands:
        return None
    cands.sort()
    return cands[0][1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--result", required=True,
                    help="Final result file <sp>_ORFsearch_SECIS.result")
    ap.add_argument("--transcripts",
                    help="gffread transcript FASTA for the species genome "
                         "(omit with --select-only)")
    ap.add_argument("--species", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--select-only", action="store_true",
                    help="only write <species>_targets.tsv (and empty "
                         "outputs if there are no targets); do not read "
                         "transcripts or translate")
    args = ap.parse_args()
    if not args.select_only:
        if not args.transcripts:
            ap.error("--transcripts is required unless --select-only is given")

    os.makedirs(args.outdir, exist_ok=True)
    fa_out = os.path.join(args.outdir, f"{args.species}_positive_proteins.fa")
    sum_out = os.path.join(args.outdir, f"{args.species}_extraction_summary.tsv")
    log_out = os.path.join(args.outdir, f"{args.species}_extraction_log.tsv")

    df = pd.read_csv(args.result, index_col=0)

    # ---- select positives with an identified SECIS element ----
    for col in ("re_score_start", "re_score_end", "re_score_score",
                "og_score_score", "TGA_site_score"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    diff = df["re_score_score"] - df["og_score_score"]
    positive = (diff >= -1.8) & (df["re_score_score"] >= -1.5)
    n_pos_total = int(positive.sum())

    sel = df[positive]
    sel = sel[sel["all_secis_start"].notna()]
    sel = sel[sel["re_score_start"].notna() & sel["re_score_end"].notna()
              & sel["TGA_site_score"].notna()]
    # the SECIS left-join fans a transcript out to one row per hit
    sel = sel[~sel.index.duplicated(keep="first")]
    print(f"[{args.species}] positives={n_pos_total}, "
          f"positives with SECIS={len(sel)}", file=sys.stderr)

    summary_cols = ["species", "transcript_name", "gene_name",
                    "protein_length", "U_position", "tga_1based",
                    "tga_offset_from_annotation", "re_score_start",
                    "re_score_end", "re_score_score", "og_score_score",
                    "score_diff", "all_secis_start", "all_secis_end",
                    "warning"]

    if args.select_only:
        targets_out = os.path.join(args.outdir, f"{args.species}_targets.tsv")
        with open(targets_out, "w") as f:
            f.write("transcript_name\n")
            for tid in sel.index:
                f.write(f"{tid}\n")
        if sel.empty:
            # write empty outputs so downstream steps stay uniform
            open(fa_out, "w").close()
            pd.DataFrame(columns=summary_cols).to_csv(sum_out, sep="\t",
                                                      index=False)
            pd.DataFrame(columns=["transcript_name", "reason"]).to_csv(
                log_out, sep="\t", index=False)
            print(f"[{args.species}] no targets; empty outputs written",
                  file=sys.stderr)
        else:
            print(f"[{args.species}] {len(sel)} targets -> {targets_out}",
                  file=sys.stderr)
        return

    # ---- load gffread transcripts, with shortened-id fallback ----
    seqs = read_fasta(args.transcripts)
    lookup = {}
    for tid, s in seqs.items():
        lookup.setdefault(tid, s)
        if len(tid) > 63:
            lookup.setdefault(shorten_id(tid), s)

    summary = []
    log = []
    records = []

    for tid, row in sel.iterrows():
        re_start = int(row["re_score_start"])
        re_end = int(row["re_score_end"])
        tga_0 = int(row["TGA_site_score"])

        seq = lookup.get(tid)
        if seq is None:
            log.append((tid, "transcript not found in gffread output"))
            continue
        if re_end > len(seq):
            log.append((tid, f"re_score_end {re_end} beyond transcript length {len(seq)}"))
            continue

        cds = seq[re_start - 1:re_end]
        if len(cds) % 3 != 0:
            log.append((tid, f"CDS length {len(cds)} not a multiple of 3"))
            continue

        p = find_tga(seq, tga_0, re_start, re_end)
        if p is None:
            log.append((tid, f"TGA not found in-frame within +/-1 of annotated position {tga_0 + 1}"))
            continue

        u_pos = (p - re_start) // 3
        prot = str(Seq(cds).translate())
        # The only in-frame stop allowed in the original CDS is the TGA
        # itself (every other in-frame stop would have terminated the
        # recoded ORF geneid predicted, since only this TGA was recoded).
        if prot[u_pos] != "*" or prot.count("*") != 1:
            log.append((tid, f"in-frame stop codon(s) in CDS {re_start}-{re_end} "
                              f"other than the TGA at position {p}"))
            continue
        warning = ""
        if prot[0] != "M":
            warning = f"does not start with M (starts {prot[0]})"
        prot = prot[:u_pos] + "U" + prot[u_pos + 1:]

        gene = "" if pd.isna(row.get("gene_name")) else str(row["gene_name"])
        header = (f"{tid}|species={args.species}|gene={gene}"
                  f"|U_pos={u_pos + 1}|re_score={row['re_score_score']:.2f}"
                  f"|og_score={row['og_score_score']:.2f}"
                  f"|secis={int(row['all_secis_start'])}-{int(row['all_secis_end'])}")
        records.append((header, prot))

        summary.append({
            "species": args.species,
            "transcript_name": tid,
            "gene_name": gene,
            "protein_length": len(prot),
            "U_position": u_pos + 1,
            "tga_1based": p,
            "tga_offset_from_annotation": p - (tga_0 + 1),
            "re_score_start": re_start,
            "re_score_end": re_end,
            "re_score_score": row["re_score_score"],
            "og_score_score": row["og_score_score"],
            "score_diff": row["re_score_score"] - row["og_score_score"],
            "all_secis_start": int(row["all_secis_start"]),
            "all_secis_end": int(row["all_secis_end"]),
            "warning": warning,
        })

    with open(fa_out, "w") as f:
        for header, prot in records:
            f.write(f">{header}\n{prot}\n")

    pd.DataFrame(summary).to_csv(sum_out, sep="\t", index=False)
    pd.DataFrame(log, columns=["transcript_name", "reason"]).to_csv(
        log_out, sep="\t", index=False)

    print(f"[{args.species}] wrote {len(records)} proteins -> {fa_out} "
          f"({len(log)} skipped)", file=sys.stderr)


if __name__ == "__main__":
    main()
