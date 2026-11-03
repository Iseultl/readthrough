#!/usr/bin/env python3
"""
Extract full protein sequences (CDS start -> stop, in-frame TGA as
selenocysteine U) for one species of the multi-species target set defined by
a windowed .faa file (e.g. /no_backup/rg/ileahy/protists/SECIS_2species.faa).

Target set: one record per target transcript, header
    <prefix>|<species>|<transcript_name>|sec_start=..|codon=TGA|score_diff=..|cds_start=..|...
The record sequence is the protein of the 603-nt window around the U site
(frame-0 translation, TGA -> U, to the next in-frame stop) and is used only as
a cross-check of the extracted protein.

CDS boundaries (re_score_start..re_score_end, 1-based inclusive, on the gffread
cDNA transcript; re_score_end already excludes the stop codon) and the TGA site
(TGA_site_score, 0-based) are read from the species'
<output_dir>/<species>_ORFsearch_SECIS.result.

The transcript FASTA must contain only the target transcripts: the driver
(bin/run_extract_faa_full_proteins.sh) filters each species' Lyric GFF down to
the target transcript_ids before gffread, so this script never touches the genome.

Extraction/translation is identical to bin/extract_positive_proteins.py:
in-frame TGA within +/-1 of the annotated position, exactly one in-frame stop
in the CDS (the TGA itself), standard-code translation, X kept where the genome
has N gaps, TGA replaced by U.

Window cross-check: the .faa window protein is aligned to the full protein by
the U position (clipped to the overlap); an X in the full protein (genome N
gap) matches any residue. Result per record: exact | xgap | none (logged,
never dropped).

Outputs (in --outdir):
    <species>_full_proteins.fa        one record per target transcript
    <species>_extraction_summary.tsv  per-protein details
    <species>_extraction_log.tsv      skipped transcripts + reasons
"""

import argparse
import csv
import os
import sys

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


def parse_faa(path):
    """Return {transcript_name: {prefix, species, tid, kv, winprot}} for the
    whole file (all species)."""
    recs = {}
    cur = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                parts = line[1:].split("|")
                kv = dict(x.split("=", 1) for x in parts[3:] if "=" in x)
                recs[parts[2]] = {
                    "prefix": parts[0],
                    "species": parts[1],
                    "tid": parts[2],
                    "kv": kv,
                    "winprot": "",
                }
                cur = parts[2]
            elif cur is not None:
                recs[cur]["winprot"] += line
    return recs


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


def _num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def load_coords(result, tids):
    """Stream the result CSV (it can be large and a transcript may fan out to
    several rows via the SECIS left-join); keep the first row per target tid.
    Returns {tid: dict(re_start, re_end, tga_0, re_score, og_score)}."""
    want = set(tids)
    found = {}
    with open(result) as f:
        reader = csv.reader(f)
        header = next(reader)
        if "transcript_name" not in header:
            header[0] = "transcript_name"
        i = {c: k for k, c in enumerate(header)}
        for row in reader:
            tid = row[i["transcript_name"]]
            if tid in want and tid not in found:
                found[tid] = {
                    "re_start": _num(row[i["re_score_start"]]),
                    "re_end": _num(row[i["re_score_end"]]),
                    "tga_0": _num(row[i["TGA_site_score"]]),
                    "re_score": _num(row[i["re_score_score"]]),
                    "og_score": _num(row[i["og_score_score"]]),
                }
                if len(found) == len(want):
                    break
    return found


def window_match(full, win):
    """Align the .faa window protein to the full protein by the U position,
    clipped to the overlap. An X in the full protein (genome N gap) matches
    any residue. Returns 'exact', 'xgap' or 'none'."""
    if not win or "U" not in win or "U" not in full:
        return "none"
    ui, uidx = win.index("U"), full.index("U")
    exact = True
    any_cmp = False
    for i, w in enumerate(win):
        j = uidx + (i - ui)
        if 0 <= j < len(full):
            any_cmp = True
            if w != full[j]:
                if full[j] == "X":
                    exact = False
                else:
                    return "none"
    if not any_cmp:
        return "none"
    return "exact" if exact else "xgap"


SUMMARY_COLS = ["species", "faa_prefix", "transcript_name", "protein_length",
                "U_position", "tga_1based", "tga_offset_from_annotation",
                "re_score_start", "re_score_end", "re_score_score",
                "og_score_score", "score_diff", "window_match", "warning"]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--faa", required=True,
                    help="multi-species window .faa (target list + window proteins)")
    ap.add_argument("--species", required=True)
    ap.add_argument("--result", required=True,
                    help="<species>_ORFsearch_SECIS.result")
    ap.add_argument("--transcripts", required=True,
                    help="gffread transcript FASTA of the filtered target transcripts")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    fa_out = os.path.join(args.outdir, f"{args.species}_full_proteins.fa")
    sum_out = os.path.join(args.outdir, f"{args.species}_extraction_summary.tsv")
    log_out = os.path.join(args.outdir, f"{args.species}_extraction_log.tsv")

    recs = {t: r for t, r in parse_faa(args.faa).items()
            if r["species"] == args.species}
    if not recs:
        print(f"[{args.species}] no target records in {args.faa}", file=sys.stderr)
        open(fa_out, "w").close()
        with open(sum_out, "w") as f:
            f.write("\t".join(SUMMARY_COLS) + "\n")
        with open(log_out, "w") as f:
            f.write("transcript_name\treason\n")
        return

    coords = load_coords(args.result, set(recs))

    # gffread headers are the transcript ids (with the >63-char shortening
    # fallback, same as extract_positive_proteins.py)
    seqs = read_fasta(args.transcripts)
    lookup = {}
    for tid, s in seqs.items():
        lookup.setdefault(tid, s)
        if len(tid) > 63:
            lookup.setdefault(shorten_id(tid), s)

    records = []
    summary = []
    log = []

    for tid, rec in recs.items():
        c = coords.get(tid)
        if c is None:
            log.append((tid, "transcript not in result file"))
            continue
        if c["re_start"] is None or c["re_end"] is None or c["tga_0"] is None:
            log.append((tid, "NA CDS coordinates in result file"))
            continue
        re_start = int(c["re_start"])
        re_end = int(c["re_end"])
        tga_0 = int(c["tga_0"])

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

        prot = str(Seq(cds).translate())
        u_pos = (p - re_start) // 3
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

        wmatch = window_match(prot, rec["winprot"])
        header = (f"{rec['prefix']}|species={args.species}|transcript={tid}"
                  f"|U_pos={u_pos + 1}|re_score={c['re_score']:.2f}"
                  f"|og_score={c['og_score']:.2f}|window={wmatch}")
        records.append((header, prot))

        summary.append({
            "species": args.species,
            "faa_prefix": rec["prefix"],
            "transcript_name": tid,
            "protein_length": len(prot),
            "U_position": u_pos + 1,
            "tga_1based": p,
            "tga_offset_from_annotation": p - (tga_0 + 1),
            "re_score_start": re_start,
            "re_score_end": re_end,
            "re_score_score": c["re_score"],
            "og_score_score": c["og_score"],
            "score_diff": (c["re_score"] - c["og_score"]
                           if c["re_score"] is not None and c["og_score"] is not None
                           else ""),
            "window_match": wmatch,
            "warning": warning,
        })

    with open(fa_out, "w") as f:
        for header, prot in records:
            f.write(f">{header}\n{prot}\n")

    with open(sum_out, "w") as f:
        f.write("\t".join(SUMMARY_COLS) + "\n")
        for row in summary:
            f.write("\t".join(str(row[c]) for c in SUMMARY_COLS) + "\n")
    with open(log_out, "w") as f:
        f.write("transcript_name\treason\n")
        for tid, reason in log:
            f.write(f"{tid}\t{reason}\n")

    print(f"[{args.species}] {len(recs)} targets -> {len(records)} proteins, "
          f"{len(log)} skipped -> {fa_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
