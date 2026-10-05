#!/usr/bin/env python3
"""Pre-flight checks for new protist species in protist_filepaths.csv.

For every CSV row whose species has no runs/params_<sp>.yaml yet (i.e. a newly
added species), checks in a single pass over the (gzipped) inputs:

  1. all 4 input files exist and are non-empty
  2. every .gz decompresses to TEXT with a SINGLE gunzip
     (double-gz = gzip magic 1f 8b visible after one zcat -> the pipeline's
      UNZIP stage would fail; fix: zcat|zcat|gzip in place, keep a backup)
  3. lyric transcript count (whitelist of transcript-level feature types,
     same as the submit scripts)
  4. per-scaffold transcript presence from the LyRic GFF + genome scaffold
     set -> "transcriptless scaffold" count (scaffolds with GFF lines but no
     transcript-level line -> 0-byte transcripts_clean files downstream)
  5. SECISSEARCH_BATCH empty-bucket risk: with scaffold_batches=24 hash
     buckets (Java String.hashCode, same expression as main_protists.nf
     batchByScaffold), flag buckets whose ONLY members are transcriptless
     scaffolds (such a bucket produces no *.gff -> missing-output failure,
     Giardia_duodenalis job 29112739, 2026-10-02)
  6. max transcript ID length (>63 chars -> Seblastian header limit; handled
     by the GFFREAD shortening fix 760ab22, reported for the record)
  7. orphan CDS lines (CDS Parent=<id> with no matching transcript-level ID;
     removed by FILTER_ORPHAN_CDS, reported for the record)
  8. genome scaffold count, genome size
  9. pre-existing /no_backup/rg/ileahy/<sp>/ORFsearch* dirs (job backs them up)

Read-only: never touches the reference files.
"""
import gzip
import os
import sys

REPO = "/users/rg/ileahy/git/gitlab/readthrough"
CSV = os.path.join(REPO, "protist_filepaths.csv")
OUTBASE = "/no_backup/rg/ileahy"
K = 24  # params.scaffold_batches


def java_hash(s: str) -> int:
    h = 0
    for ch in s:
        h = (h * 31 + ord(ch)) & 0xFFFFFFFF
    if h >= 0x80000000:
        h -= 0x100000000
    return h


def bucket(filename: str) -> int:
    # Groovy: ((f.name as String).hashCode() % K).abs(), where f.name is the
    # FULL base name (transcripts_clean_<scaffold>.fa), not the scaffold name.
    # abs(x % K) == abs(x) % K
    return abs(java_hash(filename)) % K


def gz_first_bytes(path: str, n: int = 16) -> bytes:
    with gzip.open(path, "rb") as f:
        return f.read(n)


def opens(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


TRANSCRIPT_TYPES = {
    "mRNA", "transcript", "RNA", "rRNA", "tRNA",
    "scRNA", "circRNA", "ncRNA", "lncRNA",
}


def parse_ids(line9: str):
    """Extract (transcript_id, parent) from GTF3 or GFF3 attribute string."""
    t_id, parent = None, None
    for attr in line9.split(";"):
        attr = attr.strip()
        if attr.startswith("transcript_id "):
            t_id = attr.split(None, 1)[1]
        elif attr.startswith("ID="):
            t_id = t_id or attr[3:]
        elif attr.startswith("Parent="):
            parent = parent or attr[7:]
    return t_id, parent


def main():
    rows = []
    with open(CSV) as f:
        lines = [l.rstrip("\n") for l in f if l.strip()]
    header = lines[0].split(",")
    for line in lines[1:]:
        cols = line.split(",")
        path = cols[1].strip()
        # Path = <species_dir>.<taxid>/<assembly> -> species = dir name w/o .taxid
        sp = os.path.basename(os.path.dirname(path))
        if "." in sp:
            sp = sp.rsplit(".", 1)[0]  # strip .taxid
        rows.append((sp, cols))

    # Optional explicit list: python3 preflight_new_species.py sp1 sp2 ...
    # (auto-detection flags Trichomonas_vaginalis_G3: the already-run species
    #  whose params file uses the legacy short name params_Trichomonas_vaginalis.yaml)
    if len(sys.argv) > 1:
        wanted = sys.argv[1:]
        todo = [(sp, cols) for sp, cols in rows if sp in wanted]
    else:
        todo = [(sp, cols) for sp, cols in rows
                if not os.path.exists(os.path.join(REPO, "runs", f"params_{sp}.yaml"))]
    if not todo:
        print("No new species found (all have params files).")
        return

    print(f"=== Pre-flight: {len(todo)} new species ===\n")
    ok_all = True
    for sp, cols in todo:
        path = cols[1].strip()
        param, lyric, genome, gffref = (cols[i].strip() for i in (2, 3, 4, 5))
        P, LY, GM, GR = (f"{path}/{x}" for x in (param, lyric, genome, gffref))
        print(f"--- {sp} ---")
        issues = []

        # 1) existence
        missing = [n for n, f in (("param", P), ("lyric", LY), ("genome", GM), ("gffref", GR))
                   if not (os.path.isfile(f) and os.path.getsize(f) > 0)]
        if missing:
            issues.append(f"MISSING/EMPTY: {', '.join(missing)}")
            print("  " + "; ".join(issues))
            print()
            ok_all = False
            continue

        # 2) double-gz check (the three .gz inputs)
        doublegz = [n for n, f in (("lyric", LY), ("genome", GM), ("gffref", GR))
                    if f.endswith(".gz") and gz_first_bytes(f)[:2] == b"\x1f\x8b"]
        if doublegz:
            issues.append(f"DOUBLE-GZ (gzip magic after 1 gunzip): {', '.join(doublegz)}")

        # 3-7) single pass over the lyric GFF
        lyric_count = 0
        n_lines = 0
        scaffold_any = set()      # scaffolds with any feature line
        scaffold_transcript = set()  # scaffolds with a transcript-level line
        ids = set()               # all transcript-level IDs
        max_id = (0, "")
        cds_parents = set()
        n_cds = 0
        with opens(LY) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                n_lines += 1
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                scaffold, ftype = parts[0], parts[2]
                scaffold_any.add(scaffold)
                t_id, parent = parse_ids(parts[8])
                if ftype in TRANSCRIPT_TYPES:
                    lyric_count += 1
                    scaffold_transcript.add(scaffold)
                    if t_id:
                        ids.add(t_id)
                        if len(t_id) > max_id[0]:
                            max_id = (len(t_id), t_id)
                elif ftype == "CDS":
                    n_cds += 1
                    if parent:
                        cds_parents.add(parent)
        orphan_cds = cds_parents - ids
        if lyric_count == 0:
            issues.append("NO transcript-level features in lyric GFF")
        if max_id[0] > 63:
            issues.append(f"transcript IDs >63 chars (max {max_id[0]}: {max_id[1][:40]}...) - shortened at GFFREAD (760ab22), check result names")

        # genome: scaffold set + size
        n_scaffolds = 0
        genome_bp = 0
        genome_scaffolds = set()
        with opens(GM) as f:
            for line in f:
                if line.startswith(">"):
                    n_scaffolds += 1
                    genome_scaffolds.add(line[1:].split()[0])
                else:
                    genome_bp += len(line.strip())

        # 4/5) transcriptless scaffolds + empty-bucket risk
        # A transcripts_clean file exists only for scaffolds in BOTH the
        # GFF (any line) and the genome FASTA (main_protists.nf combine by:0).
        paired = scaffold_any & genome_scaffolds
        transcriptless = [s for s in paired if s not in scaffold_transcript]
        # SECISSEARCH_BATCH (and SPLIT/RECODE) bucket by the per-scaffold
        # TRANSCRIPT file name, transcripts_clean_<scaffold>.fa (gffread_chr.nf).
        buckets_all = {}
        for s in paired:
            buckets_all.setdefault(bucket(f"transcripts_clean_{s}.fa"), []).append(s)
        risky = {b: mem for b, mem in buckets_all.items()
                 if all(s not in scaffold_transcript for s in mem)}
        only_transcriptless = len(transcriptless)
        gff_not_in_fasta = len(scaffold_any - genome_scaffolds)

        print(f"  inputs: all 4 present")
        print(f"  double-gz: {'NONE' if not doublegz else doublegz}")
        print(f"  lyric GFF: {n_lines} feature lines, {lyric_count} transcript-level "
              f"(gffref/geneidM not scanned: only used for README)")
        print(f"  genome: {n_scaffolds} scaffolds, {genome_bp:,} bp")
        print(f"  scaffolds in GFF but not in FASTA: {gff_not_in_fasta} "
              f"(dropped at the combine; never reach gffread)")
        print(f"  transcriptless paired scaffolds (GFF lines, no transcript feature): "
              f"{only_transcriptless}"
              + (f" -> BUCKETS FULLY TRANSCRIPTLESS: {sorted(risky.keys())}" if risky else ""))
        print(f"  max transcript id length: {max_id[0]}")
        print(f"  orphan CDS parents (dropped by FILTER_ORPHAN_CDS): "
              f"{len(orphan_cds)} (of {n_cds} CDS lines)")

        # 9) pre-existing outputs
        pre = sorted(d for d in os.listdir(OUTBASE)
                     if d.startswith(sp + "/") or d.startswith(sp + ".")
                     or d.startswith(f"{sp}_") or d == sp)
        pre = [d for d in pre if d.startswith(sp)]
        outdirs = [d for d in pre if os.path.isdir(os.path.join(OUTBASE, d))]
        print(f"  pre-existing /no_backup/rg/ileahy/{sp}*: "
              f"{outdirs if outdirs else 'none (first run)'}")

        if issues:
            ok_all = False
            for i in issues:
                print(f"  ** ISSUE: {i}")
        print()

    print("=== Pre-flight " + ("PASSED" if ok_all else "HAS ISSUES") + " ===")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
