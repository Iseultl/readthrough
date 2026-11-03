#!/usr/bin/env bash
#
# Extract full proteins (CDS start -> stop, in-frame TGA as U) for the
# multi-species target set defined by a windowed .faa file.
#
# Per species (species list = unique species in the .faa, first-seen order):
#   1) gunzip the original genome + Lyric GFF (paths from runs/params_*.yaml,
#      matched by species_name; never modify the references)
#   2) filter the GFF to the target transcript_ids (exact match of
#      transcript_id / ID / Parent attribute values) -- only 1-9 transcripts
#      per species, so no per-scaffold batching or AGAT chain is needed
#   3) gffread -F -w transcripts.fa -g <genome> <filtered.gff>
#   4) bin/extract_faa_full_proteins.py
#        translates the CDS (re_score_start..re_score_end from
#        <output_dir>/<sp>_ORFsearch_SECIS.result) with U at the in-frame TGA
#        and cross-checks the result against the .faa window protein
#
# Outputs (in --outdir, default /no_backup/rg/ileahy/protists/full_proteins):
#   <sp>_full_proteins.fa         per species
#   <sp>_extraction_summary.tsv   per protein
#   <sp>_extraction_log.tsv       skipped transcripts + reasons
#   all_species_full_proteins.fa  combined, species in .faa order
#   extraction_summary.tsv        combined summary
#   .failed_species               species that errored (empty = all OK)
#
# Usage:
#   run_extract_faa_full_proteins.sh [--faa FILE] [--outdir DIR] [--parallel N] [--tmp-root DIR]

set -euo pipefail

REPO="/users/rg/ileahy/git/gitlab/readthrough"
SIF_GFFREAD="${SIF_GFFREAD:-$HOME/singularities/gffread.sif}"
SIF_PYTHON="${SIF_PYTHON:-$HOME/singularities/python.sif}"
FAA="/no_backup/rg/ileahy/protists/SECIS_2species.faa"
OUTDIR="/no_backup/rg/ileahy/protists/full_proteins"
PAR=4
TMPROOT="/nfs/scratch01/rg/ileahy/faa_full_proteins_tmp"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --faa)      FAA="$2"; shift 2 ;;
    --outdir)   OUTDIR="$2"; shift 2 ;;
    --parallel) PAR="$2"; shift 2 ;;
    --tmp-root) TMPROOT="$2"; shift 2 ;;
    *)          echo "unknown option: $1" >&2; exit 1 ;;
  esac
done

[[ -s "$FAA" ]] || { echo "faa file missing or empty: $FAA" >&2; exit 1; }

mkdir -p "$OUTDIR" "$TMPROOT"
: > "$OUTDIR/.failed_species"

# species list from the .faa (first-seen order)
mapfile -t SPECIES < <(grep '^>' "$FAA" | cut -d'|' -f2 | awk '!seen[$0]++')

# species_name -> params file (one line per params file; children read this,
# since exported associative arrays do not survive `bash -c` reliably)
: > "$TMPROOT/.params_map.tsv"
for f in "$REPO"/runs/params_*.yaml; do
  printf '%s\t%s\n' "$(awk -F'"' '/^species_name:/{print $2}' "$f")" "$f" \
    >> "$TMPROOT/.params_map.tsv"
done

# ungz <in> <out> — gunzip to text, tolerating the double-gz case seen in the
# reference sets (second layer detected by the gzip magic bytes).
ungz() {
  if [[ "$1" == *.gz ]]; then
    zcat -- "$1" > "$2"
    if head -c 2 "$2" | od -An -tx1 | tr -d ' \n' | grep -q '^1f8b'; then
      mv "$2" "$2.gz"
      zcat -- "$2.gz" > "$2"
    fi
  else
    cp -- "$1" "$2"
  fi
}

# filter_gff <gff> <tids.txt> <out.gff> — keep the lines whose
# transcript_id / ID / Parent attribute value is one of the target
# transcript ids (exact match, quoted or unquoted, comma lists allowed).
# Prints the number of kept lines.
filter_gff() {
  singularity exec -B "$TMPROOT" -B /no_backup "$SIF_PYTHON" python - "$1" "$2" "$3" <<'PYEOF'
import sys
gff, tidsf, out = sys.argv[1], sys.argv[2], sys.argv[3]
want = {l.strip() for l in open(tidsf) if l.strip()}

def attr_values(field9):
    vals = set()
    for part in field9.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
        elif " " in part:
            k, v = part.split(" ", 1)
        else:
            k, v = part, ""
        vals.add((k.strip(), v.strip().strip('"')))
    return vals

kept = 0
with open(gff) as f, open(out, "w") as o:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9:
            continue
        if any(k in ("transcript_id", "ID", "Parent") and
               (v in want or any(x in want for x in v.split(",")))
               for k, v in attr_values(p[8])):
            o.write(line)
            kept += 1
print(kept)
PYEOF
}

process_species() {
  local sp="$1"
  local start; start=$(date +%s)
  local params
  params=$(awk -F'\t' -v s="$sp" '$1==s{print $2; exit}' "$TMPROOT/.params_map.tsv")
  if [[ -z "$params" ]]; then
    echo "[$sp] ERROR: no runs/params_*.yaml with species_name $sp"
    echo "$sp" >> "$OUTDIR/.failed_species"
    return 0
  fi
  local outdir lyric genome result
  outdir=$(awk -F'"' '/^output_dir:/{print $2}' "$params")
  lyric=$(awk -F'"' '/^lyric_gtf:/{print $2}' "$params")
  genome=$(awk -F'"' '/^genome_fasta:/{print $2}' "$params")
  result="$outdir/${sp}_ORFsearch_SECIS.result"
  if [[ ! -s "$result" ]]; then
    echo "[$sp] ERROR: result file missing or empty ($result)"
    echo "$sp" >> "$OUTDIR/.failed_species"
    return 0
  fi

  local t="$TMPROOT/$sp"
  mkdir -p "$t"

  # 1) gunzip the original inputs (never modify the references)
  local gname lname
  gname=$(basename "$genome" .gz)
  lname=$(basename "$lyric" .gz)
  if ! ungz "$genome" "$t/$gname" || ! ungz "$lyric" "$t/$lname"; then
    echo "[$sp] ERROR: failed to decompress inputs ($genome / $lyric)"
    echo "$sp" >> "$OUTDIR/.failed_species"
    rm -rf "$t"
    return 0
  fi

  # 2) filter the GFF to the target transcripts of this species
  grep '^>' "$FAA" | awk -F'|' -v s="$sp" '$2==s {print $3}' > "$t/tids.txt"
  local ntargets
  ntargets=$(grep -c '' "$t/tids.txt")
  local nkept
  if ! nkept=$(filter_gff "$t/$lname" "$t/tids.txt" "$t/filtered.gff"); then
    echo "[$sp] ERROR: GFF filtering failed"
    echo "$sp" >> "$OUTDIR/.failed_species"
    rm -rf "$t"
    return 0
  fi
  if [[ "$nkept" -eq 0 ]]; then
    echo "[$sp] ERROR: no GFF lines carry the target transcript_ids"
    echo "$sp" >> "$OUTDIR/.failed_species"
    rm -rf "$t"
    return 0
  fi

  # 3) gffread the target transcripts only (1-9 per species: no OOM risk)
  if ! singularity exec -B "$t" -B /no_backup "$SIF_GFFREAD" \
        gffread -F -w "$t/transcripts.fa" -g "$t/$gname" "$t/filtered.gff" \
        > "$t/gffread.out" 2>&1; then
    echo "[$sp] ERROR: gffread failed"
    tail -n 20 "$t/gffread.out" 2>/dev/null
    echo "$sp" >> "$OUTDIR/.failed_species"
    rm -rf "$t"
    return 0
  fi
  rm -f "$t/$gname.fai"

  # 4) translate selected CDS, U at the TGA
  if ! singularity exec -B "$t" -B /no_backup -B "$REPO" -B "$OUTDIR" "$SIF_PYTHON" \
        python "$REPO/bin/extract_faa_full_proteins.py" \
        --faa "$FAA" --species "$sp" --result "$result" \
        --transcripts "$t/transcripts.fa" --outdir "$OUTDIR"; then
    echo "[$sp] ERROR: extraction failed"
    echo "$sp" >> "$OUTDIR/.failed_species"
    rm -rf "$t"
    return 0
  fi

  local nprot
  nprot=$(grep -c '^>' "$OUTDIR/${sp}_full_proteins.fa" || true)
  local dur=$(( $(date +%s) - start ))
  echo "[$sp] done in ${dur}s: $ntargets targets -> $nprot proteins"
  rm -rf "$t"
}
export -f process_species ungz filter_gff
export REPO SIF_GFFREAD SIF_PYTHON FAA OUTDIR TMPROOT

# combine over the .faa species (in that order), so a subset re-run of the
# python step never truncates the combined outputs
combine_outputs() {
  local COMBINED="$OUTDIR/all_species_full_proteins.fa"
  local SUMMARY="$OUTDIR/extraction_summary.tsv"
  local f sp
  : > "$COMBINED"
  for sp in "${SPECIES[@]}"; do
    f="$OUTDIR/${sp}_full_proteins.fa"
    [[ -s "$f" ]] && cat "$f" >> "$COMBINED"
  done
  : > "$SUMMARY"
  local header_written=0
  for sp in "${SPECIES[@]}"; do
    f="$OUTDIR/${sp}_extraction_summary.tsv"
    [[ -s "$f" ]] || continue
    if [[ $header_written -eq 0 ]]; then
      head -n 1 "$f" > "$SUMMARY"
      header_written=1
    fi
    tail -n +2 "$f" >> "$SUMMARY"
  done
  local nprots nrows
  nprots=$(grep -c '^>' "$COMBINED" || true)
  nrows=$(( $(wc -l < "$SUMMARY") - 1 )); (( nrows < 0 )) && nrows=0
  echo
  echo "=== full-protein extraction complete ==="
  echo "combined: $COMBINED ($nprots proteins)"
  echo "summary:  $SUMMARY ($nrows rows)"
  if [[ -s "$OUTDIR/.failed_species" ]]; then
    echo "FAILED species:"
    cat "$OUTDIR/.failed_species"
  else
    echo "all species processed OK"
  fi
}
export -f combine_outputs
export SPECIES

printf '%s\n' "${SPECIES[@]}" | xargs -d '\n' -P "$PAR" -I {} \
  bash -c 'process_species "$1"' _ {}

combine_outputs
