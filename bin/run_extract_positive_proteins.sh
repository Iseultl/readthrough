#!/usr/bin/env bash
#
# Extract positive, SECIS-containing protein sequences for all (or a subset of)
# protist species.
#
# Per species (driven by runs/params_<sp>.yaml):
#   1) bin/extract_positive_proteins.py --select-only
#        selects positives (score_diff >= -1.8 AND re_score_score >= -1.5)
#        with at least one identified SECIS element (all_secis_start not NA)
#        from <output_dir>/<sp>_ORFsearch_SECIS.result
#   2) gffread -F -w transcripts.fa -g <genome_fasta> <lyric_gtf>
#        (same flags as the pipeline's GFFREAD_CHR, on the ORIGINAL inputs)
#   3) bin/extract_positive_proteins.py
#        translates the recoded CDS (re_score_start..re_score_end) and writes
#        the protein with selenocysteine (U) at the in-frame TGA
#
# Outputs (in --outdir, default /no_backup/rg/ileahy/positive_proteins):
#   <sp>_positive_proteins.fa        per species
#   <sp>_extraction_summary.tsv      per protein
#   <sp>_extraction_log.tsv          skipped transcripts + reasons
#   <sp>_targets.tsv                 selected transcript list
#   all_species_positive_proteins.fa combined, species in input order
#   extraction_summary.tsv           combined summary
#   .failed_species                  species that errored (empty = all OK)
#
# Usage:
#   run_extract_positive_proteins.sh [--outdir DIR] [--parallel N] [--tmp-root DIR] [species ...]
#   no species given = all runs/params_*.yaml

set -euo pipefail

REPO="/users/rg/ileahy/git/gitlab/readthrough"
SIF_GFFREAD="${SIF_GFFREAD:-$HOME/singularities/gffread.sif}"
SIF_PYTHON="${SIF_PYTHON:-$HOME/singularities/python.sif}"
OUTDIR="/no_backup/rg/ileahy/positive_proteins"
PAR=4
TMPROOT="/nfs/scratch01/rg/ileahy/positive_proteins_tmp"
COMBINE_ONLY=0
PER_SCAFFOLD=0
AGAT_CHAIN=0
SIF_AGAT="${SIF_AGAT:-$HOME/singularities/agat.sif}"
SPECIES=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir)       OUTDIR="$2"; shift 2 ;;
    --parallel)     PAR="$2"; shift 2 ;;
    --tmp-root)     TMPROOT="$2"; shift 2 ;;
    --combine-only) COMBINE_ONLY=1; shift ;;
    # Split GFF/FASTA by scaffold and gffread in batches of 50 scaffolds:
    # whole-genome gffread OOMs (login-node user memory cgroup) on the two
    # largest genomes (Haematococcus 820 scaffolds, Plasmodium_vivax 2748).
    --per-scaffold) PER_SCAFFOLD=1; shift ;;
    # Replicate the pipeline's GFF chain (orphan-CDS filter -> AGAT gff2gtf
    # -> clean_gtf -> gffread). Needed when the reference GFF is EMBL-style
    # (CDS features parented to gene, siblings of the mRNA) and AGAT
    # reconstructs the exons from the CDS, so direct gffread yields a
    # DIFFERENT transcript sequence (Heterosigma: 40,884 such CDS lines).
    --agat-chain)   AGAT_CHAIN=1; shift ;;
    *)              SPECIES+=("$1"); shift ;;
  esac
done

# Always combine over EVERY species with a params file (in that order), so a
# subset re-run never truncates the combined outputs.
combine_outputs() {
  local COMBINED="$OUTDIR/all_species_positive_proteins.fa"
  local SUMMARY="$OUTDIR/extraction_summary.tsv"
  local f sp
  local ALL=()
  for p in "$REPO"/runs/params_*.yaml; do
    ALL+=("$(basename "$p" .yaml | sed 's/^params_//')")
  done
  : > "$COMBINED"
  for sp in "${ALL[@]}"; do
    f="$OUTDIR/${sp}_positive_proteins.fa"
    [[ -s "$f" ]] && cat "$f" >> "$COMBINED"
  done
  : > "$SUMMARY"
  local header_written=0
  for sp in "${ALL[@]}"; do
    f="$OUTDIR/${sp}_extraction_summary.tsv"
    [[ -s "$f" ]] || continue
    if [[ $header_written -eq 0 ]]; then
      head -n 1 "$f" > "$SUMMARY"
      header_written=1
    fi
    tail -n +2 "$f" >> "$SUMMARY"
  done
  echo
  echo "=== extraction complete ==="
  echo "combined: $COMBINED ($(grep -c '^>' "$COMBINED" || true) proteins)"
  echo "summary:  $SUMMARY ($(( $(wc -l < "$SUMMARY") - 1 )) rows)"
  if [[ -s "$OUTDIR/.failed_species" ]]; then
    echo "FAILED species:"
    cat "$OUTDIR/.failed_species"
  else
    echo "all species processed OK"
  fi
}
export -f combine_outputs

if [[ $COMBINE_ONLY -eq 1 ]]; then
  combine_outputs
  exit 0
fi

# Default: every species that has a params file
if [[ ${#SPECIES[@]} -eq 0 ]]; then
  for p in "$REPO"/runs/params_*.yaml; do
    SPECIES+=("$(basename "$p" .yaml | sed 's/^params_//')")
  done
fi

mkdir -p "$OUTDIR" "$TMPROOT"
: > "$OUTDIR/.failed_species"

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

# reduce_gff_attrs <in.gff> <out.gff> — the attribute reduction done by
# modules/agat_gff2gtf.nf before agat_convert_sp_gff2gtf.pl (keep only
# gene_id/transcript_id or ID/Parent).
reduce_gff_attrs() {
  awk -F'\t' 'BEGIN { OFS="\t" }
  /^#/ { print; next }
  {
    original_attrs = $9
    gene_id = ""; transcript_id = ""; id = ""; parent = ""
    n = split($9, attrs, ";")
    for (i = 1; i <= n; i++) {
      attr = attrs[i]
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", attr)
      if (attr ~ /^gene_id[[:space:]]+/) gene_id = attr
      else if (attr ~ /^transcript_id[[:space:]]+/) transcript_id = attr
      else if (attr ~ /^ID=/) id = attr
      else if (attr ~ /^Parent=/) parent = attr
    }
    out = ""
    if (gene_id != "" || transcript_id != "") {
      if (gene_id != "") out = gene_id
      if (transcript_id != "") out = (out == "" ? transcript_id : out "; " transcript_id)
      out = out ";"
    }
    else if (id != "" || parent != "") {
      if (id != "") out = id
      if (parent != "") out = (out == "" ? parent : out ";" parent)
    }
    else { out = original_attrs }
    $9 = out
    print
  }' "$1" > "$2"
}

# split_scaffolds <gff> <fasta> <dir> — write per-scaffold files
# <dir>/sc_<scaffold>.gff and <dir>/pf_<scaffold>.fa
# Streaming split: the reference GFFs/FASTAs are sorted by scaffold, so at
# most one file is open at a time. (Open/close per LINE was pathologically
# slow on NFS — ~1 scaffold/min on a big genome.)
split_scaffolds() {
  singularity exec -B "$3" -B /no_backup "$SIF_PYTHON" python - "$1" "$2" "$3" <<'PYEOF'
import sys
gff, fa, t = sys.argv[1], sys.argv[2], sys.argv[3]

def split_stream(path, prefix, keyfn):
    handles, cur = {}, None
    with open(path) as f:
        for line in f:
            name = keyfn(line)
            if name is None or name.startswith("#"):
                continue
            if name != cur:
                if cur is not None and cur in handles:
                    handles.pop(cur).close()
                if name not in handles:
                    handles[name] = open(f"{t}/{prefix}{name}", "a")
                cur = name
            handles[cur].write(line)
    for h in handles.values():
        h.close()

split_stream(gff, "sc_", lambda l: l.split("\t", 1)[0] + ".gff" if not l.startswith("#") else None)

# FASTA: key = first header word of the current record
handles, cur = {}, None
with open(fa) as f:
    for line in f:
        if line.startswith(">"):
            cur = line[1:].split()[0]
        if cur is None:
            continue
        if cur not in handles:
            handles[cur] = open(f"{t}/pf_{cur}.fa", "a")
        handles[cur].write(line)
for h in handles.values():
    h.close()
PYEOF
}

# run_gffread_batch <dir> <sc_*.gff ...> — one gffread call per batch of
# scaffolds (re-paired by filename, like the pipeline's GFFREAD_CHR),
# appending to <dir>/transcripts.fa
run_gffread_batch() {
  local t="$1"; shift
  local args=() g s
  for g in "$@"; do
    s=$(basename "$g" .gff | sed 's/^sc_//')
    if [[ -f "$t/pf_$s.fa" ]]; then
      args+=(-g "$t/pf_$s.fa")
    fi
  done
  if ! singularity exec -B "$t" -B /no_backup "$SIF_GFFREAD" \
        gffread -F -w "$t/part.fa" "${args[@]}" "$@" \
        >> "$t/gffread.out" 2>&1; then
    return 1
  fi
  cat "$t/part.fa" >> "$t/transcripts.fa"
  rm -f "$t/part.fa"
}

process_species() {
  local sp="$1"
  local params="$REPO/runs/params_$sp.yaml"
  local start; start=$(date +%s)
  if [[ ! -f "$params" ]]; then
    echo "[$sp] SKIP: no runs/params_$sp.yaml"
    return 0
  fi
  local species lyric genome outdir result
  species=$(awk -F'"' '/^species_name:/{print $2}' "$params")
  outdir=$(awk -F'"' '/^output_dir:/{print $2}' "$params")
  lyric=$(awk -F'"' '/^lyric_gtf:/{print $2}' "$params")
  genome=$(awk -F'"' '/^genome_fasta:/{print $2}' "$params")
  result="$outdir/${species}_ORFsearch_SECIS.result"

  if [[ ! -s "$result" ]] || [[ $(wc -l < "$result") -le 1 ]]; then
    echo "[$sp] SKIP: result file missing or header-only ($result)"
    : > "$OUTDIR/${species}_positive_proteins.fa"
    return 0
  fi

  local t="$TMPROOT/$sp"
  mkdir -p "$t"

  # 1) fast target selection (no genome access needed)
  if ! singularity exec -B "$t" -B /no_backup -B "$REPO" -B "$OUTDIR" "$SIF_PYTHON" \
        python "$REPO/bin/extract_positive_proteins.py" \
        --result "$result" --species "$species" --outdir "$OUTDIR" --select-only; then
    echo "[$sp] ERROR: selection failed"
    echo "$sp" >> "$OUTDIR/.failed_species"
    rm -rf "$t"
    return 0
  fi
  local ntargets
  ntargets=$(( $(wc -l < "$OUTDIR/${species}_targets.tsv") - 1 ))
  if [[ $ntargets -eq 0 ]]; then
    echo "[$sp] no targets — empty outputs written"
    rm -rf "$t"
    return 0
  fi

  # 2) gunzip the original inputs (never modify the references)
  local gname lname
  gname=$(basename "$genome" .gz)
  lname=$(basename "$lyric" .gz)
  if ! ungz "$genome" "$t/$gname" || ! ungz "$lyric" "$t/$lname"; then
    echo "[$sp] ERROR: failed to decompress inputs ($genome / $lyric)"
    echo "$sp" >> "$OUTDIR/.failed_species"
    rm -rf "$t"
    return 0
  fi

  # 3) transcript extraction. Three modes:
  #    default      — gffread on the whole genome (same command as the
  #                   pipeline's GFFREAD_CHR); OK up to ~250 Mb genomes
  #    --per-scaffold — split GFF/FASTA by scaffold, gffread in batches of
  #                   50 scaffolds per call (whole-genome gffread OOMs in
  #                   the login-node user memory cgroup for Haematococcus
  #                   and Plasmodium_vivax)
  #    --agat-chain — replicate the pipeline's GFF chain (orphan-CDS filter
  #                   -> AGAT gff2gtf -> clean_gtf -> gffread); needed for
  #                   EMBL-style GFFs where AGAT reconstructs exons from
  #                   the CDS (Heterosigma)
  local step3_ok=1
  if [[ $AGAT_CHAIN -eq 1 ]]; then
    if ! singularity exec -B "$t" -B /no_backup -B "$REPO" "$SIF_PYTHON" \
          python "$REPO/bin/filter_orphan_cds.py" \
          --gff "$t/$lname" --output "$t/orphan_free.gff" \
      || ! reduce_gff_attrs "$t/orphan_free.gff" "$t/reduced.gff" \
      || ! singularity exec -B "$t" "$SIF_AGAT" bash -c \
          ". /usr/local/env-activate.sh; agat_convert_sp_gff2gtf.pl --gff $t/reduced.gff --gtf_version 3 --output $t/converted.gtf" \
          > "$t/agat.out" 2>&1 \
      || ! bash "$REPO/bin/clean_gtf.sh" "$t/converted.gtf" > /dev/null; then
      step3_ok=0
    else
      singularity exec -B "$t" -B /no_backup "$SIF_GFFREAD" \
          gffread -F -w "$t/transcripts.fa" -g "$t/$gname" "$t/converted.cleaned.gtf" \
          > "$t/gffread.out" 2>&1 || step3_ok=0
    fi
  elif [[ $PER_SCAFFOLD -eq 1 ]]; then
    if ! split_scaffolds "$t/$lname" "$t/$gname" "$t"; then
      step3_ok=0
    else
      : > "$t/transcripts.fa"
      local batch=() gff s
      for gff in "$t"/sc_*.gff; do
        batch+=("$gff")
        if (( ${#batch[@]} >= 50 )); then
          run_gffread_batch "$t" "${batch[@]}" || step3_ok=0
          batch=()
          (( step3_ok )) || break
        fi
      done
      if (( step3_ok )) && (( ${#batch[@]} > 0 )); then
        run_gffread_batch "$t" "${batch[@]}" || step3_ok=0
      fi
    fi
  else
    singularity exec -B "$t" -B /no_backup "$SIF_GFFREAD" \
        gffread -F -w "$t/transcripts.fa" -g "$t/$gname" "$t/$lname" \
        > "$t/gffread.out" 2>&1 || step3_ok=0
  fi
  if [[ $step3_ok -eq 0 ]]; then
    echo "[$sp] ERROR: transcript extraction failed"
    tail -n 20 "$t/gffread.out" 2>/dev/null
    echo "$sp" >> "$OUTDIR/.failed_species"
    rm -rf "$t"
    return 0
  fi

  # 4) translate selected CDS, U at the TGA
  if ! singularity exec -B "$t" -B /no_backup -B "$REPO" -B "$OUTDIR" "$SIF_PYTHON" \
        python "$REPO/bin/extract_positive_proteins.py" \
        --result "$result" --transcripts "$t/transcripts.fa" \
        --species "$species" --outdir "$OUTDIR"; then
    echo "[$sp] ERROR: translation failed"
    echo "$sp" >> "$OUTDIR/.failed_species"
    rm -rf "$t"
    return 0
  fi

  local nprot
  nprot=$(grep -c '^>' "$OUTDIR/${species}_positive_proteins.fa" || true)
  local dur=$(( $(date +%s) - start ))
  echo "[$sp] done in ${dur}s: $ntargets targets -> $nprot proteins"
  rm -rf "$t"
}
export -f process_species ungz reduce_gff_attrs split_scaffolds run_gffread_batch
export REPO SIF_GFFREAD SIF_PYTHON SIF_AGAT OUTDIR TMPROOT PER_SCAFFOLD AGAT_CHAIN

printf '%s\n' "${SPECIES[@]}" | xargs -d '\n' -P "$PAR" -I {} \
  bash -c 'process_species "$1"' _ {}

# ---- combine (over all species with params files, in that order) ----
combine_outputs
