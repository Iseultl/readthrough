#!/usr/bin/env bash
#SBATCH --no-requeue
#SBATCH --mem 16G
#SBATCH -p genoa64
#SBATCH --qos=pipelines
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ileahy@crg.es
#SBATCH --output=/no_backup/rg/ileahy/logs/nf_Cyanidiococcus_yangmingshanensis_%A.out
#SBATCH --error=/no_backup/rg/ileahy/logs/nf_Cyanidiococcus_yangmingshanensis_%A.err
#
# Cyanidiococcus_yangmingshanensis: SecORFsearch pipeline run + analysis + work-dir cleanup.
# Success = nextflow exit 0 AND non-empty result file.
set -euo pipefail

SP="Cyanidiococcus_yangmingshanensis"
REPO="/users/rg/ileahy/git/gitlab/readthrough"
OUT="/no_backup/rg/ileahy/${SP}/ORFsearch"
WORK="/nfs/scratch01/rg/ileahy/nf_work/${SP}"
RESULT="${OUT}/${SP}_ORFsearch_SECIS.result"
LYRIC_GTF="/no_backup/rg/references/species/Cyanidiococcus_yangmingshanensis.2690220/GCA_013995675.1/LyRic_complete.gff.gz"

module load Java
export NXF_JVM_ARGS="-Xms2g -Xmx5g"
cd "$REPO"

# 0) Back up any pre-existing outputs (stale old-layout results)
if [ -d "$OUT" ]; then
    mv "$OUT" "${OUT}_old_$(date +%Y%m%d_%H%M%S)"
    echo "Backed up existing $OUT"
fi

# 1) Run the pipeline (fresh; no -resume: input content changed).
#    Retry ONLY on the transient Nextflow session-lock race: concurrent
#    `nextflow run`s started near-simultaneously from this launch dir on NFS
#    can lose the shared .nextflow session lock ("Unable to acquire lock on
#    session" / "No locks available"). 4 jobs in the 2026-10-19 batch died
#    this way <30s after start. Any OTHER failure aborts immediately (work
#    dir kept for debugging) - do not blindly re-run a genuinely broken run.
JLOG="/no_backup/rg/ileahy/logs/nf_${SP}_${SLURM_JOBID}.out"
attempt=0
while :; do
    attempt=$((attempt+1))
    if nextflow run main_protists.nf -params-file runs/params_${SP}.yaml -profile cluster \
        --max_cpus 4 --max_memory 16GB -w "$WORK"; then
        break
    fi
    sleep 2   # let nextflow flush its final stderr/stdout to $JLOG
    if [ "$attempt" -lt 5 ] && grep -qE "Unable to acquire lock on session|No locks available" "$JLOG" 2>/dev/null; then
        echo "WARN: nextflow session-lock race (attempt $attempt) - retrying in $((attempt*120))s"
        sleep $((attempt * 120))
        continue
    fi
    echo "ERROR: nextflow failed (attempt $attempt) - keeping work dir for debugging" >&2
    exit 1
done

# 2) Sanity check: final result exists and is non-empty
if [ ! -s "$RESULT" ]; then
    echo "ERROR: $RESULT missing or empty - keeping work dir for debugging" >&2
    exit 1
fi

# 3) Transcript counts (MUST be taken before the work-dir cleanup)
#    lyric: sum of transcript-level feature types (format varies per species)
LYRIC_COUNT=$(zcat "$LYRIC_GTF" | awk -F'\t' '!/^#/ && ($3=="mRNA" || $3=="transcript" || $3=="RNA" || $3=="rRNA" || $3=="tRNA" || $3=="scRNA" || $3=="circRNA" || $3=="ncRNA" || $3=="lncRNA")' | wc -l)
# count UNIQUE headers, not raw lines: a -resume re-run may re-execute every
# task in fresh task dirs while the previous session's task dirs (identical
# outputs) remain in the work dir - a plain line count then sums both
# generations (Trichomonas_vaginalis 29112744 reported 259952 = 2x the true
# 129976). Duplicate ids cannot occur in a successful run (Seblastian rejects
# them), so dedup by id is lossless.
GFFREAD_COUNT=$(find "$WORK" -path '*/gffread_out/transcripts_clean_*.fa' -exec cat {} + | grep '^>' | awk '{print $1}' | sort -u | wc -l || true)
RESULT_COUNT=$(cut -d, -f1 "$RESULT" | tail -n +2 | sort -u | wc -l)
echo "Transcript counts: lyric=$LYRIC_COUNT gffread=$GFFREAD_COUNT result=$RESULT_COUNT"

# 4) Analyse (python in singularity)
singularity run ~/singularities/python.sif python analysis/analyse_protist.py \
    --result "$RESULT" \
    --species "$SP" \
    --outdir "$OUT/analysis" \
    --gffread-count "$GFFREAD_COUNT" \
    --lyric-count "$LYRIC_COUNT"

# 5) Cleanup: remove Nextflow work dir to save space
#    (final results live under /no_backup; run metadata lives in the launch dir)
rm -rf "$WORK"
echo "DONE: $SP (lyric=$LYRIC_COUNT gffread=$GFFREAD_COUNT result=$RESULT_COUNT)"
