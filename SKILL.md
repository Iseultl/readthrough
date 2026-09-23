# SKILL.md — Protist SecORFsearch batch runs

Operating playbook for executing `main_protists.nf` on all 46 protist species on the CRG
Genoa cluster, verifying each run, and analysing each result. Read this file at the start
of every session, and keep the **Status** section up to date as work progresses.

## 1. Mission

Run the SecORFsearch Nextflow pipeline (`main_protists.nf`) sequentially for all 46 protist
species (list + input paths in `protist_filepaths.csv`), **one species at a time**. For each
species: submit → monitor → verify success → run analysis → record metrics in
`run_tracker.tsv` → only then move to the next species.

## 2. Cluster access & division of labour

- **I have direct cluster access**: `ssh ileahy@login1.hpc.crg.es` (key-based auth, no
  password; alias `ssh login` also works — `~/.ssh/config` fixed 2026-09-21).
- **I** run everything on the cluster myself: file sync checks, sbatch submission,
  monitoring (poll `squeue` + logs), verification, analysis, and I update
  `run_tracker.tsv` locally after each species.
- **Code sync** (local → cluster), two acceptable paths:
  1. **scp (default, what I actually do)**: I `scp` each changed/new file straight into
     `/users/rg/ileahy/git/gitlab/readthrough/...` and verify it landed
     (`git status` / `grep` on the cluster). No git commit needed.
  2. **git**: user pushes to gitlab → I `git pull` in the cluster clone and verify.
     Use when the user prefers version control for a change.
- **User's roles**:
  1. Optionally push code changes (see above); otherwise just review.
  2. Review my per-species verification report before I submit the next species.
- Per-run generated files (`runs/params_<sp>.yaml`, `runs/submit_<sp>.sh`) are written
  **directly on the cluster** by me (run artifacts, not code).
- All remote commands: `ssh login '<cmd>'`. Long monitors: poll, don't block
  (use `timeout` / periodic `squeue` checks between messages).

## 3. Key paths

| What | Where |
|---|---|
| Local repo (I edit here) | `/Users/iseult/gitlab/SECIS_independent/CascadeProjects/windsurf-project` |
| Cluster repo clone (code runs from here) | `/users/rg/ileahy/git/gitlab/readthrough` |
| Species input files (source of truth) | `protist_filepaths.csv` in this repo |
| Reference data root (cluster) | `/no_backup/rg/references/species/<Species_dir>.<taxid>/<assembly>/` |
| SLURM logs (cluster) | `/no_backup/rg/ileahy/logs/nf_<sp>_<jobid>.out` / `.err` |
| Nextflow work dirs (cluster) | `/nfs/scratch01/rg/ileahy/nf_work/<sp>` |
| Per-run artifacts (cluster) | `/users/rg/ileahy/git/gitlab/readthrough/runs/` (`params_<sp>.yaml`, `submit_<sp>.sh`) |
| Per-species results (cluster) | `/no_backup/rg/ileahy/<sp>/ORFsearch/` |
| Stale-output backups (cluster) | `/no_backup/rg/ileahy/<sp>/ORFsearch_old_<ts>` |
| Python singularity (cluster) | `~/singularities/python.sif` |
| Run tracker (I maintain) | `run_tracker.tsv` in this repo |
| Analysis script | `analysis/analyse_protist.py` in this repo |

### Per-species result layout (fixed)

```
/no_backup/rg/ileahy/<sp>/ORFsearch/
├── README.txt                                  (genome length, gene count)
├── <sp>_ORFsearch_SECIS.result                 ← final pipeline result
├── secmarker/
├── sequence_logos/<sp>_sequence_logos/
└── analysis/                                   (written by analyse_protist.py)
    ├── <sp>_score_diff.png
    ├── <sp>_summary.csv
    ├── <sp>_candidates.csv
    └── <sp>_result_with_diff.csv
```

This layout requires the two `publishDir` fixes already applied in
`modules/combine_orfsecis.nf` and `modules/extract_sequence_logos.nf` — do not revert them.

## 4. Fixed decisions (agreed with user — do not re-ask)

- All 46 species run fresh, in `protist_filepaths.csv` order (alphabetical).
- `species_name` = reference dir name **without** the `.taxid` suffix
  (e.g. `Babesia_duncani.323732` → `Babesia_duncani`).
- `output_dir = /no_backup/rg/ileahy/<sp>/ORFsearch` (see layout above).
- One species at a time (no parallel runs).
- **No MFE/ML prediction merge** in the analysis — skip that step entirely.
- Selenoprotein prediction criteria: `score_diff >= -1.8 AND re_score_score >= -1.5`,
  where `score_diff = re_score_score - og_score_score`.
- Resources: `--max_cpus 4 --max_memory 16GB`, profile `cluster` (SLURM genoa64, qos=pipelines, Singularity).

## 5. Per-species workflow (the core loop)

For the current species `<sp>`, with input dir `REF` and files from `protist_filepaths.csv`:

### 5.1 Prepare & submit
**First: back up any existing outputs** (several species have stale old-layout results):
```bash
[ -d /no_backup/rg/ileahy/<sp>/ORFsearch ] && \
  mv /no_backup/rg/ileahy/<sp>/ORFsearch /no_backup/rg/ileahy/<sp>/ORFsearch_old_$(date +%Y%m%d_%H%M%S)
```
(Timestamped suffix so a same-day re-run never clobbers an earlier backup.)

**Pre-flight: verify every `.gz` input decompresses to TEXT with a single gunzip**
(Chlamydomonas 2026-09-22 had a double-gzipped GFF — one `gunzip` left binary that
would have broken AGAT). Check `zcat <f> | head -c 16 | od -An -tx1`: if the output
starts `1f 8b` the file is double-compressed — recompress in place
(`zcat <f> | zcat | gzip -c > <f>.fixed && mv <f>.fixed <f>`, keep the original as a
timestamped backup) and verify the transcript feature count.

Then I write `runs/params_<sp>.yaml` and `runs/submit_<sp>.sh` directly into the cluster
repo (`/users/rg/ileahy/git/gitlab/readthrough/runs/`).
Params template:

```yaml
lyric_gtf: "<REF>/<lyric_file>"
genome_fasta: "<REF>/<genome_file>"
genome_gtf: "<REF>/<gff_ref_file>"
species_name: "<sp>"
output_dir: "/no_backup/rg/ileahy/<sp>/ORFsearch"
geneid_param: "<REF>/<param_file>"
max_cpus: 4
max_memory: 16GB
```

Note: for Entamoeba, Hamiltosporidium, Paramecium_tetraurelia, Vairimorpha the CSV `Path`
ends in `/` (files sit directly in that dir) — full path is still `Path/<file>`.

Submit script template (`runs/submit_<sp>.sh`):

```bash
#!/usr/bin/env bash
#SBATCH --no-requeue
#SBATCH --mem 4G
#SBATCH -p genoa64
#SBATCH --qos=pipelines
#SBATCH --mail-type=ALL
#SBATCH --mail-user=iseult.leahy@crg.eu
#SBATCH --output=/no_backup/rg/ileahy/logs/nf_<sp>_%A.out
#SBATCH --error=/no_backup/rg/ileahy/logs/nf_<sp>_%A.err
set -e
module load Java
export NXF_JVM_ARGS="-Xms2g -Xmx5g"
cd /users/rg/ileahy/git/gitlab/readthrough
nextflow run main_protists.nf -params-file runs/params_<sp>.yaml -profile cluster \
    --max_cpus 4 --max_memory 16GB -w /nfs/scratch01/rg/ileahy/nf_work/<sp>
```

(The old `submit.sh` had a broken `wait $pid` with `set -u` — the template above removes it.
Always judge success from the Nextflow log, not the sbatch exit code.)

### 5.2 Monitor
**Poll every ~15 minutes** until the job leaves the queue (user preference 2026-09-22:
advance the next steps WITHOUT waiting for user approval). Practical form: a background
loop that checks `squeue` every 900 s and exits when the job is gone (or re-check
between messages at ~15-min intervals). On confirmed completion (signals below), go
straight through §5.3 → §5.4 → §5.5 and submit the next species.

I poll directly (don't block the session; check between messages):

```bash
# progress: which of my jobs are still alive + their process names
ssh login 'squeue -u ileahy --format="%.10i %.2t %.15j %.20T"'
# completion: definitive exit state of the main sbatch job
ssh login 'sacct -j <jobid> --format=JobID,State,ExitCode | head'
# tail the pipeline log
ssh login 'tail -n 40 /no_backup/rg/ileahy/logs/nf_<sp>_<jobid>.out'
```

**WARNING — do NOT use the log banner as a completion signal.**
`workflowCompletionMessage()` is called inline in the workflow body
(`main_protists.nf:213`), so `Pipeline completed successfully!` is printed at
STARTUP, not on completion. Real completion signals:
1. Main `submit_<sp>` job gone from `squeue`, AND
2. `sacct -j <jobid> --format=State,ExitCode` shows `COMPLETED` / `0:0`, AND
3. Final result file `<sp>_ORFsearch_SECIS.result` exists in the output dir.
Then run the post-run process audit (below) and move to verify.

**Nextflow 25.x layout notes (learned 2026-09-21):**
- Run metadata lives in the LAUNCH dir, not the work dir:
  `cd /users/rg/ileahy/git/gitlab/readthrough && nextflow log` lists runs
  (run name + UUID + start time + status).
- `nextflow log <run-uuid>` on an ACTIVE run fails with a lock error (session
  holds the lock) — while running, monitor via the `.out` log's progress table
  (`162 of 162 ✔` per process) instead.
- `nextflow log <run-uuid>` works AFTER completion: use it to audit that every
  process is `COMPLETED` and count `FAILED`/`ERROR` rows:
  `cd /users/rg/ileahy/git/gitlab/readthrough && nextflow log <run-uuid> | tail -40`
- Requires `module load Java` on the login node first.

### 5.3 Verify (ALL must pass before analysis)
I run directly:

```bash
ssh login 'cat /no_backup/rg/ileahy/<sp>/ORFsearch/README.txt'
ssh login 'ls -la /no_backup/rg/ileahy/<sp>/ORFsearch/'
ssh login 'wc -l /no_backup/rg/ileahy/<sp>/ORFsearch/<sp>_ORFsearch_SECIS.result'
# Transcript counts — take all three BEFORE the work-dir cleanup in 5.5:
# (1) LyRic transcript-level features. Inspect feature types first (format varies
#     per species), then sum the transcript-level types:
ssh login "zcat <lyric_gtf> | awk -F'\t' '!/^#/{print \$3}' | sort | uniq -c"
#     GFF3 merged files (e.g. Babesia_duncani): mRNA + RNA + tRNA
#     plain GTF files: transcript (+ tRNA if present)
#     DO NOT use `grep -c 'transcript'` — on GFF3 it matches the transcript_id
#     attribute on exon/CDS lines and massively over-counts.
# (2) gffread-extracted transcripts (work dir):
ssh login "find /nfs/scratch01/rg/ileahy/nf_work/<sp> -path '*/gffread_out/transcripts_clean_*.fa' -exec cat {} + | grep -c '^>'"
# (3) unique transcripts in the result (row count can exceed this for
#     SPLIT_IF_TOO_LARGE .pN parts):
ssh login "cut -d, -f1 /no_backup/rg/ileahy/<sp>/ORFsearch/<sp>_ORFsearch_SECIS.result | tail -n +2 | sort -u | wc -l"
```

Checklist:
- [ ] `sacct` shows main job `COMPLETED` exit `0:0`; `nextflow log` shows no
      `FAILED`/`ERROR` processes (see 5.2 — the log banner is unreliable)
- [ ] `README.txt` present; genome length + gene count look sane for the species
- [ ] `<sp>_ORFsearch_SECIS.result` exists, non-empty
- [ ] **Transcript counts** (record all three in the tracker):
      `gffread_transcripts` ≥ 99% of `lyric_transcripts`, and
      `result_transcripts` ≥ 99% of `gffread_transcripts` (near-1:1 expected).
      Known quirk (accepted by user, 2026-09-22): gffread occasionally emits broken
      ~70 bp fragments that `recode_any_TGA.py` then silently drops —
      Babesia_duncani: lyric 12,133 → gffread 12,133 → result 12,076 (99.5%).
      A shortfall beyond that → suspect scaffold-name mismatch in the GTF↔FASTA
      join (`main_protists.nf:164` `combine(by: 0)` silently drops non-matching
      scaffolds); investigate before accepting the run.
      Also: transcripts longer than the recode `--limit` are dropped by design
      (`modules/recode_tga.nf`; protists 100,000 bp ≈ none, mammals 8,000).
      Chaetoceros (old hardcoded 8,000) lost 1,320 (1.5%) to this — fixed in 43e4aae.
- [ ] Spot-check `head -3` of the result file: columns `og_score_*`, `re_score_*`,
      `TGA_site_score`, `all_secis_*`, `filtered_secis_*` present.

If any check fails → go to §7 (failure handling). Do NOT proceed to the next species.

### 5.4 Analyse
I run directly (pass the §5.3 counts so the script prints the transcript
reconciliation):
```bash
ssh login 'cd /users/rg/ileahy/git/gitlab/readthrough && \
  singularity run ~/singularities/python.sif python analysis/analyse_protist.py \
    --result /no_backup/rg/ileahy/<sp>/ORFsearch/<sp>_ORFsearch_SECIS.result \
    --species <sp> \
    --outdir /no_backup/rg/ileahy/<sp>/ORFsearch/analysis \
    --gffread-count <N> \
    --lyric-count <M>'
```
I check the printed summary:
- [ ] All 4 output files listed as written
- [ ] `all_secis` / `filtered_secis` counts ≥ 0 and consistent with result file
- [ ] Transcript reconciliation printed: `result ≥ 99% of gffread` and
      `gffread ≥ 99% of lyric` (see §5.3 for the known-quirk tolerance)
- [ ] Summary tables printed for overall / all_secis / filtered_secis
- [ ] Candidate list printed (may legitimately be 0)

### 5.5 Record & advance
- I fill the row for `<sp>` in `run_tracker.tsv`
  (status, job_id, submitted/finished, result_rows, lyric_transcripts,
  gffread_transcripts, all_secis, filtered_secis, predicted_* counts,
  candidates, notes).
- **Cleanup — only when ALL §5.3 checks passed**: remove the Nextflow work dir to
  save space. Final results live under `/no_backup/...`, and run history/metadata
  lives in the launch dir, so `nextflow log` still works afterwards:
  ```bash
  ssh login 'rm -rf /nfs/scratch01/rg/ileahy/nf_work/<sp>'
  ```
  Do the §5.3 gffread count and the `nextflow log <uuid>` audit BEFORE this step.
  On a FAILED run, keep the work dir for §7 debugging.
- I send the user a one-page verification report for `<sp>`, then **automatically
  submit the next species** in `protist_filepaths.csv` order — no approval needed
  (user preference 2026-09-22). Exception: if any §5.3 check fails, stop at §7 and
  report to the user before advancing.

## 6. Batch completion

After all 46: build a consolidated table from `run_tracker.tsv` (species, transcripts,
all/filtered SECIS, predicted counts, candidate count) and flag species with 0 candidates,
0 filtered SECIS, or unusual transcript counts for follow-up.

## 7. Failure handling

1. I pull the diagnostics directly:
   `ssh login 'tail -n 100 /no_backup/rg/ileahy/logs/nf_<sp>_<jobid>.out'` (and `.err`), and
   (after the run has exited; see 5.2 layout notes)
   `ssh login 'module load Java; cd /users/rg/ileahy/git/gitlab/readthrough && nextflow log <run-uuid> | tail -n 40'`
   (lists process states; find the run UUID from `nextflow log` in the same dir).
2. Common suspects:
   - `geneid_param` path wrong / file missing
   - Singularity image pull failure (network / cache `/no_backup/rg/ileahy/Mouse_Analysis/singularity_cache`)
    - Scaffold mismatch: LyRic GTF contig names ≠ genome FASTA headers →
      GFFREAD_CHR / combine steps lose data (see 5.3)
    - SLURM: out of memory / time limit (profile: 8h default per task) — check `.err`
    - Queue/QoS rejection (`qos=pipelines`)
3. Fix locally (params, or code if a bug). Code fixes → sync to cluster (see §2: scp or
   git); params/script fixes → I rewrite them on the cluster directly. Re-submit the same
   species before moving on. Record the issue + fix in the tracker `notes` column.

## 8. Known quirks (context, not action items)

- `main_protists.nf:12` default `geneid_param` is a macOS path — always supply via params file.
- There is NO global `publishDir` any more (removed 2026-09-22, bb7da54): intermediates
  live only in the work dir and are removed after a successful run (§5.5). The old shared
  `/no_backup/rg/ileahy/SecORFsearch/results` folder (268 GB) was a legacy of that block
  (user purging it). The `local` profile still has a harmless `publishDir './results'`
  (only active with `-profile local`).
- Pipeline scripts come from this repo's `bin/`, **bind-mounted into the containers**
  (the readthrough image, Aug 2026, contains none of them) — repo edits to `bin/` take
  effect on the next task run without an image rebuild (verified 2026-09-22 via
  `.command.run` `-B` lines).
- `params.yaml` in the repo root is a stale single-species example (Leishmania) —
  batch runs use `runs/params_<sp>.yaml` only.
- SPLIT_IF_TOO_LARGE (2026-09-23, after e463458): now a single
  `seqkit split2 --by-size 40000 ${input_file}` — no pre-check, no prints.
  Outputs `<input>.fa.split/<base>.part_NNN.fa` (same `part_NNN` naming as the
  old awk split, so `SELECT_INTERESTING`/`GET_ORIGINAL_PREDICTIONS` id-derivation
  from geneid file names is unaffected; empty input → empty .split folder, exit 0,
  no output — verified in the seqkit 2.10.0 container).
- **Re-runs in the same work dir need `-resume` to reuse completed task results**
  (learned 2026-09-23, Conticribra attempt 2): a plain `nextflow run -w <workdir>`
  re-executes EVERY task even when the previous run completed them. `-resume` picks
  up the LATEST session for that work dir — after a new (killed/finished) session
  exists, the older session's cache is no longer reachable. First runs per species
  use the standard template (no -resume); on a re-submit of a stalled/failed run,
  add `-resume` to the nextflow command in `submit_<sp>.sh`.
- Duplicate CSV row: `Pythium_sp._B7052-1` appears twice in `protist_filepaths.csv`;
  run it once.
- The merged gidRef GFFs can contain **orphan CDS records** (transdecoder ORF scans,
  `Parent=<transcript>.pN` with no corresponding transcript line) — gffread materializes
  each as a spurious CDS-only "transcript". Chlorella_sorokiniana: 19,907 spurious
  (46,367 gffread vs 26,460 real; 19,907 orphan parents, CDS-only, no orphan exons).
  Removed by the `FILTER_ORPHAN_CDS` process (`bin/filter_orphan_cds.py`), wired between
  UNZIP_IF_NEEDED and AGAT_SPLITGFF (2026-09-22). If a run shows gffread > lyric
  transcripts, suspect this.

## 9. Status (update every session)

- Phase: **1 — batch running autonomously** (user-approved 2026-09-22): poll every
  15 min → on completion verify → analyse → record → cleanup → next species, no approval
- Current species: species 7 — Cyanidiococcus_yangmingshanensis (job 28725643,
  submitted 2026-09-23 13:10)
- Deferred: species 5 — Conticribra_weissflogii. Attempt 2 (28713877) scancelled
  2026-09-23 12:35 at user decision (16–25 h ETA too slow; user may run it over the
  weekend when the cluster is less busy). Work dir
  `/nfs/scratch01/rg/ileahy/nf_work/Conticribra_weissflogii` KEPT — a re-run with
  `-resume` reuses the attempt-2 partial cache (AGAT stage ~done).
- Completed: 5 / 46
  - Babesia_duncani (job 28607869): lyric 12,133 → gffread 12,133 → result 12,076 unique
    (99.5%, gffread quirk accepted); 1 candidate (agat-rna-5082, score_diff 0.81)
  - Chaetoceros_neogracilis (job 28633444): lyric 88,400 → gffread 88,400 → result
    87,080 unique (98.5% — 1,320 transcripts >8 kb dropped by the old hardcoded recode
    limit; NOT re-run per user); 6 candidates, top agat-rna-20420 (score_diff 15.75)
  - Chlamydomonas_reinhardtii (job 28639838): lyric 19,527 → gffread 19,526 → result
    19,526 unique (100%); first run with --limit 100000 (long transcripts included);
    1 candidate (rna-XM_001696020.2, score_diff 0.37)
  - Chlorella_sorokiniana (job 28645869, re-run after orphan-CDS fix 51610a3):
    lyric 26,460 → gffread 26,460 → result 26,460 unique (100% three-way match;
    attempt 1, 28642581, had 19,907 spurious CDS-only transcripts from transdecoder
    ORF records — user-diagnosed, fixed by FILTER_ORPHAN_CDS; attempt-1 output
    backed up ORFsearch_old_20260922_171225); 3 candidates, top agat-rna-13006 (2.12)
  - Cryptosporidium_parvum_Iowa_II (job 28721317): lyric 17,956 → gffread 17,956 →
    result 17,956 unique (100% three-way match; first run with seqkit split2 SPLIT
    734f846); lyric input was a merged gidRef (AGAT LyRic 14,075 + RefSeq 3,881);
    8 scaffolds, 18m39s; 119 SECIS found, 0 survived the filter → 0 candidates;
    old Apr-08 output backed up ORFsearch_old_20260923_123718
- Done this session (2026-09-22):
  - [x] Pipeline fixes: sequence_logos per-scaffold overwrite (bb7da54 + e87c972:
        `collectFile` on dirs → `collect()` of files), global publishDir removed
        (bb7da54, intermediates in work dir only), logo header fix
  - [x] Chlamydomonas double-gzipped reference GFF fixed in place (backup kept)
  - [x] SKILL.md: 15-min polling + autonomous species advance (§5.2/§5.5), pre-flight
        gz-integrity check (§5.1), recode --limit note (§5.3), publishDir quirk corrected
        + bind-mounted bin/ note (§8)
  - [x] Chaetoceros closed out: analysis (6 candidates), work dir removed (42 GB)
  - [x] Orphan-CDS fix: bin/filter_orphan_cds.py + FILTER_ORPHAN_CDS wired between
        UNZIP and AGAT_SPLITGFF (51610a3); Chlorella re-run 28645869 closed out:
        100% three-way, 3 candidates
- Done this session (2026-09-23):
  - [x] Diagnosed Conticribra attempt 1 (28647231): user scancelled 10:12 after a
        16h16m stall. Root cause: `maxForks 1` in SPLIT_IF_TOO_LARGE (leftover from
        c0ffc43 "run as github actions") serialized the per-scaffold stage — 7,719
        scaffolds × ~65s at 1 job at a time ≈ 140h; only 639/7719 done. Other stages
        ran 18–38 parallel (verified via sacct child-job overlap). Previous 4 species
        had ≤162 scaffolds so the bug was masked. Fix e463458 removes the
        process-level maxForks (CI unaffected: github profile sets maxForks=1
        globally). Attempt 2 = 28713877 (submitted 11:34); attempt-1 partial output
        backed up ORFsearch_old_20260923_113417; pre-flight re-passed
        (all .gz OK, lyric 14,490 mRNA). Caveat: template used plain `nextflow run`
        → NO cache reuse from attempt 1 (would have needed -resume), all stages
        re-run but SPLIT now parallel; ETA ~16-25h. See §8 for the -resume rule.
   - [x] SPLIT_IF_TOO_LARGE simplified per user: single `seqkit split2 --by-size
         40000 ${input_file}` (no pre-check, no prints) — commit 734f846. Produces
         `<input>.fa.split/<base>.part_NNN.fa` — same `part_NNN` naming as the old
         awk split (verified in the seqkit 2.10.0 container: 85k-seq multi-part and
         empty-file cases), so downstream id derivation from geneid file names is
         unaffected.
   - [x] Conticribra attempt 2 (28713877) scancelled 12:35 at user decision;
         Conticribra deferred to the weekend (work dir kept for a -resume re-run).
   - [x] Species 6 Cryptosporidium_parvum_Iowa_II submitted: job 28721317 (12:43).
         Pre-flight: all 4 inputs OK (single-gz); merged gidRef has 17,956
         transcript-level features (AGAT LyRic 14,075 + RefSeq 3,881) and only 8
         scaffolds → run expected ~30-60 min; local re-run of the GFF chain
         (orphan-CDS 0 removed, AGAT GFF2GTF, gffread per chr) = 17,956 transcripts,
         0 warnings. Old Apr-08 output backed up ORFsearch_old_20260923_123718.
    - [x] Cryptosporidium closed out: 28721317 COMPLETED 0:0 in 18m39s (12:43→13:03);
          100% three-way match lyric=gffread=result=17,956; 119 SECIS found, 0
          survived the filter → 0 candidates; work dir removed. Batch now 5/46.
    - [x] Species 7 Cyanidiococcus_yangmingshanensis submitted: job 28725643 (13:10).
          Pre-flight: reference gidRef was DOUBLE-GZ (3rd occurrence) — fixed in
          place, orig kept *_doublegz_backup_20260923; GenBank reference annotation
          5,189 mRNA (not a LyRic transcriptome, like Conticribra); 20 scaffolds →
          fast run expected; old Apr-08 output backed up ORFsearch_old_20260923_131019.
- Next:
  - [ ] Monitor 28725643 (15-min polls) → verify → analyse → record
        → next species (Cyanidioschyzon_merolae_strain_10D)
- Notes: this session ran **directly on the cluster login node** (genoa64-05, user
  ileahy) — no `ssh login` prefix needed; from the local machine use the `ssh login`
  forms as written. Re-read this file at the start of each session and keep this
  section current.
- Session 2026-09-23 (genoa64-04): code fix e463458 was committed directly in the
  CLUSTER clone — the local Mac repo
  (`/Users/iseult/gitlab/SECIS_independent/CascadeProjects/windsurf-project`) needs a
  `git pull` to pick it up.
