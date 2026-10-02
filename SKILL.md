# SKILL.md — Protist SecORFsearch batch runs

Operating playbook for executing `main_protists.nf` on all 46 protist species on the CRG
Genoa cluster, verifying each run, and analysing each result. Read this file at the start
of every session, and keep the **Status** section up to date as work progresses.

## 1. Mission

Run the SecORFsearch Nextflow pipeline (`main_protists.nf`) for all 46 protist species
(list + input paths in `protist_filepaths.csv`). Since **2026-09-30** this is done in
**batch format**: one self-contained SLURM job per species (pipeline → verify → analyse →
work-dir cleanup) and a single wrapper (`runs/submit_all_remaining.sh`) that submits all
remaining species at once. I no longer trigger species one at a time; instead I monitor
the batch, and as each job completes I verify its log and record metrics in
`run_tracker.tsv`.

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
- Per-run generated files (`runs/params_<sp>.yaml`, `runs/submit_<sp>.sh`) are generated
  **on the cluster** for all remaining species at once (batch generator, see §5.1);
  they live in the repo `runs/` dir and count as run artifacts, not code.
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
- **Batch submission (2026-09-30, user decision)**: all remaining species are submitted
  at once via `runs/submit_all_remaining.sh` (one sbatch job per species); SLURM
  schedules the main jobs, each job's child tasks schedule themselves. Replacing the old
  "one species at a time" rule.
- **No MFE/ML prediction merge** in the analysis — skip that step entirely.
- Selenoprotein prediction criteria: `score_diff >= -1.8 AND re_score_score >= -1.5`,
  where `score_diff = re_score_score - og_score_score`.
- Resources: `--max_cpus 4 --max_memory 16GB`, profile `cluster` (SLURM genoa64, qos=pipelines, Singularity).

## 5. Per-species workflow (batch format, since 2026-09-30)

One self-contained SLURM job per species does: **backup stale outputs → run pipeline →
sanity-check result → transcript counts → analysis → remove work dir**. A single wrapper
(`runs/submit_all_remaining.sh`) submits all remaining species at once. I no longer
trigger species one at a time; instead I monitor the batch and close out each species as
its job completes (verify log → record tracker row). The old manual one-species-at-a-time
loop is kept in §5.6 for reference.

### 5.1 Artifacts (generated on the cluster, in `runs/`)

A batch generator (reads `run_tracker.tsv` status + `protist_filepaths.csv`) writes for
every remaining species:

- `params_<sp>.yaml` — same fields as the legacy template (lyric_gtf, genome_fasta,
  genome_gtf, species_name, output_dir, geneid_param, max_cpus 4, max_memory 16GB,
  scaffold_batches 24).
- `submit_<sp>.sh` — the self-contained job:

```bash
#!/usr/bin/env bash
#SBATCH --no-requeue
#SBATCH --mem 16G            # was 4G: the in-job analysis OOM-killed a 180k-row
#SBATCH -p genoa64          # result (Eimeria, 09-30); job hosts JVM (-Xmx5g) + analysis
#SBATCH --qos=pipelines
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ileahy@crg.es
#SBATCH --output=/no_backup/rg/ileahy/logs/nf_<sp>_%A.out
#SBATCH --error=/no_backup/rg/ileahy/logs/nf_<sp>_%A.err
set -euo pipefail
# 0) back up existing /no_backup/rg/ileahy/<sp>/ORFsearch -> ORFsearch_old_<ts>
# 1) module load Java; NXF_JVM_ARGS="-Xms2g -Xmx5g"; cd repo
#    nextflow run main_protists.nf -params-file runs/params_<sp>.yaml -profile cluster \
#        --max_cpus 4 --max_memory 16GB -w /nfs/scratch01/rg/ileahy/nf_work/<sp>
#        (Conticribra only: -resume — reuses the attempt-2 partial cache)
# 2) abort if <sp>_ORFsearch_SECIS.result missing/empty (work dir KEPT on failure)
# 3) print "Transcript counts: lyric=... gffread=... result=..." — gffread count from
#    work-dir transcripts_clean_*.fa BEFORE step 5; lyric count = transcript-level
#    feature types (mRNA/transcript/RNA/rRNA/tRNA/scRNA/circRNA/ncRNA/lncRNA)
# 4) singularity run ~/singularities/python.sif python analysis/analyse_protist.py \
#        --result <RESULT> --species <sp> --outdir <OUT>/analysis \
#        --gffread-count N --lyric-count M
# 5) rm -rf /nfs/scratch01/rg/ileahy/nf_work/<sp>   (reached only on success)
#    echo "DONE: <sp> (lyric=... gffread=... result=...)"
```

- `submit_all_remaining.sh` — **the wrapper**: `sbatch`es every remaining species and
  prints each job id. Re-running it re-submits EVERY listed species — scancel the
  previous batch first.

Pre-flight (double-gz) check is no longer done species-by-species up front: a
double-gzipped input now simply fails that species' job at the UNZIP stage (work dir
kept, §7). Fix the file in place (keep a `*_doublegz_backup_<date>`), then re-submit that
one species only (`sbatch runs/submit_<sp>.sh`).

### 5.1a Legacy per-species preparation (reference)
The old flow backed up existing outputs first (timestamped suffix so a same-day re-run
never clobbers an earlier backup), then pre-flight-checked every `.gz` input decompresses
to TEXT with a single gunzip (`zcat <f> | head -c 16 | od -An -tx1`; `1f 8b` =
double-compressed — recompress in place with `zcat <f> | zcat | gzip -c`), then wrote
`params_<sp>.yaml` with:

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
scaffold_batches: 24
```

`scaffold_batches` (default 24, since 2026-09-28): the 6 per-scaffold stages of
`main_protists.nf` (agat gff2gtf, clean gtf, gffread, recode TGA, split, secissearch)
run as 24 hash-bucket jobs each instead of one job per scaffold (small-scaffold
genomes like Cylindrotheca/Conticribra used to spawn 15k-46k short SLURM jobs).
Peak memory is unchanged (sequential loop inside each job). See §8 for the
NF 26.04.6 channel/staging rules the implementation depends on.

Note: for Entamoeba, Hamiltosporidium, Paramecium_tetraurelia, Vairimorpha the CSV `Path`
ends in `/` (files sit directly in that dir) — full path is still `Path/<file>`.

The legacy submit template (pipeline only, no analysis/cleanup) is obsolete — the
batch-format template above (§5.1) supersedes it. Still applies: judge success from the
Nextflow log / result file, not the sbatch exit code banner; the "Pipeline completed
successfully!" banner prints at STARTUP (`workflowCompletionMessage()` inline).

### 5.2 Monitor the batch
Poll directly, don't block the session (check between messages):

```bash
# progress: which of my jobs are still alive + their names
squeue -u ileahy --format="%.10i %.2t %.15j %.20T"
# completion: definitive exit state of a main job
sacct -j <jobid> --format=JobID,State,ExitCode | head
# tail the job log (success = the LAST lines, see below)
tail -n 40 /no_backup/rg/ileahy/logs/nf_<sp>_<jobid>.out
```

For each species job that leaves the queue, completion = ALL of:
1. `sacct` shows `COMPLETED` / `0:0`;
2. the `.out` log's final lines contain the job's own
   `Transcript counts: lyric=... gffread=... result=...` and `DONE: <sp> ...`
   (NOT the startup `Pipeline completed successfully!` banner — see §5.1a warning);
3. the 4 analysis files exist in `/no_backup/rg/ileahy/<sp>/ORFsearch/analysis/`.
The work dir has already been removed by the job itself (step 5 of the submit script);
on any failure it is KEPT for §7 debugging.

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

### 5.3 Verify (post-hoc, per completed job)
The job already did: result-file sanity check, transcript counts (printed to the log),
analysis, and work-dir cleanup. My verification now reads the evidence:

```bash
cat /no_backup/rg/ileahy/<sp>/ORFsearch/README.txt
ls /no_backup/rg/ileahy/<sp>/ORFsearch/analysis/          # 4 files expected
tail -n 20 /no_backup/rg/ileahy/logs/nf_<sp>_<jobid>.out  # counts + DONE line
```

Checklist:
- [ ] `sacct` shows main job `COMPLETED` exit `0:0` and the log has the `DONE:` line
- [ ] **Transcript counts** (read from the log, record in the tracker):
      `gffread ≥ 99% of lyric` and `result ≥ 99% of gffread` (near-1:1 expected).
      Known quirks (accepted, 2026-09-22): gffread emits broken ~70 bp fragments that
      recode drops (Babesia 99.5%); SPLIT .pN parts can give more result rows than
      unique transcripts. A shortfall beyond that → suspect scaffold-name mismatch in
      the GTF↔FASTA join (`main_protists.nf:164` `combine(by: 0)`); investigate via §7
      (work dir is still there on failure). Also: transcripts longer than the recode
      `--limit` are dropped by design (protists 100,000 bp ≈ none; fixed 43e4aae).
      Caveat: the job's automated lyric count is a whitelist of transcript-level
      feature types (mRNA/transcript/RNA/rRNA/tRNA/scRNA/circRNA/ncRNA/lncRNA) — if a
      species' GFF uses another transcript-level type the lyric count will be low;
      recount manually with the §5.1a/§8 feature-type inspection.
- [ ] `README.txt` present; genome length + gene count look sane
- [ ] Analysis outputs present: `<sp>_score_diff.png`, `<sp>_summary.csv`,
      `<sp>_candidates.csv`, `<sp>_result_with_diff.csv`
- [ ] Spot-check `head -3` of the result: columns `og_score_*`, `re_score_*`,
      `TGA_site_score`, `all_secis_*`, `filtered_secis_*` present

If any check fails → §7 (failure handling). The work dir survives only on failure, so
debug from log + `nextflow log <run-uuid>` (layout notes below) + the kept work dir.

### 5.4 Record
Fill the row for `<sp>` in `run_tracker.tsv` (status, job_id, submitted/finished,
result_rows, lyric_transcripts, gffread_transcripts, all_secis, filtered_secis,
predicted_overall/all_secis/filtered_secis, candidates, notes). Metrics come from the
analysis outputs (`<sp>_summary.csv`, `<sp>_candidates.csv`). No "submit next species"
step anymore — the whole batch is already queued.

**Nextflow layout notes (still valid, learned 2026-09-21):**
- Run metadata lives in the LAUNCH dir, not the work dir:
  `cd /users/rg/ileahy/git/gitlab/readthrough && nextflow log` lists runs (name + UUID +
  status) — works after the work dir is removed.
- `nextflow log <run-uuid>` on an ACTIVE run fails with a lock error; while running,
  monitor via the `.out` log's progress table instead.
- After completion: `nextflow log <run-uuid> | tail -40` audits process states
  (count `FAILED`/`ERROR` rows). Requires `module load Java` first.
- WARNING — do NOT use the log banner as a completion signal:
  `workflowCompletionMessage()` is called inline in the workflow body
  (`main_protists.nf:213`), so `Pipeline completed successfully!` is printed at
  STARTUP, not on completion.

### 5.5 Legacy manual loop (retired 2026-09-30, kept for reference)
The pre-batch per-species loop was: back up existing outputs → pre-flight double-gz
check → write params/submit → submit → poll every ~15 min → verify (README, result
file, three-way transcript counts taken from the work dir BEFORE cleanup) → run
`analysis/analyse_protist.py` in singularity manually → record tracker row →
`rm -rf /nfs/scratch01/rg/ileahy/nf_work/<sp>` (only after ALL checks passed) → submit
the next species without approval (user preference 2026-09-22). Every step of it now
happens inside the per-species submit job or in §5.2–5.4 above.

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
- **Seblastian rejects fasta identifiers >63 characters** (learned 2026-09-23,
  Cyanidiococcus attempt 1, 28725643): SECISSEARCH runs `Seblastian.py -t <transcript
  fasta>` and dies if any header's first word is >63 chars — GenBank `gnl|WGS` ids in
  the merged gidRef GFFs can be 64-65 chars (2 of 5,189 for Cyanidiococcus). Fixed in
  760ab22: GFFREAD/GFFREAD_CHR shorten any header >63 chars to first 55 chars + '_' +
  6-digit rolling hash of the full id (max 62, deterministic, collision-safe in
  practice); every downstream stage (recoding, split, geneid, secissearch, logos) sees
  the same shortened names, so joins stay consistent. Conticribra checked: all 14,490
  ids are 17 chars, unaffected.
- **Scaffold batching (2026-09-28, verified by differential toy test)**: the 6
  per-scaffold stages of `main_protists.nf` now run as `scaffold_batches` (default
  24) hash-bucket jobs via `modules/scaffold_batches.nf` + batched `GFFREAD_CHR`
  (`modules/gffread_chr.nf`); mammal/model pipelines untouched. NF 26.04.6 rules the
  implementation depends on (all verified empirically, 2026-09-28 toy runs):
  1. A SINGLE list `path` input without `stageAs` stages flat into the task work-dir
     root, base names preserved (production `EXTRACT_SEQUENCE_LOGOS` pattern) — the
     batch scripts loop over root globs.
  2. `stageAs` on a list input stages files as `in/1.ext`, `in/2.ext`, ... —
     NUMBERED, names lost. Never use it here.
  3. `groupTuple()` is COLUMN-WISE: for (a,b,c) tuples it emits (a, [b...], [c...]) —
     NOT (a, [(b,c),...]). Indexing a File/Path object with `it[0]`/`it[1]` returns
     its PATH COMPONENTS (getName(0)/getName(1)), which is how the "input file name
     collision: nfs, scratch01" errors arose (channels poisoned with literal
     `nfs`/`scratch01` paths — it was never a staging bug).
  4. Every batched process emits the SAME per-scaffold output names as the original
     one-job-per-scaffold module, so geneid part names and the
     SELECT_INTERESTING/GET_ORIGINAL_PREDICTIONS replaceAll id chains are unaffected.
  5. GFFREAD_CHR receives its bucket's gtf+fasta files as ONE merged list and
     re-pairs by filename (`*.part_<scaffold>.fa`) inside the job.
  Toy test (30 scaffolds / 60 transcripts, K=4): batched run 157 tasks vs 313
  unbatched, final `toy_batch_ORFsearch_SECIS.result` byte-identical (sorted diff).
  First production use: Dunaliella_salina (species 11).

## 9. Status (update every session)

- Phase: **2 — batch submitted, monitoring** (since 2026-09-30): ALL 35 remaining
  species are queued (jobs 29099233–29099267). Per-species jobs self-run pipeline →
  verify → analyse → cleanup; I poll the batch and close out each species as its job
  completes (§5.2–5.4). No more manual per-species triggering.
- Deferred→resubmitted: species 5 — Conticribra_weissflogii, job 29099233 (2026-09-30
  10:15) with `-resume` — reuses the attempt-2 partial cache in work dir
  `/nfs/scratch01/rg/ileahy/nf_work/Conticribra_weissflogii` (attempt 2, 28713877, was
  scancelled 2026-09-23 12:35 at user decision).
- Completed: 11 / 46
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
  - Cyanidiococcus_yangmingshanensis (job 28726227, attempt 2 after the 63-char
    header fix 760ab22): lyric 5,189 → gffread 5,189 → result 5,189 unique (100%
    three-way match; GenBank reference annotation, 20 scaffolds, 7m19s with -resume);
    1 candidate rna-gnl|WGS:VWRR|F1559_005062-T1_mrna (score_diff 0.84);
    old Apr-08 output backed up ORFsearch_old_20260923_131019
  - Cyanidioschyzon_merolae_strain_10D (job 28728435): lyric 5,373 → gffread 5,373
    → result 5,373 unique (100% three-way match; RefSeq reference annotation
    4803 mRNA + 521 transcript + 12 rRNA + 37 tRNA, 20 scaffolds, 36m49s);
    2 candidates, both negative score_diff (weak signal): rna-XM_005534776.1
    (-1.21), rna-XR_002461538.1 (-0.78); old Apr run backed up
    ORFsearch_old_20260923_134315
  - Cyanidium_caldarium (job 28730609): lyric 4,870 → gffread 4,870 → result
    4,870 unique (100% three-way match; GenBank reference annotation, 20
    scaffolds, 5m34s); 13 SECIS found, 0 survived the filter → 0 candidates
  - Cylindrotheca_closterium (job 28732205): lyric 67,748 → result 67,748
    unique (100% of lyric; gffread count not taken — workdir audit declined
    09-28, workdir kept); 4h12m; 1,015 SECIS found, 14 survived the filter;
    2 candidate transcripts (gene-/rna-CYCCA115_LOCUS4366, score_diff 2.44)
  - Dunaliella_salina (job 28945196): lyric 48,450 → gffread 48,450 → result
    48,450 unique (100% three-way match; 48,452 data rows — 2 transcripts
    carry 2 rows each); first production run with scaffold_batches=24, 2h12m;
    404 SECIS found, 9 survived the filter; 3 candidates, top agat-rna-6766
    (score_diff 5.19); work dir kept (first production scaffold-batched run)
  - Eimeria_necatrix (job 28964808): COMPLETED 0:0 (1h36m, 2026-09-28 21:55→23:31);
    lyric 180,391 → gffread 180,391 → result 180,391 unique (100% three-way match;
    180,438 data rows); largest transcript count so far; 2,176 SECIS found, 78
    survived the filter; 48 candidate transcripts (predicted + filtered SECIS);
    old Aug-18 output backed up ORFsearch_old_20260928_215109; work dir removed
    2026-09-30
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
    - [x] Cyanidiococcus attempt 1 (28725643) FAILED 13:14 at SECISSEARCH: Seblastian
          rejects fasta ids >63 chars (2 GenBank gnl|WGS ids were 64-65 chars). Fixed
          at the GFFREAD choke point 760ab22 (headers >63 → 55 chars + 6-digit rolling
          hash, max 62; verified 5189/5189 unique + deterministic). Attempt 2 =
          28726227 (13:20, -resume reuses UNZIP/AGAT/FILTER tasks).
- Done this session (2026-09-23):
  - [x] Cyanidiococcus closed out (attempt 2, 28726227, -resume, 7m19s): 100%
        three-way 5,189; 1 candidate rna-gnl|WGS:VWRR|F1559_005062-T1_mrna
        (score_diff 0.84); work dir removed. Batch 6/46.
  - [x] Species 8 Cyanidioschyzon_merolae_strain_10D submitted: job 28728435
        (13:45). Pre-flight: gidRef was DOUBLE-GZ (4th occurrence) — fixed in
        place, orig kept *_doublegz_backup_20260923; RefSeq reference annotation,
        5,373 transcript-level features (4803 mRNA + 521 transcript + 12 rRNA +
        37 tRNA), all 18-char ids (no >63 risk); 20 scaffolds; old Apr run backed
        up ORFsearch_old_20260923_134315.
  - [x] Species 8 closed out: 28728435 COMPLETED 0:0 (36m49s); 100% three-way
        5,373; 2 candidates both negative score_diff (weak signal); work dir
        removed. Batch 7/46.
  - [x] Species 9 Cyanidium_caldarium submitted: job 28730609 (14:41).
        Pre-flight: gidRef DOUBLE-GZ (5th occurrence) — fixed in place, orig kept
        *_doublegz_backup_20260923; GenBank reference annotation, 4,870 mRNA
        (all 33-char ids, no >63 risk); 44 orphan CDS lines for
        FILTER_ORPHAN_CDS; 20 scaffolds; no pre-existing output dir.
  - [x] Species 9 Cyanidium_caldarium closed out: 28730609 COMPLETED 0:0 (5m34s);
        100% three-way 4,870; 13 SECIS found, 0 survived the filter → 0
        candidates; work dir removed. Batch 8/46.
  - [x] Species 10 Cylindrotheca_closterium submitted: job 28732205 (15:10).
        Pre-flight: all inputs single-gz (no double-gz); lyricMerged gidRef =
        LyRic+reference merge (34,836 AGAT RNA + 8,536 AGAT mRNA + 24,187 EMBL
        mRNA + 33 rRNA + 156 tRNA = 67,748 transcript-level; 18,788 transdecoder
        CDS all with valid parents — 0 orphans); local pre-flight chain (orphan
        filter → attr-cleanup awk → AGAT GFF2GTF → gffread) = 67,748 unique
        transcripts, 0 warnings → baseline 67,748; all ids ≤26 chars; 2,534
        scaffolds (largest per-scaffold job count so far); no pre-existing output.
- Done this session (2026-09-28):
  - [x] Cylindrotheca_closterium closed out (28732205 had COMPLETED 0:0 on
        09-23, 4h12m): result 67,758 data rows / 67,748 unique (100% of
        lyric); 1,015 SECIS found, 14 survived the filter; 2 candidate
        transcripts (gene-/rna-CYCCA115_LOCUS4366, score_diff 2.44; candidates
        CSV has 4 rows = 2 transcripts x 2 all-SECIS rows each); gffread count
        not taken (workdir audit declined) — workdir still on /nfs/scratch01.
  - [x] Dunaliella_salina closed out (28945196, 16:35→18:47, 2h12m): 100%
        three-way lyric=gffread=result=48,450 (48,452 data rows); 404 SECIS
        found, 9 survived the filter; 3 candidates, top agat-rna-6766
        (score_diff 5.19); first production run with scaffold_batches=24 —
        worked cleanly; workdir KEPT (user declined audit/cleanup).
  - [x] Species 12 Eimeria_necatrix submitted: job 28964808 (21:55).
        Pre-flight: all 3 .gz single-compressed (no double-gz); merged gidRef
        180,391 transcript-level (63,554 mRNA + 116,837 RNA, no tRNA) +
        328,902 CDS lines (FILTER_ORPHAN_CDS applies); all IDs ≤19 chars (no
        >63 risk); 3,707 scaffolds; old Aug-18 output backed up
        ORFsearch_old_20260928_215109. Largest transcript count so far
        (2.1x Chaetoceros) — long run expected.
- Done this session (2026-09-30, login node):
  - [x] Workflow changed to **batch format** (user decision): per-species
        `runs/submit_<sp>.sh` now does backup → pipeline → result sanity check →
        transcript counts → analysis (singularity) → `rm -rf` work dir on success;
        `--mem 16G` (4G OOM-killed Eimeria's 180k-row analysis on the login node).
  - [x] Eimeria_necatrix closed out (28964808, 100% three-way 180,391); its
        analysis re-run as SLURM job 29099232 (16G) — COMPLETED 0:0 09-30, all 4
        analysis files written, tracker row filled; work dir removed.
  - [x] Generated params+submit for all 35 remaining species (all 140 input files
        verified to exist) and submitted them all at once via
        `runs/submit_all_remaining.sh` (2026-09-30 10:15, jobs 29099233–29099267);
        Conticribra with `-resume`.
  - [x] SKILL.md restructured for the batch format (§1, §4, §5, §9).
- Next:
  - [ ] Monitor the 35 batch jobs 29099233–29099267:
        `squeue -u ileahy --format="%.10i %.2t %.15j %.20T"`;
        per completed job → §5.3 verify → §5.4 record tracker row
  - [ ] On any FAILED job: work dir kept → §7 diagnostics (tail .out/.err,
        `nextflow log <run-uuid>`), fix, re-submit that species only
- Notes: this session ran **directly on the cluster login node** (genoa64-05, user
  ileahy) — no `ssh login` prefix needed; from the local machine use the `ssh login`
  forms as written. Re-read this file at the start of each session and keep this
  section current.
- Session 2026-09-23 (genoa64-04): code fix e463458 was committed directly in the
  CLUSTER clone — the local Mac repo
  (`/Users/iseult/gitlab/SECIS_independent/CascadeProjects/windsurf-project`) needs a
  `git pull` to pick it up.
