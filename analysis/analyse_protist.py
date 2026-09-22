#!/usr/bin/env python3
"""
Per-species analysis of a SecORFsearch result file.

Reproduces the logic of secorfsearch_data/secis_independent/Babesia_duncani.ipynb
(without the MFE merge step).

Usage:
    python analyse_protist.py \
        --result /no_backup/rg/ileahy/<species>/ORFsearch/<species>_ORFsearch_SECIS.result \
        --species <species> \
        --outdir /no_backup/rg/ileahy/<species>/ORFsearch/analysis
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PREDICTED_SCORE_DIFF_MIN = -1.8
PREDICTED_RE_SCORE_MIN = -1.5

# Transcript-count reconciliation threshold (SKILL.md §5.3): near-1:1 expected;
# a small shortfall is tolerated (known gffread ~70 bp fragment quirk, see
# Babesia_duncani: 12,076 / 12,133 = 99.5%).
RECONCILIATION_MIN = 0.99


def summary_block(df: pd.DataFrame) -> pd.DataFrame:
    score_diff = df["re_score_score"] - df["og_score_score"]
    predicted = (score_diff >= PREDICTED_SCORE_DIFF_MIN) & (
        df["re_score_score"] >= PREDICTED_RE_SCORE_MIN
    )
    summary = pd.DataFrame(
        {"is_seleno": [int(predicted.sum()), int((~predicted).sum())]},
        index=["predicted", "not predicted"],
    )
    summary["Total"] = summary.sum(axis=1)
    summary.loc["Total"] = summary.sum()
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument(
        "--gffread-count",
        type=int,
        default=None,
        help="Transcripts extracted by gffread (work dir, pre-cleanup); "
        "used to sanity-check the result transcript count",
    )
    parser.add_argument(
        "--lyric-count",
        type=int,
        default=None,
        help="Transcript-level features in the LyRic GTF/GFF; "
        "used to sanity-check the gffread extraction",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sp = args.species

    df = pd.read_csv(args.result, index_col=0)
    df["score_diff"] = df["re_score_score"] - df["og_score_score"]

    n_total = len(df)
    n_transcripts = int(df.index.nunique())
    n_all_secis = int(df["all_secis_start"].notna().sum())
    n_filtered_secis = int(df["filtered_secis_start"].notna().sum())

    print(f"Species:                 {sp}")
    print(f"Result file:             {args.result}")
    print(f"Total transcript rows:   {n_total}")
    print(f"Unique transcripts:      {n_transcripts}")
    print(f"All SECIS found:         {n_all_secis}")
    print(f"Filtered SECIS found:    {n_filtered_secis}")

    # Transcript reconciliation vs upstream counts (SKILL.md §5.3).
    recon_lines = []
    if args.gffread_count:
        r_out = n_transcripts / args.gffread_count
        flag = "" if r_out >= RECONCILIATION_MIN else "   FLAG < 99%"
        recon_lines.append(
            f"result   {n_transcripts:>7} / gffread {args.gffread_count:>7}"
            f" = {r_out:6.1%}{flag}"
        )
    if args.lyric_count:
        ref = args.gffread_count or n_transcripts
        ref_name = "gffread" if args.gffread_count else "result"
        r_up = ref / args.lyric_count
        flag = "" if r_up >= RECONCILIATION_MIN else "   FLAG < 99%"
        recon_lines.append(
            f"{ref_name}   {ref:>7} / lyric   {args.lyric_count:>7}"
            f" = {r_up:6.1%}{flag}"
        )
    if recon_lines:
        print("\n=== Transcript reconciliation (>= 99% expected) ===")
        for line in recon_lines:
            print(line)

    # score_diff bar plot (rows with both original and recoded scores)
    plot_df = df.dropna(subset=["score_diff"]).copy()
    plot_df = plot_df.sort_values("score_diff", ascending=False).reset_index(drop=True)
    plt.figure(figsize=(24, 12))
    plt.bar(plot_df.index, plot_df["score_diff"], label="Score - CDS", alpha=0.7)
    plt.axhline(0, color="black", linestyle="--")
    plt.ylabel("Difference in Coding Potential")
    plt.title(f"{sp}: Recoded Coding Potential Compared to Original Sequence")
    plt.legend()
    plt.tight_layout()
    png_path = os.path.join(args.outdir, f"{sp}_score_diff.png")
    plt.savefig(png_path, dpi=150)
    plt.close()

    # summary tables: overall / all-SECIS / filtered-SECIS
    groups = {
        "overall": df,
        "all_secis": df[df["all_secis_start"].notna()],
        "filtered_secis": df[df["filtered_secis_start"].notna()],
    }
    rows = []
    for name, sub in groups.items():
        for idx, row in summary_block(sub).iterrows():
            rows.append(
                {
                    "group": name,
                    "is_seleno": idx,
                    "count": int(row["is_seleno"]),
                    "total": int(row["Total"]),
                }
            )
    summary_long = pd.DataFrame(rows)
    summary_path = os.path.join(args.outdir, f"{sp}_summary.csv")
    summary_long.to_csv(summary_path, index=False)

    print("\n=== Summary (long) ===")
    print(summary_long.to_string(index=False))

    for name in ["overall", "all_secis", "filtered_secis"]:
        print(f"\n=== Summary: {name} ===")
        print(summary_block(groups[name]))

    # final candidates: predicted transcripts carrying a filtered SECIS
    score_diff = df["score_diff"]
    predicted = (score_diff >= PREDICTED_SCORE_DIFF_MIN) & (
        df["re_score_score"] >= PREDICTED_RE_SCORE_MIN
    )
    candidates = df[predicted & df["filtered_secis_start"].notna()]
    candidates_path = os.path.join(args.outdir, f"{sp}_candidates.csv")
    candidates.to_csv(candidates_path)
    print(f"\nFinal candidates (predicted + filtered SECIS): {len(candidates)}")
    if len(candidates):
        print(candidates.to_string())

    full_path = os.path.join(args.outdir, f"{sp}_result_with_diff.csv")
    df.to_csv(full_path)

    print(
        f"\nOutputs written to:\n  {png_path}\n  {summary_path}\n  "
        f"{candidates_path}\n  {full_path}"
    )


if __name__ == "__main__":
    main()
