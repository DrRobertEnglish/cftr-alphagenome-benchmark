"""Negative-control permutation harness.

Estimates the null distribution of ROC-AUC for a splice-effect predictor by
repeatedly shuffling the binary labels of the frozen variant list. Because
the labels are randomly permuted, any consistent lift above 0.5 is a
false-positive signal — the harness is the reference for what "no better
than chance" looks like on this specific dataset (with its specific class
balance and any structural quirks) rather than on an idealised balanced set.

Pre-registered in the OSF protocol as the negative control against which
any real predictor must be compared. Fixed seed, 1000 iterations, reported
median and 95% interval, one-sided p-value for "real AUC > null 97.5th
percentile" available if requested downstream.

Usage:

    python -m src.controls.label_permutation_null \
        --variants data/variant_list_frozen.tsv \
        --scores results/controls/splice_site_distance_scores.tsv \
        --iters 1000 --seed 20260913 \
        --output results/controls/permutation_null.tsv

The predictor's real score column is read from --scores; the harness
shuffles the label column of --variants and computes an ROC-AUC each
iteration. Output is one row per iteration with (iteration, auc_null).

Reads only frozen inputs; produces one derived output. No API calls.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


def load_variants(variants_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(variants_tsv, sep="\t")
    required = {"variant_id", "stratum", "binary_splice_affecting"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"variant TSV missing required columns: {missing}")
    return df


def load_scores(scores_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(scores_tsv, sep="\t")
    required = {"variant_id", "score_splice_site_distance"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"scores TSV missing required columns: {missing}")
    return df


def run_permutation_null(
    variants: pd.DataFrame,
    scores: pd.DataFrame,
    iters: int,
    seed: int,
    stratum_filter: str | None = "2",
) -> pd.DataFrame:
    # Only merge the score column from the scores frame to avoid duplicated
    # metadata columns (variants TSV is the source of truth for label + stratum).
    scores_slim = scores[["variant_id", "score_splice_site_distance"]]
    merged = variants.merge(scores_slim, on="variant_id", how="inner")
    if stratum_filter is not None:
        merged = merged[merged["stratum"].astype(str) == str(stratum_filter)]
    if merged.empty:
        raise ValueError("no rows after merging variants with scores + stratum filter")
    labels = merged["binary_splice_affecting"].to_numpy(dtype=int)
    if len(np.unique(labels)) < 2:
        raise ValueError("both classes must be present in the merged set for ROC-AUC")
    predictions = merged["score_splice_site_distance"].to_numpy(dtype=float)

    rng = np.random.default_rng(seed)
    aucs = np.empty(iters, dtype=float)
    perm_labels = labels.copy()
    for i in range(iters):
        rng.shuffle(perm_labels)
        # If shuffling produced a single-class draw (rare when class balance
        # is highly skewed), redraw once. AUC is undefined for one class.
        while len(np.unique(perm_labels)) < 2:
            rng.shuffle(perm_labels)
        aucs[i] = roc_auc_score(perm_labels, predictions)

    return pd.DataFrame({
        "iteration": np.arange(1, iters + 1),
        "auc_null": aucs,
    })


def summarise(null_aucs: pd.DataFrame) -> dict:
    a = null_aucs["auc_null"].to_numpy()
    return {
        "n_iterations": len(a),
        "mean_auc": float(np.mean(a)),
        "median_auc": float(np.median(a)),
        "ci_low_2p5": float(np.percentile(a, 2.5)),
        "ci_high_97p5": float(np.percentile(a, 97.5)),
        "sd_auc": float(np.std(a, ddof=1)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=Path, required=True)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--iters", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260913,
                        help="frozen seed = OSF registration date")
    parser.add_argument("--stratum", type=str, default="2",
                        help="stratum to run the null on (2 = Joynt binary)")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    variants = load_variants(args.variants)
    scores = load_scores(args.scores)
    null_aucs = run_permutation_null(
        variants, scores, iters=args.iters, seed=args.seed,
        stratum_filter=args.stratum,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    null_aucs.to_csv(args.output, sep="\t", index=False)

    summary = summarise(null_aucs)
    summary_path = args.output.with_suffix(".summary.tsv")
    with summary_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        for k, v in summary.items():
            w.writerow([k, v])
    print(f"Wrote {args.output} ({args.iters} iterations) and {summary_path}")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
