"""Primary analysis runner — hypotheses H1 (superiority) and H2 (non-inferiority).

Reads two predictor score TSVs, restricts to Stratum 2 (n=52 Joynt binary),
runs paired DeLong per Sun & Xu 2014, applies the pre-registered decision
rules, and emits both a machine-readable results TSV and a plain-language
decision block.

Pre-registered decision rules (frozen at OSF 2026-09-13):

- H1 primary (two-sided superiority, α=0.05): predictor A is superior to
  predictor B if AUC_A > AUC_B AND the two-sided paired-DeLong p-value < 0.05.
  Rejecting H0 in either direction is publishable.

- H2 secondary (non-inferiority, margin δ=0.05 on AUC scale): predictor A is
  non-inferior to predictor B if the lower bound of the 95% CI on the paired
  difference (AUC_A - AUC_B) exceeds -0.05.

- Stratified-bootstrap AUC CIs for each predictor separately are reported as
  an exploratory robustness display alongside the pre-registered results.
  See docs/protocol/WRITEUP_NOTES.md.

Usage:

    python -m src.analyse.run_primary_analysis \\
        --variants data/variant_list_frozen.tsv \\
        --scores-a results/predictions/alphagenome.tsv --label-a AlphaGenome \\
        --scores-b results/predictions/spliceai.tsv --label-b SpliceAI \\
        --score-column score \\
        --output results/primary_analysis.tsv

The predictor score TSV must contain (variant_id, <--score-column>).

Reads only frozen inputs plus the two predictor score files. Produces derived
outputs only. No API calls.
"""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

from src.analyse.io import (
    build_stratum2_frame,
    load_predictor_scores,
    load_variants,
)
from src.analyse.paired_delong import (
    non_inferiority_check,
    paired_delong_test,
    superiority_check,
)

BOOT_ITERS = 10_000  # stratified bootstrap iterations for the exploratory CI
BOOT_SEED = 20_260_913  # OSF registration date, frozen


@dataclass(frozen=True)
class BootstrapCI:
    auc_median: float
    ci_low_2p5: float
    ci_high_97p5: float
    n_iters: int
    seed: int


def stratified_bootstrap_auc_ci(
    y: np.ndarray,
    scores: np.ndarray,
    iters: int = BOOT_ITERS,
    seed: int = BOOT_SEED,
) -> BootstrapCI:
    """Class-stratified bootstrap CI on a single-predictor AUC.

    Preserves the empirical class balance in each resample (both classes are
    resampled with replacement of size n_pos and n_neg independently) so
    ROC-AUC is always defined per iteration. Exploratory robustness display,
    not a pre-registered test.
    """
    y = np.asarray(y, dtype=int)
    scores = np.asarray(scores, dtype=float)
    pos_idx = np.where(y == 1)[0]
    neg_idx = np.where(y == 0)[0]
    if pos_idx.size == 0 or neg_idx.size == 0:
        raise ValueError("stratified bootstrap needs both classes present")

    rng = np.random.default_rng(seed)
    aucs = np.empty(iters, dtype=float)
    for i in range(iters):
        p = rng.choice(pos_idx, size=pos_idx.size, replace=True)
        n = rng.choice(neg_idx, size=neg_idx.size, replace=True)
        idx = np.concatenate([p, n])
        aucs[i] = roc_auc_score(y[idx], scores[idx])
    return BootstrapCI(
        auc_median=float(np.median(aucs)),
        ci_low_2p5=float(np.percentile(aucs, 2.5)),
        ci_high_97p5=float(np.percentile(aucs, 97.5)),
        n_iters=iters,
        seed=seed,
    )


def format_decision_block(
    label_a: str,
    label_b: str,
    delong: dict,
    sup: dict,
    non_inf: dict,
    boot_a: BootstrapCI,
    boot_b: BootstrapCI,
    ni_margin: float,
    coverage: dict,
) -> str:
    """Human-readable plain-language summary of the frozen decision rules."""
    lines: list[str] = []
    lines.append(f"# Primary analysis — {label_a} vs {label_b}")
    lines.append("")
    lines.append("## Stratum 2 coverage")
    lines.append(f"- Frozen Stratum 2 variants: {coverage['n_stratum2']}")
    lines.append(f"- Scored by both predictors: {coverage['n_paired']}")
    lines.append(f"- Missing from {label_a}: {coverage['n_missing_a']}")
    lines.append(f"- Missing from {label_b}: {coverage['n_missing_b']}")
    lines.append(f"- Class balance in analysed set: {coverage['n_pos']} splice-affecting / {coverage['n_neg']} no-effect")
    lines.append("")
    lines.append("## Pre-registered paired DeLong (H1 primary, two-sided, α=0.05)")
    lines.append(f"- AUC {label_a}: {delong['auc_a']:.4f}")
    lines.append(f"- AUC {label_b}: {delong['auc_b']:.4f}")
    lines.append(f"- Paired difference (AUC_{label_a} − AUC_{label_b}): {delong['auc_diff']:+.4f}")
    lines.append(f"- Standard error: {delong['std_error']:.4f}")
    lines.append(f"- 95% CI on paired difference: [{delong['ci_lower']:+.4f}, {delong['ci_upper']:+.4f}]")
    lines.append(f"- z: {delong['z']:.3f}")
    lines.append(f"- Two-sided p: {delong['p_value']:.4g}")
    lines.append(f"- H1 decision: **{sup['decision'].upper()}** ({label_a} vs {label_b})")
    lines.append("")
    lines.append(f"## Pre-registered non-inferiority (H2 secondary, margin δ={ni_margin:.2f})")
    lines.append(f"- Lower bound of 95% CI on paired difference: {non_inf['ci_lower']:+.4f}")
    lines.append(f"- Non-inferiority margin: −{non_inf['margin']:.2f}")
    lines.append(f"- H2 decision: **{non_inf['decision'].upper()}**")
    lines.append("")
    lines.append("## Exploratory stratified-bootstrap AUC CIs (robustness display)")
    lines.append(f"- {label_a}: median {boot_a.auc_median:.4f}, 95% CI [{boot_a.ci_low_2p5:.4f}, {boot_a.ci_high_97p5:.4f}]")
    lines.append(f"- {label_b}: median {boot_b.auc_median:.4f}, 95% CI [{boot_b.ci_low_2p5:.4f}, {boot_b.ci_high_97p5:.4f}]")
    lines.append(f"- {boot_a.n_iters} iterations, seed {boot_a.seed} (OSF registration date). Stratified by outcome class.")
    lines.append("")
    lines.append("## Provenance")
    lines.append("- Test implementation: src.analyse.paired_delong (Sun & Xu 2014 fast DeLong).")
    lines.append("- Frozen inputs: data/variant_list_frozen.tsv (hash-locked via data/HASH_MANIFEST.tsv).")
    lines.append("- Pre-registration: OSF DOI 10.17605/OSF.IO/6PGX8, registered 2026-09-13.")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=Path, required=True)
    parser.add_argument("--scores-a", type=Path, required=True)
    parser.add_argument("--label-a", type=str, required=True)
    parser.add_argument("--scores-b", type=Path, required=True)
    parser.add_argument("--label-b", type=str, required=True)
    parser.add_argument("--score-column", type=str, default="score",
                        help="Score column name in both predictor TSVs (default: score)")
    parser.add_argument("--score-column-a", type=str, default=None,
                        help="Override for predictor A's score column (defaults to --score-column)")
    parser.add_argument("--score-column-b", type=str, default=None,
                        help="Override for predictor B's score column (defaults to --score-column)")
    parser.add_argument("--ni-margin", type=float, default=0.05,
                        help="Pre-registered non-inferiority margin (default 0.05)")
    parser.add_argument("--alpha", type=float, default=0.05,
                        help="Two-sided significance level (default 0.05)")
    parser.add_argument("--boot-iters", type=int, default=BOOT_ITERS)
    parser.add_argument("--boot-seed", type=int, default=BOOT_SEED)
    parser.add_argument("--output", type=Path, required=True,
                        help="Machine-readable results TSV path")
    args = parser.parse_args()

    col_a = args.score_column_a or args.score_column
    col_b = args.score_column_b or args.score_column

    variants = load_variants(args.variants)
    scores_a = load_predictor_scores(args.scores_a, col_a)
    scores_b = load_predictor_scores(args.scores_b, col_b)

    frame = build_stratum2_frame(variants, scores_a, scores_b)

    n_stratum2 = len(frame)
    n_missing_a = int(frame["score_a"].isna().sum())
    n_missing_b = int(frame["score_b"].isna().sum())
    analysed = frame.dropna(subset=["score_a", "score_b"]).reset_index(drop=True)
    n_paired = len(analysed)
    if n_paired < 10:
        raise SystemExit(
            f"Refusing to run: only {n_paired} Stratum 2 variants scored by both predictors "
            f"(missing from {args.label_a}: {n_missing_a}; from {args.label_b}: {n_missing_b})."
        )

    y = analysed["binary_splice_affecting"].to_numpy(int)
    a = analysed["score_a"].to_numpy(float)
    b = analysed["score_b"].to_numpy(float)

    result = paired_delong_test(y, a, b, alpha=args.alpha)
    sup = superiority_check(result, alpha=args.alpha)
    non_inf = non_inferiority_check(result, margin=args.ni_margin)

    boot_a = stratified_bootstrap_auc_ci(y, a, iters=args.boot_iters, seed=args.boot_seed)
    boot_b = stratified_bootstrap_auc_ci(y, b, iters=args.boot_iters, seed=args.boot_seed)

    coverage = {
        "n_stratum2": n_stratum2,
        "n_paired": n_paired,
        "n_missing_a": n_missing_a,
        "n_missing_b": n_missing_b,
        "n_pos": int((y == 1).sum()),
        "n_neg": int((y == 0).sum()),
    }

    # Machine-readable TSV
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["label_a", args.label_a])
        w.writerow(["label_b", args.label_b])
        for k, v in coverage.items():
            w.writerow([k, v])
        for k, v in result.as_dict().items():
            w.writerow([f"delong_{k}", v])
        for k, v in sup.items():
            w.writerow([f"h1_{k}", v])
        for k, v in non_inf.items():
            w.writerow([f"h2_{k}", v])
        for name, ci in (("a", boot_a), ("b", boot_b)):
            w.writerow([f"boot_{name}_auc_median", ci.auc_median])
            w.writerow([f"boot_{name}_ci_low_2p5", ci.ci_low_2p5])
            w.writerow([f"boot_{name}_ci_high_97p5", ci.ci_high_97p5])
        w.writerow(["boot_iters", args.boot_iters])
        w.writerow(["boot_seed", args.boot_seed])

    # Sidecar JSON with full nested structure
    json_path = args.output.with_suffix(".json")
    payload = {
        "label_a": args.label_a,
        "label_b": args.label_b,
        "coverage": coverage,
        "paired_delong": result.as_dict(),
        "h1_superiority": sup,
        "h2_non_inferiority": non_inf,
        "bootstrap": {
            "iters": args.boot_iters,
            "seed": args.boot_seed,
            "a": {
                "auc_median": boot_a.auc_median,
                "ci_low_2p5": boot_a.ci_low_2p5,
                "ci_high_97p5": boot_a.ci_high_97p5,
            },
            "b": {
                "auc_median": boot_b.auc_median,
                "ci_low_2p5": boot_b.ci_low_2p5,
                "ci_high_97p5": boot_b.ci_high_97p5,
            },
        },
        "provenance": {
            "test_implementation": "src.analyse.paired_delong (Sun & Xu 2014)",
            "variant_list": str(args.variants),
            "osf_doi": "10.17605/OSF.IO/6PGX8",
        },
    }
    json_path.write_text(json.dumps(payload, indent=2))

    # Plain-language decision block
    md_path = args.output.with_suffix(".md")
    md_path.write_text(format_decision_block(
        args.label_a, args.label_b,
        result.as_dict(), sup, non_inf, boot_a, boot_b,
        args.ni_margin, coverage,
    ))

    print(format_decision_block(
        args.label_a, args.label_b,
        result.as_dict(), sup, non_inf, boot_a, boot_b,
        args.ni_margin, coverage,
    ))
    print(f"Wrote {args.output}, {json_path}, {md_path}")


if __name__ == "__main__":
    main()
