"""Stratum 1 analysis — hypothesis H3 (descriptive rank correlation).

Reads one predictor score TSV, restricts to Stratum 1 variants with a
numeric % normally spliced measurement, computes Spearman ρ and a 10 000-iter
percentile bootstrap CI, and emits both a machine-readable results TSV and a
plain-language summary block.

Pre-registered decision rule (frozen at OSF 2026-09-13):

- H3 is descriptive. Report Spearman ρ between predictor score and measured
  % normally spliced, with a 10 000-iter bootstrap 95% CI. No p-value is
  reported. No formal decision rule. This exists to characterise whether a
  predictor's continuous score tracks the quantitative ground truth in
  Stratum 1, independent of the binary Stratum 2 comparison.

The Masvidal complex allele c.[2657+5G>A;2562T>G] is not in the frozen TSV
by design — see docs/protocol/WRITEUP_NOTES.md.

Usage:

    python -m src.analyse.run_stratum1 \\
        --variants data/variant_list_frozen.tsv \\
        --scores results/predictions/alphagenome.tsv \\
        --label AlphaGenome --score-column score \\
        --output results/stratum1_alphagenome.tsv

Reads only frozen inputs plus the predictor score file. Produces derived
outputs only. No API calls.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy import stats

from src.analyse.io import (
    build_stratum1_frame,
    load_predictor_scores,
    load_variants,
)

BOOT_ITERS = 10_000  # pre-registered
BOOT_SEED = 20_260_913  # OSF registration date, frozen


def bootstrap_spearman_ci(
    x: np.ndarray,
    y: np.ndarray,
    iters: int = BOOT_ITERS,
    seed: int = BOOT_SEED,
) -> dict:
    """Percentile bootstrap CI on Spearman ρ."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape or x.ndim != 1:
        raise ValueError("x and y must be 1-D arrays of equal length")
    if x.size < 4:
        raise ValueError(f"need at least 4 points for a bootstrap CI on rho; got {x.size}")
    n = x.size
    rng = np.random.default_rng(seed)
    rhos = np.empty(iters, dtype=float)
    for i in range(iters):
        idx = rng.integers(0, n, size=n)
        # Resample can produce zero-variance draws; catch and mark NaN
        xi = x[idx]
        yi = y[idx]
        if np.all(xi == xi[0]) or np.all(yi == yi[0]):
            rhos[i] = np.nan
            continue
        rhos[i] = stats.spearmanr(xi, yi).statistic
    rhos_clean = rhos[~np.isnan(rhos)]
    if rhos_clean.size < iters * 0.9:
        # Too many zero-variance resamples means the input is degenerate
        raise ValueError(
            f"only {rhos_clean.size}/{iters} bootstrap iterations produced a "
            f"defined rho — inputs are likely degenerate (many ties or constant)"
        )
    return {
        "rho_median": float(np.median(rhos_clean)),
        "ci_low_2p5": float(np.percentile(rhos_clean, 2.5)),
        "ci_high_97p5": float(np.percentile(rhos_clean, 97.5)),
        "n_iters_defined": int(rhos_clean.size),
        "n_iters_total": int(iters),
        "seed": int(seed),
    }


def format_summary(
    label: str,
    rho: float,
    boot: dict,
    n: int,
    variant_ids: list[str],
) -> str:
    lines: list[str] = []
    lines.append(f"# Stratum 1 rank correlation — {label} (H3, descriptive)")
    lines.append("")
    lines.append(f"- Variants included (with numeric % normally spliced and predictor score): n = {n}")
    lines.append(f"- Spearman ρ (point estimate): {rho:+.4f}")
    lines.append(f"- Bootstrap 95% CI on ρ: [{boot['ci_low_2p5']:+.4f}, {boot['ci_high_97p5']:+.4f}]")
    lines.append(f"- Bootstrap median ρ: {boot['rho_median']:+.4f}")
    lines.append(f"- Bootstrap iterations (defined / total): {boot['n_iters_defined']} / {boot['n_iters_total']}, seed {boot['seed']}")
    lines.append("")
    lines.append("## Variants analysed")
    for vid in variant_ids:
        lines.append(f"- {vid}")
    lines.append("")
    lines.append("## Pre-registration note")
    lines.append("H3 is descriptive per the frozen OSF protocol. No p-value reported. The pre-specified expected direction is NEGATIVE Spearman ρ: higher predictor score (more splice-affecting) should track LOWER pct_normally_spliced (less healthy transcript). The correct direction of association is fixed by the scoring convention documented in the predictor's own runner.")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=Path, required=True)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--label", type=str, required=True)
    parser.add_argument("--score-column", type=str, default="score")
    parser.add_argument("--boot-iters", type=int, default=BOOT_ITERS)
    parser.add_argument("--boot-seed", type=int, default=BOOT_SEED)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    variants = load_variants(args.variants)
    scores = load_predictor_scores(args.scores, args.score_column)

    frame = build_stratum1_frame(variants, scores)
    analysed = frame.dropna(subset=["score", "pct_normally_spliced"]).reset_index(drop=True)
    n = len(analysed)
    if n < 4:
        raise SystemExit(f"Refusing to run: only {n} Stratum 1 variants scored (need ≥ 4)")

    x = analysed["score"].to_numpy(float)
    y = analysed["pct_normally_spliced"].to_numpy(float)
    rho = float(stats.spearmanr(x, y).statistic)
    boot = bootstrap_spearman_ci(x, y, iters=args.boot_iters, seed=args.boot_seed)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["label", args.label])
        w.writerow(["n_variants", n])
        w.writerow(["spearman_rho", rho])
        for k, v in boot.items():
            w.writerow([f"boot_{k}", v])

    json_path = args.output.with_suffix(".json")
    json_path.write_text(json.dumps({
        "label": args.label,
        "n_variants": n,
        "spearman_rho": rho,
        "bootstrap": boot,
        "variant_ids": analysed["variant_id"].tolist(),
        "provenance": {
            "variant_list": str(args.variants),
            "osf_doi": "10.17605/OSF.IO/6PGX8",
        },
    }, indent=2))

    md_path = args.output.with_suffix(".md")
    md_path.write_text(format_summary(
        args.label, rho, boot, n, analysed["variant_id"].tolist(),
    ))
    print(format_summary(
        args.label, rho, boot, n, analysed["variant_id"].tolist(),
    ))
    print(f"Wrote {args.output}, {json_path}, {md_path}")


if __name__ == "__main__":
    main()
