"""Generate publication figures for the CFTR AlphaGenome benchmark.

Figure 1 — Stratum 2 (n=52) ROC curves for AlphaGenome, SpliceAI, Pangolin,
with bootstrap 95% CI bands and AUCs annotated.

Figure 2 — Stratum 1 (n=8) scatter of predictor score vs pct_normally_spliced
for each predictor, with rank-based fit line and Spearman rho annotated.

Reads the frozen variant list + hash-locked prediction TSVs; writes PNGs to
`figures/`. No random state used except the pre-registered bootstrap seed
20260913 (mirrors the analysis runners).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import auc, roc_curve

from src.analyse.io import (
    build_stratum1_frame,
    load_predictor_scores,
    load_variants,
)

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
PRED = REPO / "results" / "predictions"
FIG = REPO / "figures"
FIG.mkdir(parents=True, exist_ok=True)

BOOT_ITERS = 10_000
BOOT_SEED = 20260913

PREDICTORS = [
    ("AlphaGenome", "alphagenome.tsv", "#1f77b4"),
    ("SpliceAI", "spliceai.tsv", "#d62728"),
    ("Pangolin", "pangolin.tsv", "#2ca02c"),
]


def stratified_bootstrap_roc(y_true: np.ndarray, y_score: np.ndarray, n_iters: int, seed: int):
    """Stratified (by outcome) bootstrap; return list of (fpr, tpr) resampled curves and AUCs."""
    rng = np.random.default_rng(seed)
    pos_idx = np.flatnonzero(y_true == 1)
    neg_idx = np.flatnonzero(y_true == 0)
    aucs = []
    # Compute TPR at a common FPR grid for band construction
    fpr_grid = np.linspace(0.0, 1.0, 101)
    tprs = []
    for _ in range(n_iters):
        p = rng.choice(pos_idx, size=len(pos_idx), replace=True)
        n = rng.choice(neg_idx, size=len(neg_idx), replace=True)
        boot = np.concatenate([p, n])
        yt = y_true[boot]
        ys = y_score[boot]
        if len(np.unique(yt)) < 2:
            continue
        fpr, tpr, _ = roc_curve(yt, ys)
        aucs.append(auc(fpr, tpr))
        # Interpolate to common grid
        tpr_interp = np.interp(fpr_grid, fpr, tpr)
        tpr_interp[0] = 0.0
        tprs.append(tpr_interp)
    return fpr_grid, np.array(tprs), np.array(aucs)


def figure1_roc() -> Path:
    variants = load_variants(DATA / "variant_list_frozen.tsv")
    fig, ax = plt.subplots(figsize=(6.5, 6.0), dpi=200)

    for label, filename, color in PREDICTORS:
        scores = load_predictor_scores(PRED / filename, score_column="score")
        s2 = variants[variants["stratum"] == "2"].copy()
        frame = s2.merge(scores, on="variant_id", how="left").dropna(subset=["score"])
        y_true = frame["binary_splice_affecting"].astype(int).values
        y_score = frame["score"].values

        # Point-estimate ROC
        fpr, tpr, _ = roc_curve(y_true, y_score)
        auc_point = auc(fpr, tpr)

        # Bootstrap band
        fpr_grid, tprs, aucs = stratified_bootstrap_roc(y_true, y_score, BOOT_ITERS, BOOT_SEED)
        tpr_lo = np.percentile(tprs, 2.5, axis=0)
        tpr_hi = np.percentile(tprs, 97.5, axis=0)
        auc_lo, auc_hi = np.percentile(aucs, [2.5, 97.5])

        ax.fill_between(fpr_grid, tpr_lo, tpr_hi, color=color, alpha=0.12, linewidth=0)
        ax.plot(fpr, tpr, color=color, lw=2.2,
                label=f"{label}  AUC {auc_point:.3f}  [{auc_lo:.3f}, {auc_hi:.3f}]")

    ax.plot([0, 1], [0, 1], color="grey", ls="--", lw=1.0, label="Chance")
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.02)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("Stratum 2 (n=52) — ROC on binary splice-affecting label\n"
                 "shaded bands: 10,000-iter stratified bootstrap 95% CI",
                 fontsize=11)
    ax.legend(loc="lower right", fontsize=9, frameon=True)
    ax.grid(alpha=0.25, linewidth=0.5)
    fig.tight_layout()
    out = FIG / "figure1_roc_stratum2.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def figure2_stratum1() -> Path:
    variants = load_variants(DATA / "variant_list_frozen.tsv")
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5), dpi=200, sharey=True)

    for ax, (label, filename, color) in zip(axes, PREDICTORS):
        scores = load_predictor_scores(PRED / filename, score_column="score")
        frame = build_stratum1_frame(variants, scores)
        # Keep numeric pct_normally_spliced only
        frame = frame[pd.to_numeric(frame["pct_normally_spliced"], errors="coerce").notna()].copy()
        frame["pct_normally_spliced"] = frame["pct_normally_spliced"].astype(float)
        frame = frame.dropna(subset=["score"])
        x = frame["score"].values
        y = frame["pct_normally_spliced"].values

        rho, _ = spearmanr(x, y)
        # Bootstrap CI on rho
        rng = np.random.default_rng(BOOT_SEED)
        rhos = []
        n = len(x)
        for _ in range(BOOT_ITERS):
            idx = rng.integers(0, n, size=n)
            xb, yb = x[idx], y[idx]
            if len(np.unique(xb)) < 2 or len(np.unique(yb)) < 2:
                continue
            r, _ = spearmanr(xb, yb)
            if np.isfinite(r):
                rhos.append(r)
        rho_lo, rho_hi = np.percentile(rhos, [2.5, 97.5])

        ax.scatter(x, y, color=color, s=60, edgecolor="black", linewidth=0.6, zorder=3)
        # Rank-based (Spearman) fit line: fit OLS on ranks, then map back to x for display
        rx = pd.Series(x).rank().values
        ry = pd.Series(y).rank().values
        # OLS on ranks
        m, b = np.polyfit(rx, ry, 1)
        rx_line = np.linspace(rx.min(), rx.max(), 50)
        ry_line = m * rx_line + b
        # Convert rank-line to score/pct space via linear interpolation of the sorted pairs
        # (this is only for visual guidance; the statistic reported is the actual Spearman rho)
        order_x = np.argsort(x)
        order_y = np.argsort(y)
        # Interpolate rank -> value
        x_line = np.interp(rx_line, np.sort(rx), x[order_x])
        y_line = np.interp(ry_line, np.sort(ry), y[order_y])
        ax.plot(x_line, y_line, color=color, lw=1.5, alpha=0.7)

        ax.set_title(f"{label}\nSpearman ρ = {rho:.2f}  [{rho_lo:.2f}, {rho_hi:.2f}]",
                     fontsize=11)
        ax.set_xlabel(f"{label} score")
        ax.grid(alpha=0.25, linewidth=0.5)

    axes[0].set_ylabel("% normally spliced transcript")
    fig.suptitle("Stratum 1 (n=8) — predictor score vs quantitative residual splicing\n"
                 "pre-specified expected direction: NEGATIVE Spearman ρ",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    out = FIG / "figure2_stratum1_scatter.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    p1 = figure1_roc()
    print(f"Wrote {p1.relative_to(REPO)}")
    p2 = figure2_stratum1()
    print(f"Wrote {p2.relative_to(REPO)}")


if __name__ == "__main__":
    main()
