"""Render preprint v0.2 figures: parallel-rule ROC and Stratum-1 scatter.

Overwrites figures/figure1_roc_stratum2.png and figures/figure2_stratum1_scatter.png
with the parallel-rule versions that report AlphaGenome under both the
DeepMind-recommended composite (Rule A) and the OSF-preregistered
per-junction lung/bronchial-context score (Rule B).

- Figure 1: single ROC panel with 4 curves (AlphaGenome Rule A solid, Rule B
  dashed, SpliceAI, Pangolin) with bootstrap 95% CI bands and AUC + CI in the
  legend, matching the v0.1 style.
- Figure 2: 2x2 grid: top row AlphaGenome Rule A (blue) and Rule B (orange);
  bottom row SpliceAI (red) and Pangolin (green), each with Spearman rho +
  bootstrap 95% CI in the panel title.

Deterministic: uses seed 20260913 (the OSF pre-registration date) for every
bootstrap resample, matching the analysis scripts.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_curve, roc_auc_score

REPO = Path(__file__).resolve().parents[1]
FIG_DIR = REPO / "figures"
DATA = REPO / "data" / "variant_list_frozen.tsv"
PRED_DIR = REPO / "results" / "predictions"

SEED = 20260913
N_BOOT = 10_000


def load_stratum2() -> pd.DataFrame:
    variants = pd.read_csv(DATA, sep="\t")
    s2 = variants[variants["stratum"] == 2].copy()
    a = pd.read_csv(PRED_DIR / "alphagenome.tsv", sep="\t")[["variant_id", "score"]].rename(
        columns={"score": "score_A"}
    )
    b = pd.read_csv(PRED_DIR / "alphagenome_prereg.tsv", sep="\t")[
        ["variant_id", "prereg_score"]
    ].rename(columns={"prereg_score": "score_B"})
    sai = pd.read_csv(PRED_DIR / "spliceai.tsv", sep="\t")[["variant_id", "score"]].rename(
        columns={"score": "score_SAI"}
    )
    pgl = pd.read_csv(PRED_DIR / "pangolin.tsv", sep="\t")[["variant_id", "score"]].rename(
        columns={"score": "score_PGL"}
    )
    return (
        s2.merge(a, on="variant_id")
        .merge(b, on="variant_id")
        .merge(sai, on="variant_id")
        .merge(pgl, on="variant_id")
    )


def load_stratum1() -> pd.DataFrame:
    variants = pd.read_csv(DATA, sep="\t")
    s1 = variants[variants["stratum"] == 1].copy()
    s1["pct_normally_spliced"] = pd.to_numeric(s1["pct_normally_spliced"], errors="coerce")
    s1 = s1.dropna(subset=["pct_normally_spliced"])
    a = pd.read_csv(PRED_DIR / "alphagenome.tsv", sep="\t")[["variant_id", "score"]].rename(
        columns={"score": "score_A"}
    )
    b = pd.read_csv(PRED_DIR / "alphagenome_prereg.tsv", sep="\t")[
        ["variant_id", "prereg_score"]
    ].rename(columns={"prereg_score": "score_B"})
    sai = pd.read_csv(PRED_DIR / "spliceai.tsv", sep="\t")[["variant_id", "score"]].rename(
        columns={"score": "score_SAI"}
    )
    pgl = pd.read_csv(PRED_DIR / "pangolin.tsv", sep="\t")[["variant_id", "score"]].rename(
        columns={"score": "score_PGL"}
    )
    return (
        s1.merge(a, on="variant_id")
        .merge(b, on="variant_id")
        .merge(sai, on="variant_id")
        .merge(pgl, on="variant_id")
    )


def stratified_bootstrap_roc(y_true: np.ndarray, scores: np.ndarray, rng: np.random.Generator,
                             n_boot: int = N_BOOT, fpr_grid: np.ndarray | None = None):
    """Return AUC point estimate, AUC bootstrap CI, and TPR-band on fpr_grid."""
    if fpr_grid is None:
        fpr_grid = np.linspace(0.0, 1.0, 201)
    auc_point = roc_auc_score(y_true, scores)
    idx_pos = np.where(y_true == 1)[0]
    idx_neg = np.where(y_true == 0)[0]
    tpr_boots = np.empty((n_boot, len(fpr_grid)))
    auc_boots = np.empty(n_boot)
    for i in range(n_boot):
        rs_pos = rng.choice(idx_pos, size=len(idx_pos), replace=True)
        rs_neg = rng.choice(idx_neg, size=len(idx_neg), replace=True)
        idx = np.concatenate([rs_pos, rs_neg])
        yt = y_true[idx]
        sc = scores[idx]
        # Skip pathological resamples (should not occur with stratified sampling)
        if len(np.unique(yt)) < 2:
            tpr_boots[i, :] = np.nan
            auc_boots[i] = np.nan
            continue
        fpr_b, tpr_b, _ = roc_curve(yt, sc)
        tpr_boots[i, :] = np.interp(fpr_grid, fpr_b, tpr_b)
        auc_boots[i] = roc_auc_score(yt, sc)
    auc_ci = np.nanpercentile(auc_boots, [2.5, 97.5])
    tpr_lo = np.nanpercentile(tpr_boots, 2.5, axis=0)
    tpr_hi = np.nanpercentile(tpr_boots, 97.5, axis=0)
    return auc_point, auc_ci, fpr_grid, tpr_lo, tpr_hi


def render_figure1(df: pd.DataFrame, out: Path) -> None:
    y = df["binary_splice_affecting"].to_numpy().astype(int)
    fig, ax = plt.subplots(figsize=(9, 8))

    predictors = [
        ("AlphaGenome Rule A (composite)", "score_A", "#1f77b4", "-"),
        ("AlphaGenome Rule B (per-junction)", "score_B", "#ff7f0e", "--"),
        ("SpliceAI", "score_SAI", "#d62728", "-"),
        ("Pangolin", "score_PGL", "#2ca02c", "-"),
    ]

    for label, col, colour, linestyle in predictors:
        rng = np.random.default_rng(SEED)
        auc_pt, auc_ci, fpr_grid, tpr_lo, tpr_hi = stratified_bootstrap_roc(
            y, df[col].to_numpy(), rng
        )
        fpr, tpr, _ = roc_curve(y, df[col].to_numpy())
        legend_label = f"{label}  AUC {auc_pt:.3f}  [{auc_ci[0]:.3f}, {auc_ci[1]:.3f}]"
        ax.plot(fpr, tpr, color=colour, linewidth=2.2, linestyle=linestyle, label=legend_label)
        ax.fill_between(fpr_grid, tpr_lo, tpr_hi, color=colour, alpha=0.12)

    ax.plot([0, 1], [0, 1], "--", color="grey", linewidth=1.0, label="Chance")
    ax.set_xlabel("False positive rate", fontsize=12)
    ax.set_ylabel("True positive rate", fontsize=12)
    ax.set_title(
        "Stratum 2 (n=52) — ROC on binary splice-affecting label\n"
        "AlphaGenome scored under two pre-declared rules (Amendment #2)\n"
        "shaded bands: 10,000-iter stratified bootstrap 95% CI",
        fontsize=12,
    )
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.02)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower right", fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _monotone_fit(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rank-based monotone fit: sort by predictor score, running mean over
    a small window on the corresponding %-normally-spliced values.
    Purely for visual reference; not an inferential fit."""
    order = np.argsort(x)
    x_sorted = x[order]
    y_sorted = y[order]
    n = len(x_sorted)
    # window of 3 for n=8 gives a visually smooth monotone-ish line
    window = max(2, min(3, n - 1))
    y_smoothed = np.array(
        [
            np.mean(y_sorted[max(0, i - window // 2) : min(n, i + window // 2 + 1)])
            for i in range(n)
        ]
    )
    return x_sorted, y_smoothed


def render_figure2(df: pd.DataFrame, out: Path) -> None:
    rng = np.random.default_rng(SEED)
    y = df["pct_normally_spliced"].to_numpy()
    n = len(y)

    panels = [
        ("AlphaGenome Rule A (composite)", "score_A", "#1f77b4"),
        ("AlphaGenome Rule B (per-junction)", "score_B", "#ff7f0e"),
        ("SpliceAI", "score_SAI", "#d62728"),
        ("Pangolin", "score_PGL", "#2ca02c"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharey=True)
    axes = axes.ravel()

    for ax, (label, col, colour) in zip(axes, panels):
        x = df[col].to_numpy()
        rho, _ = spearmanr(x, y)
        # Deterministic bootstrap using an independent RNG per panel seeded from SEED
        panel_rng = np.random.default_rng(SEED + hash(col) % (2**32))
        rho_boots = np.empty(N_BOOT)
        for i in range(N_BOOT):
            idx = panel_rng.integers(0, n, size=n)
            if len(np.unique(x[idx])) < 2 or len(np.unique(y[idx])) < 2:
                rho_boots[i] = np.nan
                continue
            rho_boots[i], _ = spearmanr(x[idx], y[idx])
        ci = np.nanpercentile(rho_boots, [2.5, 97.5])
        ax.scatter(x, y, s=90, edgecolor="black", facecolor=colour, linewidths=0.6, zorder=3)
        fit_x, fit_y = _monotone_fit(x, y)
        ax.plot(fit_x, fit_y, color=colour, linewidth=1.7, alpha=0.75)
        ax.set_title(
            f"{label}\nSpearman ρ = {rho:+.3f}  [{ci[0]:+.2f}, {ci[1]:+.2f}]",
            fontsize=11,
        )
        ax.set_xlabel(f"{label.split(' (')[0]} score", fontsize=11)
        ax.set_ylim(-5, 110)
        ax.grid(True, alpha=0.25)
    axes[0].set_ylabel("% normally spliced transcript", fontsize=11)
    axes[2].set_ylabel("% normally spliced transcript", fontsize=11)

    fig.suptitle(
        "Stratum 1 (n=8) — predictor score vs quantitative residual splicing\n"
        "AlphaGenome scored under two pre-declared rules (Amendment #2)\n"
        "pre-specified expected direction: NEGATIVE Spearman ρ",
        fontsize=13,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    s2 = load_stratum2()
    print(f"Stratum 2: {len(s2)} variants, {int(s2['binary_splice_affecting'].sum())} positives")
    render_figure1(s2, FIG_DIR / "figure1_roc_stratum2.png")
    print(f"Wrote {FIG_DIR / 'figure1_roc_stratum2.png'}")

    s1 = load_stratum1()
    print(f"Stratum 1: {len(s1)} variants")
    render_figure2(s1, FIG_DIR / "figure2_stratum1_scatter.png")
    print(f"Wrote {FIG_DIR / 'figure2_stratum1_scatter.png'}")


if __name__ == "__main__":
    main()
