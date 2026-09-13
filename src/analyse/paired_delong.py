"""
Paired DeLong test on ROC-AUC differences, following Sun & Xu (2014):
"Fast implementation of DeLong's algorithm for comparing the areas under
correlated receiver operating characteristic curves."
IEEE Signal Processing Letters 21(11):1389-1393.

Serves as the Python cross-check implementation for the primary confirmatory
analysis. The R `pROC::roc.test(method="delong", paired=TRUE)` is the primary
implementation. Any disagreement > 0.005 in an AUC point estimate or > 0.01
in the two-sided p-value triggers a disclosed manual investigation, per the
pre-registered analysis plan.

Deterministic. No random component. Fixed test suite in tests/test_paired_delong.py
locks the expected outputs against synthetic and reference datasets.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats


@dataclass(frozen=True)
class DeLongResult:
    auc_a: float
    auc_b: float
    auc_diff: float          # auc_a - auc_b
    variance: float          # variance of the paired difference
    std_error: float
    ci_lower: float          # two-sided 95% CI on auc_diff
    ci_upper: float
    z: float                 # test statistic under H0: auc_diff = 0
    p_value: float           # two-sided p-value
    n_pos: int
    n_neg: int

    def as_dict(self) -> dict:
        return {
            "auc_a": self.auc_a,
            "auc_b": self.auc_b,
            "auc_diff": self.auc_diff,
            "variance": self.variance,
            "std_error": self.std_error,
            "ci_lower": self.ci_lower,
            "ci_upper": self.ci_upper,
            "z": self.z,
            "p_value": self.p_value,
            "n_pos": self.n_pos,
            "n_neg": self.n_neg,
        }


def _compute_midrank(x: np.ndarray) -> np.ndarray:
    """Compute midranks of a 1D array (tie-averaged ranks), 1-based."""
    order = np.argsort(x, kind="mergesort")
    x_sorted = x[order]
    n = len(x)
    ranks = np.empty(n, dtype=np.float64)

    i = 0
    while i < n:
        j = i
        while j < n and x_sorted[j] == x_sorted[i]:
            j += 1
        # midrank for ties spanning [i, j)
        mid = 0.5 * (i + j - 1) + 1.0  # +1 for 1-based ranks
        ranks[i:j] = mid
        i = j

    out = np.empty(n, dtype=np.float64)
    out[order] = ranks
    return out


def _fast_delong(
    predictions_matrix: np.ndarray,
    label_1_count: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Fast DeLong AUC + covariance computation for k predictors.

    Sun & Xu 2014, Algorithm 2.

    Args:
        predictions_matrix: shape (k, n). Row m holds predictor m's scores;
            the first `label_1_count` columns are positive cases, the rest
            are negative cases.
        label_1_count: number of positive cases (label = 1).

    Returns:
        aucs: shape (k,)
        cov: shape (k, k)
    """
    k = predictions_matrix.shape[0]
    m = label_1_count                          # positives
    n = predictions_matrix.shape[1] - m        # negatives

    positive = predictions_matrix[:, :m]
    negative = predictions_matrix[:, m:]

    tx = np.empty((k, m), dtype=np.float64)
    ty = np.empty((k, n), dtype=np.float64)
    tz = np.empty((k, m + n), dtype=np.float64)

    for r in range(k):
        tx[r] = _compute_midrank(positive[r])
        ty[r] = _compute_midrank(negative[r])
        tz[r] = _compute_midrank(predictions_matrix[r])

    aucs = (tz[:, :m].sum(axis=1) / m - (m + 1) / 2.0) / n

    v01 = (tz[:, :m] - tx) / n           # shape (k, m)
    v10 = 1.0 - (tz[:, m:] - ty) / m     # shape (k, n)

    sx = np.cov(v01, ddof=1) if k > 1 else np.array([[v01.var(ddof=1)]])
    sy = np.cov(v10, ddof=1) if k > 1 else np.array([[v10.var(ddof=1)]])

    # np.cov returns scalar (as 0-d array) when only one row is passed; ensure 2D
    sx = np.atleast_2d(sx)
    sy = np.atleast_2d(sy)

    cov = sx / m + sy / n
    return aucs, cov


def paired_delong_test(
    y_true: np.ndarray,
    scores_a: np.ndarray,
    scores_b: np.ndarray,
    alpha: float = 0.05,
) -> DeLongResult:
    """Paired DeLong test comparing two AUCs on the same labels.

    Two-sided test. Handles ties via midranks.

    Args:
        y_true: shape (n,), binary labels {0, 1}.
        scores_a, scores_b: shape (n,), continuous scores for the two predictors.
        alpha: two-sided significance level for the CI. Default 0.05.

    Returns:
        DeLongResult with paired difference auc_a - auc_b, its variance, CI,
        z statistic, and two-sided p-value.
    """
    y_true = np.asarray(y_true).astype(np.int64)
    scores_a = np.asarray(scores_a, dtype=np.float64)
    scores_b = np.asarray(scores_b, dtype=np.float64)

    if not (y_true.shape == scores_a.shape == scores_b.shape):
        raise ValueError("y_true, scores_a, scores_b must all have the same shape")
    if y_true.ndim != 1:
        raise ValueError("inputs must be 1-D arrays")
    unique_labels = np.unique(y_true)
    if not np.array_equal(unique_labels, np.array([0, 1])):
        raise ValueError(f"y_true must contain both 0 and 1; got {unique_labels}")

    # Reorder so positives come first
    order = np.argsort(-y_true, kind="mergesort")
    y = y_true[order]
    a = scores_a[order]
    b = scores_b[order]
    n_pos = int(y.sum())
    n_neg = int(len(y) - n_pos)

    preds = np.stack([a, b], axis=0)
    aucs, cov = _fast_delong(preds, n_pos)

    auc_a, auc_b = float(aucs[0]), float(aucs[1])
    diff = auc_a - auc_b
    var = float(cov[0, 0] + cov[1, 1] - 2.0 * cov[0, 1])
    if var < 0:
        var = 0.0
    se = float(np.sqrt(var))

    if se > 0:
        z = diff / se
        p_two_sided = float(2.0 * (1.0 - stats.norm.cdf(abs(z))))
    else:
        z = float("nan") if diff != 0 else 0.0
        p_two_sided = float("nan") if diff != 0 else 1.0

    z_crit = float(stats.norm.ppf(1.0 - alpha / 2.0))
    ci_lower = diff - z_crit * se
    ci_upper = diff + z_crit * se

    return DeLongResult(
        auc_a=auc_a,
        auc_b=auc_b,
        auc_diff=diff,
        variance=var,
        std_error=se,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        z=z,
        p_value=p_two_sided,
        n_pos=n_pos,
        n_neg=n_neg,
    )


def non_inferiority_check(
    result: DeLongResult, margin: float = 0.05,
) -> dict:
    """Apply the pre-registered H2 non-inferiority decision rule.

    Passes if the lower bound of the CI on (auc_a - auc_b) exceeds -margin.
    Margin is on the AUC scale.
    """
    lower = result.ci_lower
    passed = lower > -abs(margin)
    return {
        "margin": abs(margin),
        "ci_lower": lower,
        "ci_upper": result.ci_upper,
        "non_inferior": bool(passed),
        "decision": "non-inferior" if passed else "not non-inferior",
    }


def superiority_check(
    result: DeLongResult, alpha: float = 0.05,
) -> dict:
    """Apply the pre-registered H1 superiority decision rule.

    Passes if auc_a point estimate > auc_b AND two-sided p < alpha.
    """
    passed = (result.auc_diff > 0.0) and (result.p_value < alpha)
    return {
        "alpha": alpha,
        "auc_diff": result.auc_diff,
        "p_value": result.p_value,
        "superior": bool(passed),
        "decision": "superior" if passed else "not superior",
    }
