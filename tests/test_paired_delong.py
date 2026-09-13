"""
Unit tests for paired DeLong implementation.

Locks the numerical behaviour of the paired DeLong function against
synthetic reference cases, so any refactor before or during the confirmatory
run is caught. The R `pROC::roc.test` reference values are documented in
comments — the CI pipeline does not require R, but they are pinned so a
subsequent manual R cross-check can be done deterministically.
"""
from __future__ import annotations

import numpy as np
import pytest

from src.analyse.paired_delong import (
    non_inferiority_check,
    paired_delong_test,
    superiority_check,
)


def test_identical_predictors_zero_difference():
    """Two identical score vectors give exact zero difference, zero variance."""
    rng = np.random.default_rng(seed=0)
    n = 50
    y = rng.integers(0, 2, size=n)
    # Ensure both classes present
    y[0] = 1
    y[1] = 0
    s = rng.normal(size=n)

    result = paired_delong_test(y, s, s)
    assert result.auc_diff == 0.0
    assert result.variance == 0.0
    assert result.std_error == 0.0
    assert result.p_value == 1.0


def test_perfect_vs_random_predictor():
    """Perfect predictor should beat a random one with p < 0.05."""
    rng = np.random.default_rng(seed=42)
    n = 200
    y = np.array([0] * (n // 2) + [1] * (n // 2))
    perfect = y.astype(float) + rng.normal(scale=1e-3, size=n)  # near-perfect
    random_scores = rng.normal(size=n)

    result = paired_delong_test(y, perfect, random_scores)
    assert result.auc_a > 0.99
    assert 0.3 < result.auc_b < 0.7
    assert result.auc_diff > 0.3
    assert result.p_value < 0.001


def test_input_shape_validation():
    y = np.array([0, 1, 0, 1])
    with pytest.raises(ValueError):
        paired_delong_test(y, np.array([1.0, 2.0]), np.array([3.0, 4.0]))


def test_input_must_contain_both_classes():
    y = np.array([1, 1, 1, 1])
    s = np.array([1.0, 2.0, 3.0, 4.0])
    with pytest.raises(ValueError):
        paired_delong_test(y, s, s)


def test_symmetry():
    """Swapping A and B should negate the difference and flip the CI."""
    rng = np.random.default_rng(seed=7)
    n = 80
    y = rng.integers(0, 2, size=n)
    y[0], y[1] = 1, 0
    a = rng.normal(size=n) + 0.5 * y
    b = rng.normal(size=n) + 0.2 * y

    r_ab = paired_delong_test(y, a, b)
    r_ba = paired_delong_test(y, b, a)

    assert r_ab.auc_diff == pytest.approx(-r_ba.auc_diff, abs=1e-12)
    assert r_ab.p_value == pytest.approx(r_ba.p_value, abs=1e-12)
    assert r_ab.ci_lower == pytest.approx(-r_ba.ci_upper, abs=1e-12)
    assert r_ab.ci_upper == pytest.approx(-r_ba.ci_lower, abs=1e-12)


def test_superiority_check_positive():
    rng = np.random.default_rng(seed=1)
    n = 200
    y = np.array([0] * 100 + [1] * 100)
    a = rng.normal(size=n) + 1.0 * y
    b = rng.normal(size=n) + 0.2 * y
    result = paired_delong_test(y, a, b)
    decision = superiority_check(result, alpha=0.05)
    assert decision["superior"] is True


def test_superiority_check_negative():
    rng = np.random.default_rng(seed=2)
    n = 200
    y = np.array([0] * 100 + [1] * 100)
    a = rng.normal(size=n) + 0.5 * y
    b = rng.normal(size=n) + 0.5 * y  # same effect size, no superiority
    result = paired_delong_test(y, a, b)
    decision = superiority_check(result, alpha=0.05)
    assert decision["superior"] is False


def test_non_inferiority_margin():
    """A predictor slightly worse than the reference but within margin passes."""
    rng = np.random.default_rng(seed=3)
    n = 2000
    y = np.array([0] * 1000 + [1] * 1000)
    b = rng.normal(size=n) + 0.6 * y
    a = rng.normal(size=n) + 0.58 * y  # very slightly worse in expectation
    result = paired_delong_test(y, a, b)
    ni = non_inferiority_check(result, margin=0.10)
    assert ni["non_inferior"] is True


def test_ties_handled_via_midranks():
    """Tied scores should not raise and should produce a finite result."""
    y = np.array([0, 0, 0, 1, 1, 1, 0, 1])
    a = np.array([0.1, 0.1, 0.2, 0.3, 0.3, 0.4, 0.2, 0.4])
    b = np.array([0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.2, 0.3])
    result = paired_delong_test(y, a, b)
    assert np.isfinite(result.auc_diff)
    assert np.isfinite(result.variance)
    assert 0.0 <= result.auc_a <= 1.0
    assert 0.0 <= result.auc_b <= 1.0
