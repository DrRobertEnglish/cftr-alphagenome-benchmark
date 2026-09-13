"""Tests for the negative-control label-permutation harness."""
from __future__ import annotations

import numpy as np
import pandas as pd

from src.controls.label_permutation_null import (
    run_permutation_null,
    summarise,
)


def _make_variants(n_pos: int = 43, n_neg: int = 9) -> pd.DataFrame:
    ids = [f"v{i:03d}" for i in range(n_pos + n_neg)]
    return pd.DataFrame({
        "variant_id": ids,
        "stratum": ["2"] * (n_pos + n_neg),
        "binary_splice_affecting": [1] * n_pos + [0] * n_neg,
    })


def _make_scores(variants: pd.DataFrame, seed: int = 1) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "variant_id": variants["variant_id"],
        "score_splice_site_distance": rng.random(len(variants)),
    })


def test_null_mean_close_to_half():
    v = _make_variants()
    s = _make_scores(v)
    aucs = run_permutation_null(v, s, iters=500, seed=42)
    assert 0.45 < aucs["auc_null"].mean() < 0.55


def test_reproducible_seed():
    v = _make_variants()
    s = _make_scores(v, seed=7)
    a1 = run_permutation_null(v, s, iters=100, seed=99)
    a2 = run_permutation_null(v, s, iters=100, seed=99)
    assert np.allclose(a1["auc_null"], a2["auc_null"])


def test_different_seeds_produce_different_draws():
    v = _make_variants()
    s = _make_scores(v, seed=7)
    a1 = run_permutation_null(v, s, iters=100, seed=1)
    a2 = run_permutation_null(v, s, iters=100, seed=2)
    assert not np.allclose(a1["auc_null"], a2["auc_null"])


def test_summarise_returns_expected_keys():
    v = _make_variants()
    s = _make_scores(v)
    aucs = run_permutation_null(v, s, iters=200, seed=0)
    summary = summarise(aucs)
    for k in ["n_iterations", "mean_auc", "median_auc",
              "ci_low_2p5", "ci_high_97p5", "sd_auc"]:
        assert k in summary
    assert summary["n_iterations"] == 200
    assert summary["ci_low_2p5"] < summary["median_auc"] < summary["ci_high_97p5"]
