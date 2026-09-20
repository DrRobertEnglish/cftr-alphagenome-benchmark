"""Tests for the Stratum 1 (H3) rank-correlation runner."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.analyse.run_stratum1 import bootstrap_spearman_ci

REPO_ROOT = Path(__file__).resolve().parents[1]
FROZEN_TSV = REPO_ROOT / "data" / "variant_list_frozen.tsv"


def _stratum1_ids() -> list[str]:
    df = pd.read_csv(FROZEN_TSV, sep="\t")
    return df[df["stratum"].astype(str) == "1"]["variant_id"].tolist()


def test_bootstrap_spearman_ci_perfect_positive():
    x = np.arange(20, dtype=float)
    y = 2 * x + 1.0
    boot = bootstrap_spearman_ci(x, y, iters=500, seed=1)
    assert boot["rho_median"] == pytest.approx(1.0, abs=1e-9)
    assert boot["ci_low_2p5"] == pytest.approx(1.0, abs=1e-9)


def test_bootstrap_spearman_ci_negative_association():
    rng = np.random.default_rng(0)
    x = np.arange(50, dtype=float)
    y = -x + rng.normal(0, 1.0, size=50)
    boot = bootstrap_spearman_ci(x, y, iters=500, seed=1)
    assert boot["rho_median"] < -0.8
    assert boot["ci_high_97p5"] < 0


def test_bootstrap_rejects_tiny_input():
    with pytest.raises(ValueError, match="at least 4"):
        bootstrap_spearman_ci(np.array([1.0, 2.0]), np.array([3.0, 4.0]))


def test_pipeline_end_to_end_runs_against_stratum1(tmp_path):
    ids = _stratum1_ids()
    scores = tmp_path / "scores.tsv"
    # Use a monotone-in-percent predictor: higher score = less normally spliced
    df = pd.read_csv(FROZEN_TSV, sep="\t")
    s1 = df[df["stratum"].astype(str) == "1"]
    fake_scores = 1.0 - s1["pct_normally_spliced"].astype(float) / 100.0
    pd.DataFrame({
        "variant_id": s1["variant_id"].tolist(),
        "score": fake_scores.tolist(),
    }).to_csv(scores, sep="\t", index=False)

    out = tmp_path / "s1.tsv"
    cmd = [
        sys.executable, "-m", "src.analyse.run_stratum1",
        "--variants", str(FROZEN_TSV),
        "--scores", str(scores),
        "--label", "Monotone",
        "--boot-iters", "200",
        "--output", str(out),
    ]
    subprocess.run(cmd, check=True, cwd=REPO_ROOT,
                   env={"PYTHONPATH": str(REPO_ROOT), "PATH": ""})
    payload = json.loads(out.with_suffix(".json").read_text())
    # Higher score = less normally spliced => Spearman ρ(score, %normal) must be strongly negative
    assert payload["spearman_rho"] < -0.9
    assert payload["n_variants"] == 8
    assert len(payload["variant_ids"]) == 8
    _ = ids  # silence lint
