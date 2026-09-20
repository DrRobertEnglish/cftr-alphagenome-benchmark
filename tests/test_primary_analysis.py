"""End-to-end tests for the primary analysis runner (H1 + H2)."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.analyse.run_primary_analysis import stratified_bootstrap_auc_ci

REPO_ROOT = Path(__file__).resolve().parents[1]
FROZEN_TSV = REPO_ROOT / "data" / "variant_list_frozen.tsv"


def _write_predictor_scores(path: Path, variant_ids: list[str], scores: list[float]) -> None:
    pd.DataFrame({"variant_id": variant_ids, "score": scores}).to_csv(
        path, sep="\t", index=False,
    )


def _stratum2_ids() -> list[str]:
    df = pd.read_csv(FROZEN_TSV, sep="\t")
    return df[df["stratum"].astype(str) == "2"]["variant_id"].tolist()


def _stratum2_labels() -> np.ndarray:
    df = pd.read_csv(FROZEN_TSV, sep="\t")
    return df[df["stratum"].astype(str) == "2"]["binary_splice_affecting"].to_numpy(int)


def test_stratified_bootstrap_ci_covers_true_auc():
    rng = np.random.default_rng(0)
    n = 52
    y = np.array([1] * 43 + [0] * 9)
    # Moderately informative scores: label + noise wide enough to give
    # non-degenerate bootstrap CIs on a 43:9 class-imbalanced set.
    scores = y + rng.normal(0, 0.8, size=n)
    boot = stratified_bootstrap_auc_ci(y, scores, iters=500, seed=1)
    assert 0.6 < boot.auc_median < 1.0
    # CI must be a proper interval, not a degenerate point
    assert boot.ci_low_2p5 < boot.ci_high_97p5
    assert boot.ci_low_2p5 <= boot.auc_median <= boot.ci_high_97p5


def test_pipeline_end_to_end_identical_predictors_gives_zero_diff(tmp_path):
    """If both predictors give the same scores, paired difference must be
    exactly 0 and H1 must return NOT SUPERIOR."""
    ids = _stratum2_ids()
    rng = np.random.default_rng(42)
    scores = rng.random(len(ids))
    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    _write_predictor_scores(a, ids, scores.tolist())
    _write_predictor_scores(b, ids, scores.tolist())
    out = tmp_path / "primary.tsv"

    cmd = [
        sys.executable, "-m", "src.analyse.run_primary_analysis",
        "--variants", str(FROZEN_TSV),
        "--scores-a", str(a), "--label-a", "A",
        "--scores-b", str(b), "--label-b", "B",
        "--score-column", "score",
        "--boot-iters", "50",  # keep the test fast
        "--output", str(out),
    ]
    subprocess.run(cmd, check=True, cwd=REPO_ROOT, env={"PYTHONPATH": str(REPO_ROOT), "PATH": ""})
    assert out.exists()
    payload = json.loads(out.with_suffix(".json").read_text())
    assert payload["paired_delong"]["auc_diff"] == pytest.approx(0.0, abs=1e-12)
    assert payload["h1_superiority"]["superior"] is False
    assert payload["h2_non_inferiority"]["non_inferior"] is True
    # Coverage: all 52 Stratum 2 variants scored by both
    assert payload["coverage"]["n_stratum2"] == 52
    assert payload["coverage"]["n_paired"] == 52
    assert payload["coverage"]["n_missing_a"] == 0
    assert payload["coverage"]["n_missing_b"] == 0


def test_pipeline_end_to_end_perfect_vs_random_gives_superior(tmp_path):
    """Perfectly informative predictor A vs random predictor B on real
    Stratum 2 labels must give H1 superior with p < 0.05."""
    ids = _stratum2_ids()
    y = _stratum2_labels()
    rng = np.random.default_rng(0)
    a_scores = y.astype(float) + rng.normal(0, 0.05, size=len(y))  # near-perfect
    b_scores = rng.random(len(y))                                   # random

    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    _write_predictor_scores(a, ids, a_scores.tolist())
    _write_predictor_scores(b, ids, b_scores.tolist())
    out = tmp_path / "primary.tsv"

    cmd = [
        sys.executable, "-m", "src.analyse.run_primary_analysis",
        "--variants", str(FROZEN_TSV),
        "--scores-a", str(a), "--label-a", "Perfect",
        "--scores-b", str(b), "--label-b", "Random",
        "--boot-iters", "50",
        "--output", str(out),
    ]
    subprocess.run(cmd, check=True, cwd=REPO_ROOT, env={"PYTHONPATH": str(REPO_ROOT), "PATH": ""})
    payload = json.loads(out.with_suffix(".json").read_text())
    assert payload["paired_delong"]["auc_a"] > 0.95
    assert payload["h1_superiority"]["superior"] is True
    assert payload["paired_delong"]["p_value"] < 0.05


def test_pipeline_refuses_when_predictor_missing_most_variants(tmp_path):
    """If one predictor covers <10 Stratum 2 variants, the runner must
    refuse rather than emit a false-confidence result."""
    ids = _stratum2_ids()
    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    rng = np.random.default_rng(0)
    _write_predictor_scores(a, ids[:5], rng.random(5).tolist())  # only 5 scored
    _write_predictor_scores(b, ids, rng.random(len(ids)).tolist())
    out = tmp_path / "primary.tsv"

    cmd = [
        sys.executable, "-m", "src.analyse.run_primary_analysis",
        "--variants", str(FROZEN_TSV),
        "--scores-a", str(a), "--label-a", "Partial",
        "--scores-b", str(b), "--label-b", "Full",
        "--boot-iters", "50",
        "--output", str(out),
    ]
    proc = subprocess.run(
        cmd, check=False, cwd=REPO_ROOT,
        env={"PYTHONPATH": str(REPO_ROOT), "PATH": ""},
        capture_output=True, text=True,
    )
    assert proc.returncode != 0
    assert "Refusing to run" in (proc.stdout + proc.stderr)
