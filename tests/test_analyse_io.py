"""Tests for the shared IO helpers used by every analysis runner."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from src.analyse.io import (
    build_stratum1_frame,
    build_stratum2_frame,
    load_predictor_scores,
    load_variants,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
FROZEN_TSV = REPO_ROOT / "data" / "variant_list_frozen.tsv"


def test_load_variants_schema():
    df = load_variants(FROZEN_TSV)
    assert len(df) == 60
    assert (df["stratum"].isin({"1", "2"})).all()
    assert df["binary_splice_affecting"].dtype.kind in ("i", "u")


def test_load_variants_missing_column_raises(tmp_path):
    bad = tmp_path / "bad.tsv"
    pd.DataFrame({"variant_id": ["v1"]}).to_csv(bad, sep="\t", index=False)
    with pytest.raises(ValueError, match="missing required columns"):
        load_variants(bad)


def test_load_predictor_scores_renames_to_score(tmp_path):
    p = tmp_path / "preds.tsv"
    pd.DataFrame({
        "variant_id": ["v1", "v2"],
        "my_score": [0.1, 0.9],
    }).to_csv(p, sep="\t", index=False)
    out = load_predictor_scores(p, "my_score")
    assert list(out.columns) == ["variant_id", "score"]
    assert out["score"].tolist() == [0.1, 0.9]


def test_load_predictor_scores_missing_score_col_raises(tmp_path):
    p = tmp_path / "preds.tsv"
    pd.DataFrame({"variant_id": ["v1"], "wrong_col": [0.5]}).to_csv(p, sep="\t", index=False)
    with pytest.raises(ValueError, match="missing required score column"):
        load_predictor_scores(p, "score")


def test_load_predictor_scores_rejects_nan(tmp_path):
    p = tmp_path / "preds.tsv"
    pd.DataFrame({
        "variant_id": ["v1", "v2"],
        "score": [0.5, None],
    }).to_csv(p, sep="\t", index=False)
    with pytest.raises(ValueError, match="missing score values"):
        load_predictor_scores(p, "score")


def test_build_stratum2_frame_uses_left_join(tmp_path):
    v = load_variants(FROZEN_TSV)
    # Two predictors each score only a subset of variants
    ids_a = v[v["stratum"] == "2"]["variant_id"].iloc[:40].tolist()
    ids_b = v[v["stratum"] == "2"]["variant_id"].iloc[10:].tolist()
    a = pd.DataFrame({"variant_id": ids_a, "score": 0.5})
    b = pd.DataFrame({"variant_id": ids_b, "score": 0.7})
    frame = build_stratum2_frame(v, a, b)
    assert len(frame) == 52  # all Stratum 2 rows kept
    assert frame["score_a"].isna().sum() == 12  # 52 - 40 missing from A
    assert frame["score_b"].isna().sum() == 10  # first 10 missing from B


def test_build_stratum1_frame_keeps_only_stratum1():
    v = load_variants(FROZEN_TSV)
    s1_ids = v[v["stratum"] == "1"]["variant_id"].tolist()
    scores = pd.DataFrame({"variant_id": s1_ids, "score": 0.4})
    frame = build_stratum1_frame(v, scores)
    assert len(frame) == 8
    assert (frame["stratum"] == "1").all()
