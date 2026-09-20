"""Tests for the predictor parsers that don't require the model weights.

Full-model smoke tests live outside pytest (they require the .venv-predictors
environment with tensorflow / pytorch installed and take minutes to run).
Here we exercise the pure-python parsing and schema logic against captured
fixtures so future edits don't silently break the output contract that
run_primary_analysis.py depends on.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.predict.run_pangolin import extract_max_abs_score, parse_pangolin_vcf

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_pangolin_extract_single_gene_gain_and_loss() -> None:
    """Standard single-gene block: max abs across gain/loss tokens wins."""
    info = "ENSG00000001626.18|9:0.03|-1:-0.84|Warnings:"
    score, genes = extract_max_abs_score(info)
    assert score == 0.84
    assert genes == "ENSG00000001626.18"


def test_pangolin_extract_multi_gene() -> None:
    """Multi-gene block: max abs across ALL genes and ALL tokens wins."""
    info = "ENSG00000001626.18|9:0.03|-1:-0.10|Warnings:,ENSG00000000971.15|9:0.55|-1:-0.05|Warnings:"
    score, genes = extract_max_abs_score(info)
    assert score == 0.55
    assert genes == "ENSG00000001626.18,ENSG00000000971.15"


def test_pangolin_extract_missing_info() -> None:
    """Absent Pangolin annotation returns zero-score."""
    for info in ("", "."):
        score, genes = extract_max_abs_score(info)
        assert score == 0.0
        assert genes == ""


def test_pangolin_parse_synthetic_vcf(tmp_path: Path) -> None:
    """End-to-end VCF parser test with a synthetic 2-variant fixture."""
    fixture = tmp_path / "pangolin_raw.vcf"
    fixture.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=Pangolin,Number=.,Type=String,Description="Pangolin scores">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "7\t117530975\tJoynt_I01\tA\tG\t.\t.\tPangolin=ENSG00000001626.18|9:0.03|-1:-0.84|Warnings:\n"
        "7\t117531012\tJoynt_I20\tG\tA\t.\t.\tPangolin=ENSG00000001626.18|-2:0.0|-3:-0.07|Warnings:\n"
    )
    rows = parse_pangolin_vcf(fixture)
    assert len(rows) == 2
    assert rows[0]["variant_id"] == "Joynt_I01"
    assert rows[0]["score"] == 0.84
    assert rows[1]["variant_id"] == "Joynt_I20"
    assert rows[1]["score"] == 0.07
    for row in rows:
        assert row["symbol"] == "CFTR"


def test_predictor_output_schema_contract() -> None:
    """SpliceAI and Pangolin outputs must share the columns run_primary_analysis needs.

    Skips gracefully if the full predictor outputs don't exist yet (CI doesn't
    have the model weights). If they do exist, verify the contract.
    """
    spliceai_path = REPO_ROOT / "results" / "predictions" / "spliceai.tsv"
    pangolin_path = REPO_ROOT / "results" / "predictions" / "pangolin.tsv"

    for path in (spliceai_path, pangolin_path):
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        # The columns run_primary_analysis / run_stratum1 depend on
        assert "variant_id" in df.columns, f"{path} missing variant_id"
        assert "score" in df.columns, f"{path} missing score"
        # Score must be numeric in [0, 1]
        assert df["score"].dtype.kind in "fi", f"{path} score not numeric"
        assert df["score"].min() >= 0.0, f"{path} has negative score"
        assert df["score"].max() <= 1.0, f"{path} has score > 1"
        # Should cover all 60 frozen variants
        assert len(df) == 60, f"{path} has {len(df)} rows, expected 60"
        # Variant IDs must be unique
        assert df["variant_id"].is_unique, f"{path} has duplicate variant_ids"
