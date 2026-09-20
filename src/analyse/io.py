"""Shared IO helpers for reading the frozen variant list and predictor score
files, and merging them into an analysis frame.

Every analysis script in src/analyse/ reads from these same functions so the
pre-registered inclusion/exclusion rules are enforced in one place.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

VARIANT_TSV_COLUMNS_REQUIRED = {
    "variant_id",
    "stratum",
    "binary_splice_affecting",
    "pct_normally_spliced",
    "chrom",
    "grch38_pos",
}


def load_variants(variants_tsv: Path) -> pd.DataFrame:
    """Load the frozen variant list and sanity-check its schema."""
    df = pd.read_csv(variants_tsv, sep="\t")
    missing = VARIANT_TSV_COLUMNS_REQUIRED - set(df.columns)
    if missing:
        raise ValueError(f"variant TSV missing required columns: {sorted(missing)}")
    df["stratum"] = df["stratum"].astype(str)
    df["binary_splice_affecting"] = df["binary_splice_affecting"].astype(int)
    return df


def load_predictor_scores(scores_tsv: Path, score_column: str) -> pd.DataFrame:
    """Load a predictor's per-variant scores.

    Args:
        scores_tsv: TSV with at least (variant_id, <score_column>) columns.
        score_column: which column holds the predictor's continuous score.

    Returns:
        DataFrame with exactly columns (variant_id, score).
    """
    df = pd.read_csv(scores_tsv, sep="\t")
    if "variant_id" not in df.columns:
        raise ValueError(f"{scores_tsv}: missing required column 'variant_id'")
    if score_column not in df.columns:
        raise ValueError(
            f"{scores_tsv}: missing required score column '{score_column}'; "
            f"available: {list(df.columns)}"
        )
    out = df[["variant_id", score_column]].rename(columns={score_column: "score"})
    if out["score"].isna().any():
        n = int(out["score"].isna().sum())
        raise ValueError(f"{scores_tsv}: {n} rows have missing score values")
    return out


def build_stratum2_frame(
    variants: pd.DataFrame,
    scores_a: pd.DataFrame,
    scores_b: pd.DataFrame,
) -> pd.DataFrame:
    """Build the frame the paired-DeLong test consumes for H1 / H2.

    Pre-registered inclusion:
      - Stratum 2 only (n=52 Joynt binary).
      - Both predictors must have scored every included variant.
    Any variant missing from either predictor file drops out and is reported.
    """
    s2 = variants[variants["stratum"] == "2"].copy()
    merged = (
        s2.merge(scores_a.rename(columns={"score": "score_a"}), on="variant_id", how="left")
          .merge(scores_b.rename(columns={"score": "score_b"}), on="variant_id", how="left")
    )
    return merged


def build_stratum1_frame(
    variants: pd.DataFrame,
    scores: pd.DataFrame,
) -> pd.DataFrame:
    """Build the frame the Spearman-rho analysis consumes for H3."""
    s1 = variants[variants["stratum"] == "1"].copy()
    merged = s1.merge(scores.rename(columns={"score": "score"}), on="variant_id", how="left")
    return merged
