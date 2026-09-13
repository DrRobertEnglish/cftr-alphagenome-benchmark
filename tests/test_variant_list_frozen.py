"""Integrity tests for the frozen variant list.

These are frozen-input checks: they lock in the specific counts, class
balance, and coordinate range that the OSF-registered protocol depends on.
Any change here is a scientific claim, not a code refactor, and must be
accompanied by an entry in docs/protocol/AMENDMENTS.md.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
TSV = REPO_ROOT / "data" / "variant_list_frozen.tsv"

# CFTR gene span on GRCh38 chr7, from Ensembl release 116 lookup
CFTR_START = 117_480_025
CFTR_END = 117_668_665


def test_frozen_tsv_exists():
    assert TSV.is_file()


def test_row_count_matches_registered_harvest():
    df = pd.read_csv(TSV, sep="\t")
    # 52 Joynt Stratum 2 + 8 Stratum 1 quantitative (the Masvidal complex
    # allele is excluded — documented in build_variant_list.py)
    assert len(df) == 60


def test_stratum_2_counts_match_registered_harvest():
    df = pd.read_csv(TSV, sep="\t")
    s2 = df[df["stratum"].astype(str) == "2"]
    assert len(s2) == 52
    n_splice_affecting = int(s2["binary_splice_affecting"].sum())
    n_no_effect = int((s2["binary_splice_affecting"] == 0).sum())
    assert n_splice_affecting == 43, "Joynt Stratum 2 must have 43 splice-affecting"
    assert n_no_effect == 9, "Joynt Stratum 2 must have 9 no-splice-effect"


def test_stratum_1_row_count():
    df = pd.read_csv(TSV, sep="\t")
    s1 = df[df["stratum"].astype(str) == "1"]
    assert len(s1) == 8, "Stratum 1 excludes the Masvidal complex allele"


def test_all_positions_inside_cftr_gene_span():
    df = pd.read_csv(TSV, sep="\t")
    pos = df["grch38_pos"].astype(int)
    assert (pos >= CFTR_START).all()
    assert (pos <= CFTR_END).all()


def test_all_variants_on_chr7():
    df = pd.read_csv(TSV, sep="\t")
    assert (df["chrom"].astype(str) == "7").all()


def test_stratum_1_has_quantitative_measurements():
    df = pd.read_csv(TSV, sep="\t")
    s1 = df[df["stratum"].astype(str) == "1"]
    assert s1["pct_normally_spliced"].notna().all()
    pct = s1["pct_normally_spliced"].astype(float)
    assert (pct >= 0).all() and (pct <= 100).all()


def test_no_duplicate_variant_ids():
    df = pd.read_csv(TSV, sep="\t")
    assert df["variant_id"].is_unique
