"""
Smoke tests for the positive-control splice-site-distance benchmark.

Uses the real Ensembl-derived CFTR exon annotation frozen at
data/reference/CFTR_ENST00000003084_11_exons_GRCh38.tsv but a small synthetic
variant list, so the end-to-end I/O and score direction are exercised without
touching the OSF-registered variant TSV.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from src.controls.splice_site_distance_benchmark import (
    canonical_splice_sites,
    distance_to_nearest_splice_site,
    load_exon_boundaries,
    score_variants,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
EXON_TSV = REPO_ROOT / "data" / "reference" / "CFTR_ENST00000003084_11_exons_GRCh38.tsv"


def test_exon_annotation_loads():
    exons = load_exon_boundaries(EXON_TSV)
    # ENST00000003084.11 has 27 exons on chr7 forward strand
    assert len(exons) == 27
    # First exon starts at 117,480,025 per Ensembl release 116
    assert exons[0] == (117480025, 117480147)
    # Last exon ends at 117,668,665
    assert exons[-1] == (117666908, 117668665)


def test_splice_sites_are_between_exons():
    exons = load_exon_boundaries(EXON_TSV)
    sites = canonical_splice_sites(exons)
    # 27 exons → 26 introns → 52 canonical splice sites
    assert len(sites) == 26 * 2

    # First donor is 1 nt after the first exon's end
    assert sites[0] == exons[0][1] + 1
    # First acceptor is 1 nt before the second exon's start
    assert sites[1] == exons[1][0] - 1


def test_variant_on_splice_site_has_max_score():
    exons = load_exon_boundaries(EXON_TSV)
    sites = canonical_splice_sites(exons)
    first_donor = sites[0]

    variants = pd.DataFrame({
        "variant_id": ["donor_hit", "far_intronic"],
        "grch38_pos": [first_donor, first_donor + 10_000],
    })
    scored = score_variants(variants, exons)

    donor_row = scored[scored["variant_id"] == "donor_hit"].iloc[0]
    far_row = scored[scored["variant_id"] == "far_intronic"].iloc[0]
    assert donor_row["distance_to_splice_site"] == 0
    assert donor_row["score_splice_site_distance"] == 1.0
    assert far_row["distance_to_splice_site"] > 0
    assert far_row["score_splice_site_distance"] < 1.0
    assert donor_row["score_splice_site_distance"] > far_row["score_splice_site_distance"]


def test_score_monotonic_in_distance():
    exons = load_exon_boundaries(EXON_TSV)
    sites = canonical_splice_sites(exons)
    donor = sites[0]

    variants = pd.DataFrame({
        "variant_id": ["a", "b", "c", "d"],
        "grch38_pos": [donor, donor + 1, donor + 100, donor + 1000],
    })
    scored = score_variants(variants, exons).sort_values("distance_to_splice_site")
    # Scores must strictly decrease as distance increases
    scores = scored["score_splice_site_distance"].tolist()
    assert all(scores[i] > scores[i + 1] for i in range(len(scores) - 1))


def test_missing_required_columns_raises(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("start\tend\n1\t10\n")
    with pytest.raises(ValueError, match="missing required columns"):
        load_exon_boundaries(bad)


def test_distance_symmetric_around_site():
    sites = [100]
    assert distance_to_nearest_splice_site(90, sites) == 10
    assert distance_to_nearest_splice_site(110, sites) == 10
    assert distance_to_nearest_splice_site(100, sites) == 0
