"""
Positive-control benchmark: raw nucleotide distance from each variant to the
nearest canonical CFTR splice site.

Any predictor whose AUC lower CI does not exceed the AUC upper CI of this
distance benchmark on Stratum 2 is flagged in the results, per the
pre-registered analysis plan (statistical_models.pdf §4, Positive controls).

Deterministic, no API calls, no external model dependency beyond a CFTR
transcript annotation. Serves two purposes:
  1. Proves the analysis pipeline end-to-end before the first live AlphaGenome
     query.
  2. Provides the floor benchmark that any useful predictor must beat.

Score convention: distance is a *negative* predictor of splice effect
(variants closer to a splice site are more likely to be splice-affecting), so
the score S_dist(v) is defined as the reciprocal of (1 + distance in nt) to
give a non-negative, non-inferior-larger score, matching the direction of the
DL/SpliceAI/Pangolin scores.

Usage:
    python -m src.controls.splice_site_distance_benchmark \\
        --variants data/variant_list_frozen.tsv \\
        --transcript data/reference/CFTR_NM_000492.4_exons.tsv \\
        --output results/controls/splice_site_distance_scores.tsv
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def load_exon_boundaries(transcript_tsv: Path) -> list[tuple[int, int]]:
    """Load CFTR exon start/end genomic coordinates (GRCh38, 1-based inclusive).

    Expected columns: exon_number, start, end (1-based inclusive, GRCh38 forward-strand).
    """
    df = pd.read_csv(transcript_tsv, sep="\t")
    required = {"exon_number", "start", "end"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"transcript TSV missing required columns: {missing}")
    return list(zip(df["start"].astype(int), df["end"].astype(int)))


def canonical_splice_sites(exons: list[tuple[int, int]]) -> list[int]:
    """Return the list of canonical splice-site genomic coordinates.

    For each intron there are two canonical splice sites: the donor at the 5'
    end of the intron (exon.end + 1) and the acceptor at the 3' end of the
    intron (next_exon.start - 1). First-exon 5' and last-exon 3' UTRs have no
    associated splice site.
    """
    sites: list[int] = []
    for i in range(len(exons) - 1):
        donor = exons[i][1] + 1        # first intronic nt after exon i
        acceptor = exons[i + 1][0] - 1  # last intronic nt before exon i+1
        sites.extend([donor, acceptor])
    return sites


def distance_to_nearest_splice_site(pos: int, sites: list[int]) -> int:
    """Minimum absolute nucleotide distance from `pos` to any site in `sites`."""
    return min(abs(pos - s) for s in sites)


def score_variants(
    variants: pd.DataFrame,
    exons: list[tuple[int, int]],
    pos_column: str = "grch38_pos",
) -> pd.DataFrame:
    """Compute the reciprocal-distance splice-site score for each variant.

    Args:
        variants: DataFrame with at least a `grch38_pos` column (int).
        exons: exon boundary list from `load_exon_boundaries`.
        pos_column: name of the genomic-position column.

    Returns:
        A copy of `variants` with two added columns:
          - `distance_to_splice_site` (int, nucleotides)
          - `score_splice_site_distance` (float, in (0, 1])
    """
    sites = canonical_splice_sites(exons)
    if not sites:
        raise ValueError("no splice sites derived from the exon annotation")

    df = variants.copy()
    df["distance_to_splice_site"] = df[pos_column].astype(int).apply(
        lambda p: distance_to_nearest_splice_site(p, sites)
    )
    df["score_splice_site_distance"] = 1.0 / (1.0 + df["distance_to_splice_site"])
    return df


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Positive-control splice-site distance benchmark for the "
                    "CFTR AlphaGenome benchmark."
    )
    p.add_argument("--variants", required=True, type=Path,
                   help="Frozen variant list TSV (must include grch38_pos column).")
    p.add_argument("--transcript", required=True, type=Path,
                   help="CFTR exon coordinates TSV (columns: exon_number, start, end).")
    p.add_argument("--output", required=True, type=Path,
                   help="Output TSV path.")
    p.add_argument("--pos-column", default="grch38_pos",
                   help="Genomic-position column name in the variants TSV.")
    p.add_argument("--log-level", default="INFO")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=args.log_level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    logger.info("Loading variants from %s", args.variants)
    variants = pd.read_csv(args.variants, sep="\t")
    logger.info("Loaded %d variants", len(variants))

    logger.info("Loading exon boundaries from %s", args.transcript)
    exons = load_exon_boundaries(args.transcript)
    logger.info("Loaded %d exons", len(exons))

    logger.info("Scoring variants against canonical splice sites")
    scored = score_variants(variants, exons, pos_column=args.pos_column)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    scored.to_csv(args.output, sep="\t", index=False)
    logger.info("Wrote %d rows to %s", len(scored), args.output)


if __name__ == "__main__":
    main()
