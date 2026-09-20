#!/usr/bin/env bash
# End-to-end smoke test of the analysis pipeline.
#
# Uses the splice-site-distance floor benchmark as a stand-in "Predictor A",
# and a deterministically jittered copy of it as "Predictor B", to exercise
# the full paired-DeLong + non-inferiority + bootstrap-CI + H3 pipeline on
# real Stratum 2 and Stratum 1 rows.
#
# This is a mechanics smoke test only, NOT a scientific comparison. The two
# predictors here are trivially correlated by construction. The purpose is to
# prove the runners work end-to-end and that the output files are in the
# expected shape before real predictor scores (SpliceAI, Pangolin,
# AlphaGenome) are dropped in.
#
# Idempotent: overwrites results/smoke/ on every run.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT"

OUT_DIR="results/smoke"
mkdir -p "$OUT_DIR"

# Build a "Predictor A" TSV from the floor-benchmark scores in the expected schema
python - <<'PY'
import pandas as pd, numpy as np
df = pd.read_csv("results/controls/splice_site_distance_scores.tsv", sep="\t")
# Stratum 2 rows have binary_splice_affecting; drop Stratum 1 rows here
# (their scores get written out separately for the H3 runner).
out_a = df[["variant_id", "score_splice_site_distance"]].rename(
    columns={"score_splice_site_distance": "score"}
)
out_a.to_csv("results/smoke/predictor_a_floor.tsv", sep="\t", index=False)

# Predictor B: floor benchmark + tiny deterministic jitter
rng = np.random.default_rng(20260913)
jitter = rng.normal(0, 0.05, size=len(out_a))
out_b = out_a.copy()
out_b["score"] = out_a["score"] + jitter
out_b.to_csv("results/smoke/predictor_b_jittered.tsv", sep="\t", index=False)

# Stratum 1 predictor scores: reuse floor scores for the Stratum 1 subset
out_a.to_csv("results/smoke/predictor_a_stratum1.tsv", sep="\t", index=False)
print("Wrote predictor TSVs to results/smoke/")
PY

echo "--- Primary analysis (H1 + H2) ---"
python -m src.analyse.run_primary_analysis \
  --variants data/variant_list_frozen.tsv \
  --scores-a results/smoke/predictor_a_floor.tsv --label-a FloorBenchmark \
  --scores-b results/smoke/predictor_b_jittered.tsv --label-b FloorJittered \
  --score-column score \
  --boot-iters 1000 \
  --output "$OUT_DIR/primary_smoke.tsv"

echo
echo "--- Stratum 1 (H3) ---"
python -m src.analyse.run_stratum1 \
  --variants data/variant_list_frozen.tsv \
  --scores results/smoke/predictor_a_stratum1.tsv \
  --label FloorBenchmark \
  --score-column score \
  --boot-iters 10000 \
  --output "$OUT_DIR/stratum1_smoke.tsv"

echo
echo "--- Output files ---"
ls -1 "$OUT_DIR"
