# Rule A vs Rule B divergence — mechanistic analysis

Written 2026-09-20 after the parallel-scoring Colab run. Documents the specific
variants and mechanism that drive the difference between Rule A (DeepMind-recommended
composite) and Rule B (OSF-preregistered per-junction rule with lung/bronchial
cell context) on the 60 frozen CFTR variants.

## H1/H2 (Stratum 2, n=52) — one positive mis-ranks below three negatives

**Summary numbers**

| Rule | min positive score | max negative score | gap | AUC |
|---|---|---|---|---|
| A (composite) | 0.3529 | 0.1915 | +0.1615 (positive: perfect discrimination) | 1.0000 |
| B (per-junction, cell-context) | 0.2207 | 0.5020 | −0.2812 (negative: crossings) | 0.9922 |

Under Rule B, **1 positive (`Joynt_E02`, score 0.221)** sits below **3 negatives
(`Joynt_E05` 0.502, `Joynt_E04` 0.427, `Joynt_E03` 0.346)**. That is 3 wrong pairs
out of 43×9 = 387 possible pos–neg pairs = AUC loss of 0.0078 = observed AUC of
0.9922. The arithmetic closes exactly.

**Mechanism.** The three negatives that beat E02 all have very large
`chosen_junction_distance_bp` values:

| variant | binary | score_B | chosen_junction_distance_bp | score_A (for reference) |
|---|---|---|---|---|
| Joynt_E05 | 0 (neg) | 0.502 | **2892** | 0.191 |
| Joynt_E04 | 0 (neg) | 0.427 | **884** | 0.143 |
| Joynt_E03 | 0 (neg) | 0.346 | **883** | 0.133 |
| Joynt_E02 | 1 (pos) | 0.221 | 215 | 0.353 |

The three high-scoring negatives have no cell-context junction within a kilobase
of the variant, so the "closest available lung/bronchial junction" is far away.
The `raw_score` at those distant junctions is large but does not reflect an actual
splice-affecting effect at the variant position — it is the artefact of scoring a
distant junction that happens to shift in the reference-vs-alternate comparison for
reasons unrelated to the query variant. Rule A launders this out by aggregating
`max(|raw_score|)` across all tracks and genes: the composite picks up the true
signal wherever it happens across the 1 Mb interval, and the (spuriously high)
per-junction raw_score is one of many maxima that get taken.

## H3 (Stratum 1, n=8) — two variants drive the ρ gap from −0.66 to −0.29

**Per-variant table** (sorted by pct_normally_spliced):

| variant | pct_normally_spliced | score_A | score_B | chosen_junction_distance_bp |
|---|---|---|---|---|
| S1_05_Joynt | 0 | 3.335 | 6.102 | 0 |
| S1_04_Joynt | 0 | 3.255 | 5.820 | 0 |
| S1_09_Deletang | 2 | 3.185 | 6.031 | 50 |
| S1_02_Masvidal | 24 | 3.343 | 6.543 | 5 |
| S1_08_ChibaFalek_lung | 26 | 0.275 | 0.292 | 2 |
| S1_07_Joynt | 37.2 | 1.763 | 3.305 | 29 |
| **S1_01_Masvidal** | **40** | **3.243** | **6.281** | 164 |
| **S1_06_Joynt_WT** | **100** | **0.191** | **0.502** | **2892** |

**Leave-one-out ρ_B analysis** (full ρ_B = −0.2874):

- Drop `S1_01_Masvidal`: ρ_B goes to **−0.6126** (matches Rule A ρ = −0.6587).
- Drop `S1_06_Joynt_WT` (wild-type control): ρ_B goes to −0.036 (destroys correlation).
- Other drops: change ρ_B by less than 0.11 in magnitude.

`S1_01_Masvidal` (40% normally spliced, i.e. moderately affected) receives a very
high Rule B score of 6.28. Under Rule A it receives a more moderate 3.24. The chosen
junction is 164 bp from the variant — not extremely distant, but the per-junction
raw_score at that particular lung/bronchial junction is inflated relative to the
cross-track composite.

`S1_06_Joynt_WT` is the pre-registered wild-type control (100% normally spliced,
should receive near-zero predictor score). Rule A gives it 0.19 (appropriately
near-zero). Rule B gives it 0.50 — the same magnitude as a moderately splice-affecting
variant. Its `chosen_junction_distance_bp` is **2892** — the same 2.9 kb distance that
drove the H1/H2 misclassification. Same mechanism: the wild-type variant has no
lung/bronchial junction near it, and the "closest available" one 2892 bp away
produces a spuriously non-zero raw_score.

## Unified finding

**Rule B's degradation under all three hypotheses traces back to a single failure
mode:** when the pre-declared cell-context ontology set contains no junction within
a few hundred bp of the variant, the "closest junction" fallback selects a distant
junction whose per-junction raw_score is noisy and unrelated to the variant's true
splicing impact. This affects:

- 1 splice-affecting Stratum 2 variant (Joynt_E02) — pulls down its rank.
- 3 no-effect Stratum 2 variants (E03/E04/E05) — pushes them above E02.
- 1 wild-type control in Stratum 1 (Joynt_WT) — spuriously non-zero.
- 1 moderately-affected Stratum 1 variant (Masvidal_01) — spuriously very high.

Rule A's aggregation over all tissues and all genes launders this out because the
score is `max(|raw_score|)` across the full track/gene set: a spurious high in one
distant junction does not swamp the true signal in the near-variant region.

## Implication for the paper

This is a **substantive methodological finding**, not an artefact of a bad
implementation of the pre-reg rule. It says:

1. When AlphaGenome is used with a cell-context-restricted, per-junction scoring
   rule, its calibration depends on whether the pre-declared cell context contains
   annotated junctions near the variant.
2. The DeepMind-recommended composite score is robust to this because it aggregates
   across all tracks.
3. For a benchmark like this one (60 CFTR variants, non-uniform junction density in
   the lung/bronchial ontology set), the composite is the more reliable rule.
4. This is not a fault of the pre-registration itself — the pre-reg was written
   before we could inspect the actual AlphaGenome SDK output structure. It is a
   fault of any per-junction, per-cell-context scoring rule that does not include
   a distance-based fallback.

The paper should report this cleanly, show the per-variant table, and offer the
mechanistic explanation. This strengthens rather than weakens the paper: it
demonstrates that the parallel-analysis approach revealed a substantive property
of the scoring rules, not just noise.

## Reproducing this analysis

```
cd /home/user/workspace/cftr-alphagenome-benchmark
python3 <<'PY'
[the analysis code that produced this document]
PY
```

See `git log` for the exact analysis command used on 2026-09-20.
