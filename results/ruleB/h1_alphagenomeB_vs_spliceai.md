# Primary analysis — AlphaGenome-RuleB vs SpliceAI

## Stratum 2 coverage
- Frozen Stratum 2 variants: 52
- Scored by both predictors: 52
- Missing from AlphaGenome-RuleB: 0
- Missing from SpliceAI: 0
- Class balance in analysed set: 43 splice-affecting / 9 no-effect

## Pre-registered paired DeLong (H1 primary, two-sided, α=0.05)
- AUC AlphaGenome-RuleB: 0.9922
- AUC SpliceAI: 0.9819
- Paired difference (AUC_AlphaGenome-RuleB − AUC_SpliceAI): +0.0103
- Standard error: 0.0142
- 95% CI on paired difference: [-0.0176, +0.0382]
- z: 0.726
- Two-sided p: 0.4678
- H1 decision: **NOT SUPERIOR** (AlphaGenome-RuleB vs SpliceAI)

## Pre-registered non-inferiority (H2 secondary, margin δ=0.05)
- Lower bound of 95% CI on paired difference: -0.0176
- Non-inferiority margin: −0.05
- H2 decision: **NON-INFERIOR**

## Exploratory stratified-bootstrap AUC CIs (robustness display)
- AlphaGenome-RuleB: median 0.9948, 95% CI [0.9690, 1.0000]
- SpliceAI: median 0.9845, 95% CI [0.9432, 1.0000]
- 10000 iterations, seed 20260913 (OSF registration date). Stratified by outcome class.

## Provenance
- Test implementation: src.analyse.paired_delong (Sun & Xu 2014 fast DeLong).
- Frozen inputs: data/variant_list_frozen.tsv (hash-locked via data/HASH_MANIFEST.tsv).
- Pre-registration: OSF DOI 10.17605/OSF.IO/6PGX8, registered 2026-09-13.
