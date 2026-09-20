# Primary analysis — AlphaGenome vs SpliceAI

## Stratum 2 coverage
- Frozen Stratum 2 variants: 52
- Scored by both predictors: 52
- Missing from AlphaGenome: 0
- Missing from SpliceAI: 0
- Class balance in analysed set: 43 splice-affecting / 9 no-effect

## Pre-registered paired DeLong (H1 primary, two-sided, α=0.05)
- AUC AlphaGenome: 1.0000
- AUC SpliceAI: 0.9819
- Paired difference (AUC_AlphaGenome − AUC_SpliceAI): +0.0181
- Standard error: 0.0151
- 95% CI on paired difference: [-0.0115, +0.0477]
- z: 1.198
- Two-sided p: 0.2309
- H1 decision: **NOT SUPERIOR** (AlphaGenome vs SpliceAI)

## Pre-registered non-inferiority (H2 secondary, margin δ=0.05)
- Lower bound of 95% CI on paired difference: -0.0115
- Non-inferiority margin: −0.05
- H2 decision: **NON-INFERIOR**

## Exploratory stratified-bootstrap AUC CIs (robustness display)
- AlphaGenome: median 1.0000, 95% CI [1.0000, 1.0000]
- SpliceAI: median 0.9845, 95% CI [0.9432, 1.0000]
- 10000 iterations, seed 20260913 (OSF registration date). Stratified by outcome class.

## Provenance
- Test implementation: src.analyse.paired_delong (Sun & Xu 2014 fast DeLong).
- Frozen inputs: data/variant_list_frozen.tsv (hash-locked via data/HASH_MANIFEST.tsv).
- Pre-registration: OSF DOI 10.17605/OSF.IO/6PGX8, registered 2026-09-13.
