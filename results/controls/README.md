# Controls — results

## Positive control: distance-to-canonical-splice-site

- **Input:** `data/variant_list_frozen.tsv` (60 rows; 52 Joynt Stratum 2 + 8 Stratum 1) resolved to GRCh38 via VariantValidator against NM_000492.4 (CFTR MANE Select).
- **Reference:** `data/reference/CFTR_ENST00000003084_11_exons_GRCh38.tsv` (27 exons; 52 canonical splice sites, forward strand).
- **Scoring rule:** `1 / (1 + d)` where `d` is the absolute bp distance from a variant's genomic position to the nearest canonical exon boundary.
- **Command:** `python -m src.controls.splice_site_distance_benchmark --variants data/variant_list_frozen.tsv --transcript data/reference/CFTR_ENST00000003084_11_exons_GRCh38.tsv --output results/controls/splice_site_distance_scores.tsv`
- **Result on Joynt Stratum 2 (n=52; 43 splice-affecting / 9 no-effect):** ROC-AUC = **0.6835**.

## Negative control: label-permutation null

- **Command:** `python -m src.controls.label_permutation_null --variants data/variant_list_frozen.tsv --scores results/controls/splice_site_distance_scores.tsv --iters 1000 --seed 20260913 --stratum 2 --output results/controls/permutation_null.tsv`
- **Seed:** 20260913 (OSF registration date).
- **Iterations:** 1000.
- **Null distribution of ROC-AUC on Stratum 2:** median 0.500, 95% CI [0.296, 0.698], SD 0.106.

## Interpretation

The floor-benchmark AUC of 0.6835 sits at approximately the 97th percentile of the null distribution — inside the 95% null envelope. On this particular dataset, with 52 variants weighted toward canonical splice sites (43/52 splice-affecting), raw distance-to-nearest-canonical-splice-site is barely better than chance.

This is the desired shape for a floor benchmark: any real predictor (SpliceAI, Pangolin, AlphaGenome) that clears clinical utility should beat it comfortably, and the paired-DeLong test in `src/analyse/paired_delong.py` will read exactly this scoring file when comparing predictors on the same 52 variants.

Wide null tails (SD 0.106) reflect the small effective sample size and the 43:9 class imbalance: any single-variant swap between classes moves the AUC substantially. The pre-registered analysis plan already accounts for this by using paired DeLong (each predictor scored on the identical variant set) rather than independent comparisons.
