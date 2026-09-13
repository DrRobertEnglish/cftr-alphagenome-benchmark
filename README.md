# cftr-alphagenome-benchmark

A pre-registered computational benchmark of three in silico splice-prediction models — [AlphaGenome](https://www.alphagenomedocs.com/), [SpliceAI](https://github.com/Illumina/SpliceAI), and [Pangolin](https://github.com/tkzeng/Pangolin) — on published CFTR splice variants with experimentally measured ground-truth splicing outcomes.

## Pre-registration

- **OSF project:** [osf.io/5hvgf](https://osf.io/5hvgf/)
- **OSF registration DOI:** [10.17605/OSF.IO/6PGX8](https://doi.org/10.17605/OSF.IO/6PGX8)
- **Registration date:** 2026-09-13 (public, no embargo)

The variant list and full methodology were frozen and hashed on OSF **before any predictor query**. This repository holds the analysis code that will be executed against that frozen protocol. Any change to the variant list or methodology after the registration date requires a versioned OSF amendment with a fresh SHA-256 hash before further predictor queries.

## Scope

- **What this is:** a paired discrimination and rank-correlation benchmark of three publicly available in silico splice predictors on a fixed set of published CFTR variants.
- **What this is not:** a clinical decision-support system. No prescribing, patient-selection, treatment-eligibility, procurement, or commissioning claim is made. See the pre-registration for the banned-word list enforced across the manuscript.

## Hypotheses (pre-registered)

- **H1 (primary):** AlphaGenome per-junction usage predictions discriminate splice-affecting from no-effect variants better than SpliceAI Δ-scores on the Joynt et al. 2020 Stratum 2 (n = 52).
  Test: paired DeLong on ROC-AUC, two-sided α = 0.05.
- **H2 (secondary, non-inferiority):** AlphaGenome is not meaningfully worse than Pangolin on the same stratum.
  Test: paired DeLong 95% CI on the AUC difference, non-inferiority margin δ = 0.05.
- **H3 (secondary, descriptive):** AlphaGenome per-junction usage delta correlates monotonically with published percentage residual normally spliced transcript on Stratum 1 (n ≈ 9).
  Statistic: Spearman ρ with 10,000-iteration bootstrap 95% percentile CI. No p-value, no pass/reject decision.

## Repository layout

```
cftr-alphagenome-benchmark/
├── README.md                       — this file
├── LICENSE                         — MIT
├── pyproject.toml                  — Python environment specification
├── .gitignore
├── data/
│   ├── variant_list_frozen.tsv     — the frozen CFTR variant list (SHA-256 hashed on OSF)
│   └── ground_truth/               — extracted per-variant ground-truth measurements
├── src/
│   ├── harvest/                    — code that produced the frozen variant list (audit trail)
│   ├── predict/                    — per-predictor batch scripts (AlphaGenome, SpliceAI, Pangolin)
│   ├── controls/                   — positive and negative control benchmarks
│   ├── analyse/                    — paired DeLong, Spearman bootstrap, permutation tests
│   └── figures/                    — figure generation
├── results/                        — populated at runtime, never committed with API-returned scores until predictions are frozen
├── docs/
│   └── protocol/                   — mirror of the OSF pre-registration files
└── tests/                          — unit tests for the analysis code, run before the confirmatory batch
```

## Model licence positioning

AlphaGenome is used under its [non-commercial academic research licence](https://storage.googleapis.com/alphagenome/terms/AlphaGenome-Model-Terms-of-Use.pdf). This repository operates within that scope:

- No downstream ML model is trained on AlphaGenome outputs.
- No clinical-decision framing is asserted in code or documentation.
- No patient genomes are queried against AlphaGenome; all inputs are peer-reviewed public literature variants.

A courtesy notification has been sent to the AlphaGenome team at `alphagenome@google.com`. If a licence-scope response before submission indicates that rank-order AUC computation on AlphaGenome outputs falls outside permitted use, the pre-registered contingency is invoked: all AlphaGenome AUC and paired-test results are withheld and the manuscript is reduced to an expressivity tally.

## Reproducibility

- All random seeds are fixed and documented in the analysis code.
- All per-variant predictor outputs will be released as a supplementary table on submission, so any reader can re-derive every reported statistic without re-querying the models.
- Two independent implementations of the paired DeLong test are run (R `pROC` and a Python port of Sun & Xu 2014); any disagreement > 0.005 in AUC or > 0.01 in the p-value is disclosed.

## Status

- [x] OSF pre-registration submitted and public (2026-09-13)
- [x] Repository created (2026-09-13)
- [x] MIT licence
- [ ] Python environment locked
- [ ] Positive and negative control benchmarks
- [ ] Predictor batches (AlphaGenome via Atlas + API, SpliceAI, Pangolin)
- [ ] Statistical analysis
- [ ] Manuscript preprint
- [ ] Journal submission

## Citation

If you use this benchmark or its methodology, please cite the OSF pre-registration:

> English R (2026). *Ranking residual normal splicing in CFTR variants: a benchmark of sequence-model junction predictions against experimental transcript measurements.* OSF Pre-registration. https://doi.org/10.17605/OSF.IO/6PGX8

A preprint DOI will be added to this README on bioRxiv submission.

## Contact

Author: Robert English ([@DrRobertEnglish](https://github.com/DrRobertEnglish))
Correspondence: via the [OSF project page](https://osf.io/5hvgf/) or by opening an issue on this repository.
