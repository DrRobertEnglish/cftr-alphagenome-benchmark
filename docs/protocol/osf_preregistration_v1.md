# CFTR AlphaGenome splice benchmark — pre-registration v1.0

**Timestamp (UTC):** 2026-09-11T15:44:07Z
**Timestamp (Europe/London):** 2026-09-11 16:44 BST
**Author:** Robert English
**Status:** Registered before any live AlphaGenome API call

---

## Study Information

### Title

Ranking residual normal splicing in CFTR variants: a benchmark of sequence-model junction predictions against experimental transcript measurements

### Research question

In CFTR variants with published experimental measurements of splicing outcome, do [AlphaGenome](https://www.alphagenomedocs.com/) per-junction usage predictions rank-order variants by residual normal splicing better than [SpliceAI](https://github.com/Illumina/SpliceAI) Δ-scores, and no worse than [Pangolin](https://github.com/tkzeng/Pangolin) usage predictions?

### Hypotheses

- **H1 (primary):** AlphaGenome per-junction usage predictions rank-order CFTR splice variants by residual normal splicing better than SpliceAI Δ-scores.
- **H2 (secondary):** AlphaGenome's ranking performance is not statistically worse than Pangolin's.
- **H0 (falsifier):** overlapping AUC 95% CIs between AlphaGenome and SpliceAI → H1 rejected. Both directions publishable.

---

## Design plan

### Study type

Retrospective computational benchmark of three published in silico splice-prediction tools against published experimental splicing measurements. No new experimental data are collected.

### Blinding

Not applicable — all inputs and all comparator outputs are computed programmatically on a frozen variant list.

### Study design

Three tools (AlphaGenome, SpliceAI, Pangolin) applied to the same frozen list of CFTR variants. Per-variant model outputs compared against per-variant published experimental measurements using pre-registered discrimination and rank-correlation statistics. All statistics pre-declared before any live prediction is run.

---

## Sampling plan

### Existing data

Inputs are published CFTR variants with experimentally measured splicing outcomes. Data existed prior to this pre-registration.

### Existing-data safeguards

- No AlphaGenome, SpliceAI, or Pangolin prediction has been executed by the author on the variant list prior to this registration.
- The frozen variant list is stored with SHA-256 hash (below) at registration timestamp.
- No summary statistic on the model outputs has been computed prior to this registration.

### Data source

- [Masvidal 2014, EJHG](https://www.nature.com/articles/ejhg2013238)
- [Joynt 2020, PLOS Genetics](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009100)
- [Chiba-Falek 1999, AJRCCM](https://academic.oup.com/ajrccm/article/159/6/1998/8530084)
- [Deletang 2022, Gene Ther](https://www.nature.com/articles/s41434-022-00347-0)
- Additional sources (Rincon, Sharma, Sosnay, CFTR2) may be added as a formal amendment before predictions are run.

### Sample size

- **Stratum 1** (quantitative % residual): n ≈ 8–9 unique variants (+ 5 additional multi-organ data points for c.3717+12191C>T)
- **Stratum 2** (binary splice-affecting): n = 52 variants ([Joynt 2020](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009100))

Pre-registered harvest threshold: n ≥ 20 continuous OR n ≥ 30 binary. Both met.

---

## Variables

### Predictors

1. **AlphaGenome** per-junction usage delta = |usage(REF) − usage(ALT)| at the canonical junction closest to the variant, default lung/bronchial epithelial cell context, retrieved via [AlphaGenome Atlas](https://blog.google/innovation-and-ai/models-and-research/google-deepmind/alphagenome-atlas/) (SNVs) and live API (indels, complex alleles).
2. **SpliceAI** Δ-max across the four output channels (donor gain/loss, acceptor gain/loss) at ±500 nt from variant, standard Illumina release, GRCh38.
3. **Pangolin** usage delta at the reference splice site closest to the variant, standard installation.

### Outcomes

- **Primary (inferential):** binary label 1 = splice-affecting (any documented deviation from normal splicing), 0 = no splice effect.
- **Secondary (descriptive):** continuous % of wild-type-length normally spliced transcript (0–100).
- **Expressivity (descriptive):** whether each tool emits a residual-usage estimate vs a binary spliceogenicity call only.

---

## Analysis plan

### Primary inferential analysis (Stratum 2, n = 52)

- Statistic: AUROC with DeLong 95% CI
- Comparison: paired DeLong test for pairwise AUC difference (AlphaGenome vs SpliceAI; AlphaGenome vs Pangolin)
- **H1 pass criterion:** AlphaGenome AUC point estimate > SpliceAI AUC AND paired DeLong p < 0.05
- **H2 pass criterion:** 95% CI for (AlphaGenome AUC − Pangolin AUC) lower bound > −0.05 (non-inferiority margin)

### Secondary descriptive analysis (Stratum 1, n ≈ 9)

- Statistic: Spearman ρ with 10,000-iteration bootstrap 95% percentile CI
- Reported descriptively. Given n ≈ 9, no significance test is claimed; point estimate + CI only.

### Expressivity tally (descriptive)

Count and enumerate which tools emit a continuous residual-usage estimate versus a binary spliceogenicity call.

### Software

- [`pROC`](https://cran.r-project.org/web/packages/pROC/) (R) or hand-implemented DeLong in Python
- [`scipy.stats.spearmanr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html) + [`scipy.stats.bootstrap`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bootstrap.html)
- Analysis code published alongside preprint.

### Inference criteria

- All statistics pre-declared. No post-hoc stratification. Any stratification added after data inspection is labelled exploratory and reported separately.
- No p-hacking; no re-runs with different cell contexts.
- Falsification / stopping rules stated in the accompanying methodology document (§9).

### Data leakage

- AlphaGenome training cutoff logged from [model card](https://huggingface.co/google/alphagenome-all-folds).
- Per-variant publication date vs training cutoff logged; % of Stratum 2 published post-cutoff reported.
- ClinVar-annotated variants recorded as a separate exploratory stratum.

### Missing data

Variants with unresolvable HGVS c. notation, or with no traceable primary source, are excluded before the freeze. After the freeze, no further exclusions permitted.

### Exploratory analyses

CFTR2 non-canonical splice regulator adjudications (n ≈ 19) — Stratum 3 exploratory only. Reported separately from confirmatory results.

---

## Frozen protocol files (SHA-256)

| File | SHA-256 |
|---|---|
| variant_harvest_v1_public.md | `51dcde049fdb943abf7caad3a4bc35fd0694d48dbc8d555d5b35c5c0ba011334` |
| methodology_v1_public.md | `32280a7a725b5cd25defa5fc5cab5ec7b66d008b1a06a7ee32dfd0bdc88f9893` |

Any change to either file after this timestamp constitutes a protocol amendment and must be registered separately with a new hash.

---

## Model-licence firewall

- Raw AUC and rank statistics only; no downstream model trained on AlphaGenome outputs.
- Banned-word list applied to full manuscript to prevent clinical-decision framing.
- Google's non-clinical statement quoted verbatim in Methods.
- Non-commercial academic use only.
- DeepMind courtesy notification filed regarding rank-order AUC computation on AlphaGenome outputs. If the response indicates that benchmark metrics fall outside permitted non-commercial research use, the paper is descriptive-only (expressivity tally, no AUC).

---

## Registration metadata

- OSF project URL: https://osf.io/5hvgf/
- OSF registration DOI: [10.17605/OSF.IO/6PGX8](https://doi.org/10.17605/OSF.IO/6PGX8)
- Registration date: 2026-09-13
- Frozen variant list attached as supplementary file
- Full methodology attached as supplementary file

---

*Register this file on OSF before running any live AlphaGenome API call. Once registered, log the OSF DOI at the top of the manuscript's Methods section.*
