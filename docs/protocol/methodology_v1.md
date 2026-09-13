# Ranking residual normal splicing in CFTR variants: a benchmark of sequence-model junction predictions against experimental transcript measurements

**Version:** 1.0
**Date:** 2026-09-11
**Author:** Robert English
**Status:** Pre-registration methodology — frozen before first live AlphaGenome API call

---

## 1. Hypothesis

**H1 (primary):** In CFTR variants with published experimental measurements of residual normally-spliced transcript, [AlphaGenome](https://www.alphagenomedocs.com/) per-junction usage predictions rank-order variants by residual normal splicing better than [SpliceAI](https://github.com/Illumina/SpliceAI) Δ-scores.

**H2 (secondary):** The ranking performance of AlphaGenome is not statistically worse than [Pangolin](https://github.com/tkzeng/Pangolin) usage predictions.

**H0 (falsifier):** AlphaGenome AUC 95% CI overlaps SpliceAI AUC 95% CI → H1 rejected. Pangolin AUC point estimate exceeds AlphaGenome AUC and CIs do not overlap → H2 rejected.

Both directions publishable — this is an evidence-appraisal benchmark, not a positive-result-only study.

## 2. Scope

- The paper reports discrimination and rank-correlation performance only.
- No claim about clinical utility, prescribing, eligibility, modulator response, or patient management.
- Intended re-use: prioritisation of variants for functional experimental testing (theratyping pre-screen).

## 3. Population

**Inclusion:** CFTR variants (any coding class: canonical splice, splice-region, deep intronic, synonymous, missense-with-splice-effect) for which at least one of the following exists in a peer-reviewed publication:
- Numeric % of wild-type-length normally spliced transcript, with SD/SEM
- Binary classification (splice-affecting yes/no) from a systematic experimental characterisation
- Multi-organ or multi-cell-type residual measurement

**Exclusion:**
- Variants introduced by cell-line engineering that do not exist in patients
- Variants without HGVS c. notation traceable to GRCh38
- Variants whose only classification comes from in silico prediction

**Pre-declared strata:**
- **Stratum 1** (Spearman ρ, n ≈ 9–14): variants with quantitative % residual
- **Stratum 2** (AUC, n = 52): variants with binary residual/no-residual call
- **Stratum 3 (exploratory only):** CFTR2 non-canonical splice regulator adjudications (n ≈ 19)

**Effective n verified 2026-09-10:** Stratum 1 = 8–9 variants; Stratum 2 = 52 (Joynt) + ~15 (other sources). Harvest gate passed.

## 4. Predictors

### 4.1 AlphaGenome (index model)

- **Query mode:** [AlphaGenome Atlas](https://blog.google/innovation-and-ai/models-and-research/google-deepmind/alphagenome-atlas/) (precomputed SNVs) as first pass; live API for indels and complex alleles.
- **Output extracted:** per-junction usage predictions for the two junctions flanking each variant (donor and acceptor).
- **Score definition:** absolute change in predicted junction usage = |usage(REF) − usage(ALT)| at the canonical junction closest to the variant.
- **Cell context:** default lung/bronchial epithelial context; documented at pre-registration.
- **No parameter fitting:** raw predicted usage delta used directly; no calibration curve, no downstream model trained on outputs.

### 4.2 SpliceAI (comparator 1)

- Standard Illumina release, GRCh38; Δ-max across the four output channels (donor gain/loss, acceptor gain/loss) at ±500 nt from variant.

### 4.3 Pangolin (comparator 2)

- Standard installation; usage delta at the reference splice site closest to variant.
- Included as a demanding baseline: Pangolin is the harder comparison than SpliceAI, so an honest test requires it.

## 5. Outcome definitions

### 5.1 Primary outcome — Stratum 2 binary AUC (inferential)

- **Label:** 1 = splice-affecting (any documented deviation from normal splicing), 0 = no splice effect
- **Statistic:** AUROC with DeLong 95% CI
- **Comparison:** three models (AlphaGenome, SpliceAI, Pangolin) computed independently, then paired DeLong test for each pairwise AUC difference
- **H1 pass:** AlphaGenome AUC point estimate > SpliceAI AUC AND paired DeLong p < 0.05
- **H2 pass:** 95% CI for (AlphaGenome AUC − Pangolin AUC) lower bound > −0.05 (non-inferiority margin)

### 5.2 Secondary outcome — Stratum 1 rank correlation (descriptive)

- **Label:** % of wild-type-length normally spliced transcript (0–100)
- **Statistic:** Spearman ρ with 10,000-iteration bootstrap 95% percentile CI
- **Reported descriptively.** Given n ≈ 9, no significance test is claimed; point estimate + CI only. This descriptive-only framing is fixed in advance of any live prediction.

### 5.3 Expressivity outcome (descriptive)

Tally: which tools emit a residual-usage estimate at all, versus which emit only a binary spliceogenicity call.

## 6. Pre-registered analysis plan

1. Freeze variant list at OSF with SHA-256 hash before any prediction query.
2. Query all three models with identical variant list in a single batch.
3. Compute AUCs on Stratum 2 with paired DeLong tests.
4. If AlphaGenome AUC CI overlaps SpliceAI AUC CI: **H1 rejected** — report null.
5. If AlphaGenome AUC CI does not overlap SpliceAI AUC CI AND AlphaGenome > SpliceAI: **H1 supported**.
6. Test H2 against Pangolin (paired DeLong + non-inferiority margin).
7. Report Spearman ρ on Stratum 1 descriptively with bootstrap CI.
8. Report expressivity tally.
9. All figures pre-declared. No post-hoc stratification. Any added stratification is labelled exploratory and reported separately.

## 7. Data-leakage and training-contamination check

- **AlphaGenome training cutoff:** documented from [model card](https://huggingface.co/google/alphagenome-all-folds).
- **Contamination flag:** for each variant, publication date vs training cutoff logged. Report % of Stratum 2 published post-cutoff.
- **ClinVar overlap:** AlphaGenome training may have seen ClinVar annotations. ClinVar-annotated splice-region variants recorded as a separate exploratory stratum.
- **CFTR is a heavily studied gene:** training contamination is a genuine risk. The expressivity outcome (§5.3) partially insulates against this because usage magnitude is not directly annotated in most benchmark corpora.

## 8. Model-licence firewall

- **F1 — no downstream ML on outputs:** raw AUC and rank statistics only; no parameter fitting on AlphaGenome outputs.
- **F2 — no clinical-decision framing:** banned-word list applied to the full manuscript (`actionable`, `informs prescribing`, `eligible for`, `should receive`, `clinical utility`, `decision support`, `recommendation for`). Adversarial language audit before submission. Google's non-clinical statement quoted verbatim in Methods.
- **F3 — non-commercial academic use.**
- **F4 — no patient genome data:** published variants only.
- **F5 — licence clarification:** DeepMind courtesy notification filed regarding rank-order AUC computation on AlphaGenome outputs. If the clarification indicates that computed benchmark metrics fall outside permitted non-commercial research use, the paper is descriptive-only (expressivity tally in §5.3, no AUC in §5.1).

## 9. Falsification and stopping rules

- **Publish null if:** AlphaGenome AUC < SpliceAI AUC with non-overlapping CIs.
- **Publish expressivity-only report if:** DeepMind licence clarification denies rank-statistic use.
- **Publish feasibility census if:** Stratum 2 n < 30 after harvest freeze (n = 52 currently: not expected to trigger).
- **No p-hacking.** All statistics pre-declared. No re-runs with different cell contexts. No swapping of primary and secondary outcomes.

## 10. Timeline

| Week | Task |
|---|---|
| 1 | OSF pre-registration lodged (this document + variant harvest, SHA-256 hashed). Optional: extend harvest to Rincon, Sharma, Sosnay, CFTR2 records — additions logged as amendment before predictions run. |
| 2 | Query AlphaGenome Atlas + API; run SpliceAI and Pangolin locally. |
| 3 | Compute all statistics per §6. Generate figures per pre-declared plan. |
| 4 | Contact [October 2025 bioRxiv preprint authors](https://www.biorxiv.org/content/10.1101/2025.10.30.685113v1) — offer collaboration on residual-transcript characterisation. |
| 5–6 | Manuscript draft: Methods first, then Results, Introduction, Discussion. |
| 7 | Internal adversarial language audit; reviewer-hostile edit pass. |
| 8 | Preprint to bioRxiv/medRxiv + submission to *Pharmacogenomics J* or *Pharmacogenomics*. |

## 11. Novelty positioning

- **vs [Hajto 2026](https://www.nature.com/articles/s41397-026-00399-0) (541 PharmVar alleles, UK Biobank exomes):** does not cover CFTR splicing at all. Non-overlapping.
- **vs [Vitraag blog demo](https://www.vitraag.com/2025/06/26/pharmacogenomics-with-alphagenome/):** single CYP2D6\*4 case study, no cohort, no comparator, not peer-reviewed. Not prior art.
- **vs [Lin et al. bioRxiv Oct 2025](https://www.biorxiv.org/content/10.1101/2025.10.30.685113v1):** discusses two CFTR GT>GC variants as worked examples of graded residual, does not evaluate any AI model. Complementary — collaboration to be offered.
- **vs [APF2 2024](https://www.nature.com/articles/s41397-024-00338-x):** missense-focused, uses AlphaMissense not AlphaGenome, star-allele-focused. Non-overlapping.
- **vs [Joynt 2020](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009100):** provides the ground truth this paper benchmarks against. Foundation, not competition.

## 12. Target journals

1. *Pharmacogenomics Journal* (Nature) — best fit; publishes AI-variant benchmarks with clinical grounding
2. *Pharmacogenomics* (Future Medicine)
3. *Clinical Pharmacology & Therapeutics*
4. *Journal of Cystic Fibrosis* — CF-genetics specialist alternative
5. *Human Mutation* / *Human Genetics and Genomics Advances* — clinical genetics fallback
