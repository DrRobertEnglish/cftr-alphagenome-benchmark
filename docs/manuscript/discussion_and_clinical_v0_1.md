# Discussion, clinical implications, and next steps (draft v0.1)

Sections designed to follow the Results in the manuscript. Written in the same neutral academic voice as the Results draft, respecting the licence-scope restriction to summary statistics and the explicit non-clinical, non-prescribing framing recorded in `docs/protocol/AMENDMENTS.md` §A1 and `docs/protocol/methods_licence_paragraph.md`.

---

## Discussion

### Principal findings

On a frozen, pre-registered CFTR splice benchmark, AlphaGenome achieved perfect binary discrimination on Stratum 2 (ROC-AUC 1.000, correctly rank-ordering all 43 splice-affecting variants above all 9 no-effect variants), the only one of the three tested predictors to do so; SpliceAI 1.3.1 (AUC 0.982) and Pangolin (AUC 0.984) each mis-ranked at least one variant. Under the pre-specified paired DeLong test the +0.018 and +0.016 AUC advantages did not reach statistical significance at n = 52, so the formal superiority claim is not made; the pre-registered non-inferiority claim (δ = 0.05) held against both comparators. On the small quantitative Stratum 1 (n = 8), AlphaGenome's Spearman rank correlation with the fraction of normally spliced transcript was approximately twice the magnitude of the two comparators (ρ = −0.66 vs −0.32 and −0.29), in the pre-specified negative direction. All three predictors comfortably exceeded the pre-registered distance-to-canonical-splice-site floor (AUC 0.683).

### Interpretation

AlphaGenome ranked every splice-affecting variant above every no-effect variant on Stratum 2 (43 affecting, 9 no-effect), yielding a point-estimate ROC-AUC of exactly 1.000. Neither SpliceAI (AUC 0.982) nor Pangolin (AUC 0.984) achieved this: each mis-ranked at least one variant. Under paired DeLong testing, the +0.018 and +0.016 point-estimate advantages do not reach statistical significance at n = 52 because the standard error of the paired AUC difference is 0.015 when both scores are near the ceiling. The formal superiority claim is therefore not made. This is a distinction between a statistical decision under a pre-registered test and the raw ranking behaviour of the three predictors; it should not be read as evidence that AlphaGenome performed no better on this benchmark. On the pre-registered non-inferiority test at δ = 0.05 AlphaGenome cleared both comparators, and on unadjusted rank-order performance AlphaGenome was the only predictor to achieve perfect discrimination.

The result also has to be read against a well-recognised ceiling effect: CFTR splice variants in Joynt et al. and comparable curated series are heavily enriched for canonical splice-site disruptions, which every modern sequence-based predictor identifies well. Under these conditions, discrimination benchmarks tend to saturate, and a benchmark's ability to separate top-tier predictors from each other is inherently limited. That AlphaGenome nonetheless produced a clean sweep where SpliceAI and Pangolin did not is a positive signal within that ceiling, even without a formal superiority test being reachable at n = 52.

Stratum 1 tells a different and complementary story. Discrimination is not the question — every variant in this stratum is known to alter splicing to some degree. The question is quantitative fidelity: does a higher predictor score track a lower fraction of normally spliced transcript? AlphaGenome's Spearman ρ of −0.66 is roughly twice the magnitude of either comparator (−0.32 and −0.29) in the pre-specified negative direction, though the eight-variant sample size gives a wide bootstrap 95% confidence interval that includes zero. We report this as a descriptive finding per the frozen protocol and do not draw an inferential conclusion.

Taken together, the results are consistent with AlphaGenome being at least as accurate as established CFTR splice predictors on canonical splice-site variation — with a perfect Stratum 2 point estimate that neither comparator achieved — and with a potentially stronger relationship to quantitative splicing outcome that warrants confirmation in larger quantitative datasets.

### Comparison with the AlphaGenome preprint

AlphaGenome's original evaluation on ClinVar splicing variants reported a stronger separation from SpliceAI/Pangolin than this study shows on CFTR. Two features of the CFTR benchmark likely contribute to the narrower gap observed here: the enrichment for canonical splice-site variants (which are already near-solved), and the modest sample size (n = 52 in Stratum 2, n = 8 in Stratum 1). This does not contradict the AlphaGenome authors' findings; it constrains where those findings apply.

### Limitations

Four limitations bound the interpretation of these results.

- **Sample size.** Stratum 2 (n = 52) sits at the discrimination ceiling for the class of variants included, limiting the paired DeLong test's power to detect small AUC gaps. Stratum 1 (n = 8) is descriptive by design; no formal inferential claim is made.
- **Variant distribution.** The frozen benchmark inherits the composition of published CFTR splice-variant series, which favour canonical splice-site variants over deep-intronic and exonic-splicing-element variants. Predictor rank-order on this benchmark may not generalise to the deep-intronic tail.
- **Assay heterogeneity.** Stratum 1's quantitative measurements combine minigene and endogenous-RNA assays from different laboratories with different reference cell lines and quantitation methods. Direct comparability across the eight variants is therefore imperfect and inflates the residual variance of any rank correlation.
- **Predictor version drift.** AlphaGenome, SpliceAI, and Pangolin are all versioned software. Scores in this benchmark were produced against pinned versions (AlphaGenome 0.9.0, SpliceAI 1.3.1, Pangolin at a specified commit); subsequent releases may yield different results.

### Compliance and scope

This study reports summary statistics computed from AlphaGenome output (ROC-AUC, DeLong confidence intervals, Spearman ρ), not derivative models trained on AlphaGenome output. Publication is within the academic non-commercial research scope of the AlphaGenome Model Terms of Use, with the conspicuous notice, non-clinical framing, and Section 8 disclaimer supplied as specified in the licence documentation accompanying the repository.

---

## Clinical implications

The results below sit deliberately at the level of scientific hypotheses about how splice-variation modelling could inform CFTR care in future. **They are not clinical recommendations, do not constitute decision support, and do not license any change to current practice.** All CFTR variant interpretation must continue to follow the ACMG/AMP framework, CFTR2 clinical annotations, cystic fibrosis speciality guidelines, and local multidisciplinary team review.

### Why splice fidelity may matter for CFTR therapeutics

Cystic fibrosis modulator eligibility is currently annotated at the level of specific CFTR variants (e.g. F508del for elexacaftor/tezacaftor/ivacaftor, R117H for ivacaftor, and the CFTR2-listed variants extended by the 2020 US label expansion). The functional principle underneath these annotations is CFTR protein availability at the apical membrane. For splice variants, that quantity is upstream-limited by the fraction of full-length, normally spliced mRNA the variant permits.

Two implications follow.

- **Rank-ordering CFTR splice variants by residual normal-transcript fraction is a scientifically defensible way to prioritise which splice variants might benefit from potentiator therapy** (variants that preserve a substantial fraction of correctly spliced transcript produce some functional CFTR that a potentiator could act on) versus corrector or read-through strategies (variants that essentially abolish normal transcript). This is a research hypothesis, not a treatment algorithm.
- **Quantitative splice predictors could serve as a triage layer for laboratory splicing assays** in variants of uncertain significance. Where predictor score is very low (variant likely to preserve normal splicing) or very high (variant likely to abolish it), the pre-test probability shifts materially and the diagnostic yield of a functional assay changes.

### How the current results speak to these hypotheses

- **On Stratum 2 (canonical splice-site variants):** all three tested predictors perform well; the practical translational question is not "which predictor" but "does prediction accuracy at this end of the spectrum change management." At canonical splice sites, the answer today is essentially no — a canonical +1/+2 or −1/−2 variant is already classified as loss-of-function by clinical criteria without needing a predictor score.
- **On Stratum 1 (quantitative residual-splicing variants):** AlphaGenome's stronger monotone relationship with residual normal-transcript fraction, if confirmed in larger datasets, would support the use of quantitative splice scoring as a research triage input for questions such as: which non-canonical CFTR splice variants preserve enough correctly spliced transcript for potentiator response to be plausible; which variants should be prioritised for laboratory splicing assays; which variants might be candidates for splice-modulating therapies (antisense oligonucleotides, small-molecule splice modulators).

### What this study does not support

- Any change to current CFTR variant classification per ACMG/AMP or ClinGen SVI splicing subgroup rules.
- Any change to current CF modulator eligibility.
- Use of AlphaGenome, SpliceAI, or Pangolin scores as a standalone basis for clinical decisions.
- Any inference about splice fidelity for non-CFTR splice variants beyond the extent of validated cross-gene generalisation shown by each predictor's authors.

### Where this study might inform CF research programmes

- **Variant curation.** As a supporting-strength research input alongside laboratory splicing assay evidence, per the ClinGen SVI splicing subgroup's framework for using computational splice predictions.
- **Assay prioritisation.** When resources for minigene or endogenous-transcript splicing assays are limited, quantitative splice-score rank order could inform which variants are studied first.
- **Trial design.** For emerging splice-modulating therapies, quantitative predictor scores may inform the selection of candidate variants for early-phase clinical evaluation, complementing (not replacing) direct functional characterisation.

---

## Next steps

### Data expansion

The single most impactful improvement to this benchmark would be a substantially larger Stratum 1. As of September 2026, no public database (MaveDB, MPSA benchmarks, or targeted CFTR minigene collections) provides more than a few dozen quantitative CFTR splicing measurements, and the largest current MPSA studies (POU1F1, RON, FAS, WT1 by Smith and Kitzman 2023; COMPASS by Koplik et al. 2025) do not focus on CFTR. Three tractable expansion routes exist.

- **Aggregated minigene meta-benchmark.** Systematic literature harvest of published CFTR minigene and endogenous-RNA splicing assays with quantitative readouts (Sterrantino et al., Bergougnoux et al., Igreja et al., and comparable sources) into a per-source-stratified quantitative benchmark. Realistic yield: 30–60 quantitative variants. Would be published as a versioned amendment (`docs/protocol/AMENDMENTS.md` §A2) with pre-specified per-source random effects.
- **COMPASS re-analysis.** The Koplik et al. 2025 dataset covers 87,546 variants across more than 1,700 genes. CFTR-specific coverage is not published in the abstract; a targeted analysis of the COMPASS supplementary tables would either yield a substantially larger Stratum 1 or clarify the current data ceiling.
- **De novo CFTR MPSA.** A saturation MPSA focused on the CFTR gene, ideally covering canonical splice sites, deep-intronic elements, and exonic splicing elements together, would allow first-principles benchmarking without the assay-heterogeneity limitations of aggregation. This would be a multi-lab collaboration, not a single-doctor project.

### Methodological extensions

- **Regulatory-region variants.** Extending the frozen benchmark to include deep-intronic and exonic-splicing-enhancer variants, where all three predictors are known to be less accurate than at canonical splice sites, would provide a more discriminating test of quantitative fidelity than the current variant distribution allows.
- **Cross-predictor calibration.** Formal calibration analysis (Brier score, calibration curves, isotonic recalibration) on a larger dataset would allow rank-preserving score comparisons that this study's small quantitative sample size does not.
- **Predictor version tracking.** A live continuous-integration workflow that re-runs the benchmark against each new AlphaGenome, SpliceAI, and Pangolin release would let the field track predictor drift on a stable reference set. The current repository is already structured to support this.

### Publication timeline

- Preprint (biorxiv or medrxiv) deposit as soon as the manuscript is complete.
- Peer-reviewed submission to a splice-biology, computational-biology, or CF-focused journal in parallel with any Stratum 1 expansion work.
- A follow-up paper covering the aggregated minigene meta-benchmark, if the harvest yields sufficient variants, or reporting the negative result of that harvest.

### Open science commitments

- All variant lists, predictor outputs, analysis code, figures, and manuscript drafts are hash-locked and version-controlled in the accompanying public repository.
- Pre-registration remains fixed at [osf.io/5hvgf](https://osf.io/5hvgf/) (DOI [10.17605/OSF.IO/6PGX8](https://doi.org/10.17605/OSF.IO/6PGX8)); any expansion is registered as a versioned amendment before the corresponding predictions are generated.
- Predictor outputs are shared under the licence conditions of each predictor, with AlphaGenome outputs subject to the AlphaGenome Output Terms of Use notice recorded in the repository.
