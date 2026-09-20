# Discussion, clinical implications, and next steps (draft v0.2)

Sections designed to follow the Results in the manuscript. Written in the same neutral academic voice as the Results draft, respecting the terms-of-use restriction to summary statistics and the explicit non-clinical, non-prescribing framing recorded in `docs/protocol/AMENDMENTS.md` §A1 and `docs/protocol/methods_licence_paragraph.md`.

---

## Discussion

### Principal findings

On a frozen, pre-registered CFTR splice benchmark, AlphaGenome was non-inferior to SpliceAI 1.3.1 and Pangolin on Stratum 2 (n = 52) under the pre-specified paired DeLong non-inferiority test at margin δ = 0.05, and the non-inferiority conclusion held under both AlphaGenome scoring rules run as parallel primary analyses (protocol Amendment #2): the DeepMind-recommended merged splicing composite (Rule A) and the OSF-preregistered per-junction score in a lung/bronchial cell context (Rule B). The pre-registered superiority test did not clear its threshold under either rule (ΔAUC vs SpliceAI +0.018, p = 0.231 under Rule A; +0.010, p = 0.468 under Rule B). Discrimination itself differed by rule: under Rule A, AlphaGenome achieved perfect Stratum 2 rank separation (ROC-AUC = 1.000), the only one of the three tested predictors to correctly rank-order every splice-affecting variant above every no-effect variant; under Rule B, AlphaGenome achieved ROC-AUC 0.9922, its single AUC loss driven by one splice-affecting variant misranking below three no-effect variants (see "Divergence between scoring rules" below). On the small quantitative Stratum 1 (n = 8), all point estimates matched the pre-specified negative direction (Rule A ρ = −0.66; Rule B ρ = −0.29; SpliceAI ρ = −0.32; Pangolin ρ = −0.29). All three predictors comfortably exceeded the pre-registered distance-to-canonical-splice-site floor (AUC 0.683) under either scoring rule.

### Interpretation

The pre-registered inferential claim on this benchmark is non-inferiority to two established predictors, and that claim is robust across both AlphaGenome scoring rules. The superiority test is not met at n = 52; the point-estimate advantages of +0.018 (Rule A) and +0.010 (Rule B) sit inside the SE band of the paired AUC difference when all scores are near the ceiling. The two scoring rules therefore agree on the paper's inferential core: AlphaGenome performs at least as well as SpliceAI and Pangolin at classifying canonical CFTR splice variants on this benchmark, but not measurably better under the pre-registered test at this sample size.

Rank behaviour on Stratum 2 does differ between the rules. Rule A achieved perfect rank separation: every splice-affecting variant received a higher composite score than every no-effect variant. Rule B did not: one splice-affecting variant (Joynt_E02) received a lower per-junction score than three no-effect variants (Joynt_E03/E04/E05). We report both without collapsing them. Perfect rank separation is a Rule A finding; it is not a general property of AlphaGenome on this benchmark. The distinction between "a specific scoring rule achieves perfect rank separation on this benchmark" and "AlphaGenome discriminates perfectly on this benchmark" is important, and only the former is supported by these results.

The result also has to be read against a well-recognised ceiling effect: CFTR splice variants in Joynt et al. and comparable curated series are heavily enriched for canonical splice-site disruptions, which every modern sequence-based predictor identifies well. Under these conditions, discrimination benchmarks tend to saturate, and a benchmark's ability to separate top-tier predictors from each other is inherently limited. Sample size (n = 52 on Stratum 2, n = 8 on Stratum 1) was not selected by a formal power calculation; it is the full set of published quantitatively-annotated CFTR variants that met our harvest criteria and could be normalised to GRCh38 at the pre-registration date.

Stratum 1 tells a different and complementary story. Discrimination is not the question — every variant in this stratum is known to alter splicing to some degree. The question is quantitative fidelity: does a higher predictor score track a lower fraction of normally spliced transcript? Under Rule A, AlphaGenome's Spearman ρ of −0.66 is roughly twice the magnitude of either comparator (−0.32 and −0.29) in the pre-specified negative direction. Under Rule B, ρ is −0.29, comparable to the comparators. The eight-variant sample size gives wide bootstrap 95% confidence intervals that include zero under all four measurements, and the rule-specific point estimates are not statistically distinguishable from one another at this sample size. We report this as a descriptive finding per the frozen protocol and do not draw an inferential conclusion.

### Divergence between scoring rules and what it says about per-junction scoring in a restricted cell context

The two AlphaGenome scoring rules disagree on this benchmark in a systematic and mechanistically interpretable way. The Stratum 2 AUC gap between the rules (1.000 versus 0.9922) is arithmetically explained by exactly three wrong pos–neg pairs out of 43 × 9 = 387: one splice-affecting variant (Joynt_E02, Rule B score 0.221) ranked below three no-effect variants (Joynt_E03/E04/E05, Rule B scores 0.346, 0.427, 0.502). Those three no-effect variants share a distinguishing feature: the closest canonical junction inside the pre-declared cell-context ontology set (`{UBERON:0002048, UBERON:0002185, CL:0002145, CL:1000271, CL:0002632}`) lies 883, 884, and 2 892 base pairs from the variant respectively, whereas Joynt_E02's closest in-context junction is 215 bp away. When the SpliceJunctionScorer scores a distant junction that has no direct relationship to the query variant, the resulting per-junction magnitude is not a reliable measure of the variant's true splicing impact.

Leave-one-out analysis on Stratum 1 identifies the same pattern. Two variants drive the Rule A versus Rule B difference in ρ. Dropping the wild-type control S1_06_Joynt_WT (100% normally spliced, chosen-junction distance 2 892 bp, Rule B score 0.502) reduces Rule B's ρ magnitude from −0.29 to −0.04; the wild-type control receives a spuriously non-zero per-junction score because its "closest available in-context junction" is very far away. Dropping S1_01_Masvidal (40% normally spliced, chosen-junction distance 164 bp, Rule B score 6.28) raises Rule B's ρ magnitude from −0.29 to −0.61, matching Rule A. Rule A launders both problems out because it takes the maximum absolute score across all tracks and all genes in the 1 Mb interval; a spurious high in one distant junction does not swamp the true signal in the near-variant region.

The unified pattern is that a per-junction score restricted to a pre-declared cell-context ontology set is unreliable when no annotated in-context junction lies near the variant. This is not a defect in the pre-registered rule as such; the pre-registered rule was written before we could inspect the actual AlphaGenome SDK output structure and was intended to enforce the cell-context specificity that the AlphaGenome model provides. It is, however, a substantive methodological finding on this dataset: for benchmarks in which the pre-declared ontology set does not annotate a canonical junction near every variant, a cross-tissue composite score of the kind AlphaGenome's documentation recommends is the more reliable scoring rule. We report this transparently rather than choosing between the rules retrospectively.

### Comparison with the AlphaGenome preprint

AlphaGenome's original evaluation on ClinVar splicing variants ([Avsec et al., 2025, bioRxiv](https://www.biorxiv.org/content/10.1101/2025.06.25.661532v2)) reported stronger separation from SpliceAI and Pangolin than this study shows on CFTR. Two features of the CFTR benchmark likely contribute to the narrower gap observed here: the enrichment for canonical splice-site variants (which are already near-solved by all three predictors), and the modest sample size (n = 52 in Stratum 2, n = 8 in Stratum 1). This does not contradict the AlphaGenome authors' findings; it constrains where those findings apply.

### Limitations

Five limitations bound the interpretation of these results.

- **Sample size.** Stratum 2 (n = 52) sits at the discrimination ceiling for the class of variants included, limiting the paired DeLong test's power to detect small AUC gaps. Stratum 1 (n = 8) is descriptive by design; no formal inferential claim is made. Sample size was set by the available published quantitative literature at the pre-registration date, not by a prospective power calculation.
- **Variant distribution.** The frozen benchmark inherits the composition of published CFTR splice-variant series, which favour canonical splice-site variants over deep-intronic and exonic-splicing-element variants. Predictor rank-order on this benchmark may not generalise to the deep-intronic tail.
- **Assay heterogeneity.** Stratum 1's quantitative measurements combine minigene and endogenous-RNA assays from different laboratories with different reference cell lines and quantitation methods. Direct comparability across the eight variants is therefore imperfect and inflates the residual variance of any rank correlation.
- **Predictor version drift.** AlphaGenome, SpliceAI, and Pangolin are all versioned software. Scores in this benchmark were produced against pinned versions (AlphaGenome 0.9.0, SpliceAI 1.3.1, Pangolin at a specified commit); subsequent releases may yield different results. The repository is structured so that a version bump re-runs the full analysis under Amendment #2's two-rule scoring.
- **Scoring-rule sensitivity.** As documented in the Results and the "Divergence" subsection above, the choice of AlphaGenome scoring rule materially affects rank-based statistics on this benchmark. The pre-registered per-junction cell-context-restricted rule (Rule B) is not reliable on variants whose closest in-context junction lies more than a few hundred base pairs away. Reports of AlphaGenome discrimination on similarly-composed benchmarks should specify the scoring rule used.

### Compliance and scope

This study reports summary statistics computed from AlphaGenome output (ROC-AUC, DeLong confidence intervals, Spearman ρ) under two pre-declared scoring rules, and does not fit any learned function of AlphaGenome output. Publication is within the academic non-commercial research scope of the AlphaGenome Terms of Use, with the conspicuous notice, non-clinical framing, and disclaimer supplied as specified in the terms documentation accompanying the repository.

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
- **On Stratum 1 (quantitative residual-splicing variants):** the rule-sensitive AlphaGenome quantitative signal (Rule A ρ = −0.66, Rule B ρ = −0.29) is descriptive at n = 8 and does not on its own support a change in clinical practice. If confirmed in larger datasets under a robust scoring rule, quantitative splice scoring might serve as a research triage input for questions such as: which non-canonical CFTR splice variants preserve enough correctly spliced transcript for potentiator response to be plausible; which variants should be prioritised for laboratory splicing assays; which variants might be candidates for splice-modulating therapies (antisense oligonucleotides, small-molecule splice modulators). None of these applications is supported by the present sample size.

### What this study does not support

- Any change to current CFTR variant classification per ACMG/AMP or ClinGen SVI splicing subgroup rules.
- Any change to current CF modulator eligibility.
- Use of AlphaGenome, SpliceAI, or Pangolin scores as a standalone basis for clinical decisions.
- Any inference about splice fidelity for non-CFTR splice variants beyond the extent of validated cross-gene generalisation shown by each predictor's authors.
- Use of the AlphaGenome per-junction cell-context-restricted score (Rule B) on variants whose closest canonical junction inside the chosen ontology set lies more than a few hundred base pairs away, without an explicit assessment of the chosen-junction distance for each such variant.

### Where this study might inform CF research programmes

- **Variant curation.** As a supporting-strength research input alongside laboratory splicing assay evidence, per the ClinGen SVI splicing subgroup's framework for using computational splice predictions.
- **Assay prioritisation.** When resources for minigene or endogenous-transcript splicing assays are limited, quantitative splice-score rank order (from a scoring rule robust to variant-junction distance) could inform which variants are studied first.
- **Trial design.** For emerging splice-modulating therapies, quantitative predictor scores may inform the selection of candidate variants for early-phase clinical evaluation, complementing (not replacing) direct functional characterisation.

---

## Next steps

### Data expansion

The single most impactful improvement to this benchmark would be a substantially larger Stratum 1. As of September 2026, no public database (MaveDB, MPSA benchmarks, or targeted CFTR minigene collections) provides more than a few dozen quantitative CFTR splicing measurements, and the largest current MPSA studies (POU1F1, RON, FAS, WT1 by Smith and Kitzman 2023; COMPASS by Koplik et al. 2025) do not focus on CFTR. Three tractable expansion routes exist.

- **Aggregated minigene meta-benchmark.** Systematic literature harvest of published CFTR minigene and endogenous-RNA splicing assays with quantitative readouts (Sterrantino et al., Bergougnoux et al., Igreja et al., and comparable sources) into a per-source-stratified quantitative benchmark. Realistic yield: 30–60 quantitative variants. Would be published as a versioned amendment (`docs/protocol/AMENDMENTS.md` §A3) with pre-specified per-source random effects.
- **COMPASS re-analysis.** The Koplik et al. 2025 dataset covers 87 546 variants across more than 1 700 genes. CFTR-specific coverage is not published in the abstract; a targeted analysis of the COMPASS supplementary tables would either yield a substantially larger Stratum 1 or clarify the current data ceiling.
- **De novo CFTR MPSA.** A saturation MPSA focused on the CFTR gene, ideally covering canonical splice sites, deep-intronic elements, and exonic splicing elements together, would allow first-principles benchmarking without the assay-heterogeneity limitations of aggregation. This would be a multi-lab collaboration, not a single-doctor project.

### Methodological extensions

- **Scoring-rule characterisation.** A follow-up methodological analysis on a larger and more junction-dense benchmark could characterise the operating range of the per-junction cell-context-restricted scoring rule as a function of chosen-junction distance, and identify the distance beyond which the rule is not reliable. The present benchmark's chosen-junction distances range from 0 to 4 668 bp; a larger benchmark would allow this relationship to be estimated directly rather than described qualitatively.
- **Regulatory-region variants.** Extending the frozen benchmark to include deep-intronic and exonic-splicing-enhancer variants, where all three predictors are known to be less accurate than at canonical splice sites, would provide a more discriminating test of quantitative fidelity than the current variant distribution allows.
- **Cross-predictor calibration.** Formal calibration analysis (Brier score, calibration curves, isotonic recalibration) on a larger dataset would allow rank-preserving score comparisons that this study's small quantitative sample size does not.
- **Predictor version tracking.** A live continuous-integration workflow that re-runs the benchmark against each new AlphaGenome, SpliceAI, and Pangolin release would let the field track predictor drift on a stable reference set. The current repository is already structured to support this and to run each new AlphaGenome release under both scoring rules.

### Publication timeline

- Preprint (biorxiv or medrxiv) deposit as soon as the manuscript is complete.
- Peer-reviewed submission to a splice-biology, computational-biology, or CF-focused journal in parallel with any Stratum 1 expansion work.
- A follow-up paper covering the aggregated minigene meta-benchmark and the scoring-rule characterisation, if the harvest yields sufficient variants, or reporting the negative result of that harvest.

### Open science commitments

- All variant lists, predictor outputs under both scoring rules, analysis code, figures, and manuscript drafts are hash-locked and version-controlled in the accompanying public repository.
- Pre-registration remains fixed at [osf.io/5hvgf](https://osf.io/5hvgf/) (DOI [10.17605/OSF.IO/6PGX8](https://doi.org/10.17605/OSF.IO/6PGX8)); the two amendments filed on 2026-09-20 are recorded in `docs/protocol/AMENDMENTS.md` with old and new hashes and cross-linked OSF amendment identifiers.
- Predictor outputs are shared under the terms of each predictor, with AlphaGenome outputs subject to the AlphaGenome Output Terms notice recorded in the repository.
