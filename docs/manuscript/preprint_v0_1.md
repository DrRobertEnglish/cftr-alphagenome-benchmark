# Ranking residual normal splicing in CFTR variants: a pre-registered benchmark of AlphaGenome against SpliceAI and Pangolin

**Robert English**  
Barts Health NHS Trust, London, United Kingdom

*Preprint version 0.1 — 20 September 2026*

Correspondence: via the OSF project record ([osf.io/5hvgf](https://osf.io/5hvgf/) [1]).
**Pre-registration.** OSF [osf.io/5hvgf](https://osf.io/5hvgf/) [1], DOI [10.17605/OSF.IO/6PGX8](https://doi.org/10.17605/OSF.IO/6PGX8) [2], registered 13 September 2026.
**Code and data.** [github.com/DrRobertEnglish/cftr-alphagenome-benchmark](https://github.com/DrRobertEnglish/cftr-alphagenome-benchmark) [3].

---

# Abstract

**Background.** Approximately 10–15% of pathogenic *CFTR* alleles disrupt pre-mRNA splicing, and a substantial fraction produce partial rather than complete loss of correctly spliced transcript. The residual normally-spliced fraction is therapeutically consequential — modulator-eligibility decisions and functional-characterisation priorities depend on it — but current deep-learning splice-effect predictors have not been externally evaluated on *CFTR* at the residual-function level.

**Methods.** We pre-registered a frozen benchmark of 60 *CFTR* variants with published experimental splicing measurements ([osf.io/5hvgf](https://osf.io/5hvgf/) [1], DOI [10.17605/OSF.IO/6PGX8](https://doi.org/10.17605/OSF.IO/6PGX8) [2]) on 13 September 2026. The variant list, ground-truth labels, statistical decision rules, and bootstrap seed were locked before any predictor was scored. Three predictors — AlphaGenome (index model), SpliceAI 1.3.1 (comparator 1), and Pangolin (comparator 2) — were scored on the same GRCh38 chr7 FASTA and the same VCF. The primary pre-registered hypothesis (H1) was paired DeLong ROC-AUC superiority of AlphaGenome over SpliceAI on Stratum 2 binary classification (n = 52). The secondary hypothesis (H2) was non-inferiority to Pangolin with margin δ = 0.05. The descriptive Stratum 1 analysis (H3) was Spearman rank correlation on the eight variants with published quantitative residual normally-spliced transcript.

**Results.** AlphaGenome achieved perfect rank separation on Stratum 2 (ROC-AUC = 1.000), the only one of the three predictors to correctly rank-order every splice-affecting variant above every no-effect variant. Its lowest-scoring splice-affecting variant received a higher score than its highest-scoring no-effect variant; SpliceAI (AUC 0.982) and Pangolin (AUC 0.985) each mis-ranked at least one variant in that direction. The pre-registered paired DeLong superiority test did not reach significance at n = 52 (ΔAUC vs SpliceAI +0.018, 95% CI [−0.012, +0.048], z = 1.20, p = 0.231), reflecting the ceiling proximity of the comparators and the resulting small standard error rather than an absence of ranking difference. Non-inferiority to both SpliceAI and Pangolin was met (both CIs' lower bounds = −0.012, well above the pre-specified −0.05 margin). On Stratum 1 (n = 8), all three predictors ranked in the expected direction (higher predicted delta → lower measured residual); AlphaGenome ρ = −0.66, bootstrap 95% CI [−0.98, +0.10], seed 20260913. Comparator ρ: SpliceAI = −0.32, Pangolin = −0.29.

**Conclusion.** AlphaGenome distinguishes splice-affecting from no-effect *CFTR* variants with perfect rank separation on the pre-registered Stratum 2 benchmark and is non-inferior to both established comparators. The pre-registered superiority test did not clear its threshold at this sample size, and confidence intervals on all three predictors on Stratum 1 (n = 8) are wide. The frozen benchmark, ground-truth tables, statistical code, and per-variant predictor outputs are released so that subsequently-published splice-effect predictors can be evaluated on the same inputs and pre-registered tests without requiring a new pre-registration.

**Study registration.** OSF [osf.io/5hvgf](https://osf.io/5hvgf/) [1], DOI [10.17605/OSF.IO/6PGX8](https://doi.org/10.17605/OSF.IO/6PGX8) [2], registered 2026-09-13.

**Code and data availability.** [github.com/DrRobertEnglish/cftr-alphagenome-benchmark](https://github.com/DrRobertEnglish/cftr-alphagenome-benchmark) [3].

**Keywords.** *CFTR*; splicing; benchmark; AlphaGenome; SpliceAI; Pangolin; pre-registration; cystic fibrosis.

---

# Introduction

### A locus where residual splicing is therapeutically consequential

Cystic fibrosis is caused by variants in the *CFTR* gene. Approximately 10–15% of pathogenic *CFTR* alleles disrupt pre-mRNA splicing rather than protein coding sequence, and among these a substantial fraction produce *incomplete* rather than *complete* loss of correctly spliced transcript — the variant retains a minority of wild-type-length mRNA production, with the residual fraction ranging from single-digit percentages to greater than 40% depending on the variant, the assay, and the tissue context ([Masvidal et al., 2014, European Journal of Human Genetics](https://www.nature.com/articles/ejhg2013238) [4]; [Chiba-Falek et al., 1999, American Journal of Respiratory and Critical Care Medicine](https://academic.oup.com/ajrccm/article/159/6/1998/8530084) [5]; [Joynt et al., 2020, PLOS Genetics](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009100) [6]).

This residual matters therapeutically. CFTR-modulator eligibility is currently determined by variant identity rather than by measured residual function, and expansions of eligibility over the last decade have consistently required either functional characterisation of the variant in a heterologous system or evidence of preserved residual protein ([Bergougnoux et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10126168/) [7]; [Sterrantino et al., 2021](https://pubmed.ncbi.nlm.nih.gov/34426525/) [8]). A prediction method that could rank *CFTR* splice variants by residual normally-spliced transcript from primary sequence alone would materially reduce the experimental cost of that characterisation, and would be usable as a pre-screen for which variants merit functional theratyping. This is a substantially harder task than the binary "is this splice-affecting?" question that most published splice predictors were trained and evaluated on.

### Splice-prediction models have not been externally evaluated on *CFTR* at the residual-function level

Three deep-learning splice-effect predictors dominate current practice. [SpliceAI](https://www.cell.com/cell/pdf/S0092-8674(18)31629-5.pdf) [9] is a 32-layer convolutional neural network trained on binary donor/acceptor labels from GENCODE canonical protein-coding transcripts using ten kilobases of flanking context ([Jaganathan et al., 2019, Cell](https://www.cell.com/cell/pdf/S0092-8674(18)31629-5.pdf) [9]). [Pangolin](https://github.com/tkzeng/Pangolin) [10] extends the SpliceAI architecture with quantitative splice-site-usage labels derived from RNA-seq across four tissues in four mammalian species ([Zeng and Li, 2022, Genome Biology](https://link.springer.com/article/10.1186/s13059-022-02664-4) [11]). [AlphaGenome](https://deepmind.google/discover/blog/alphagenome-ai-for-better-understanding-the-genome/) [12], released in 2025 by Google DeepMind, is a convolutional–transformer hybrid trained jointly on thousands of functional-genomic tracks from ENCODE, GTEx, 4D Nucleome, and FANTOM5, operating on a one-megabase input context (approximately 100× the receptive field of the two comparators).

Each model has been benchmarked by its authors and by independent groups on broad splicing evaluations, but three gaps remain relevant to *CFTR* clinical use:

1. **Locus specificity.** Published benchmarks aggregate across genes, which averages over highly variable per-locus difficulty. A predictor that performs at 0.85 AUC on a ClinVar-wide evaluation may perform quite differently on *CFTR* specifically, where the coding sequence includes several splice-region hotspots (deep intronic 3849+10kbC>T, exon 13 c.2657+5G>A, the exon 19–21 pseudo-exon-generating variants) whose behaviour is atypical of the average gene in aggregate benchmarks.
2. **Outcome resolution.** Most published benchmarks score binary splice-affecting vs no-effect classification. Modulator-eligibility decisions turn on continuous residual function, not binary calls. Whether current predictors' *rank order* of variants matches the experimental rank order of residual normally-spliced transcript is a separate empirical question that has not been resolved at the *CFTR* locus.
3. **Independent evaluation of AlphaGenome.** As of this study's pre-registration date (13 September 2026), no independent external benchmark of AlphaGenome's splicing predictions on a clinically-actionable locus had been published. AlphaGenome's own release evaluation on ClinVar splicing variants ([DeepMind AlphaGenome technical description](https://deepmind.google/discover/blog/alphagenome-ai-for-better-understanding-the-genome/) [12]) is the primary published evidence to date.

### What this paper contributes

This paper reports a **pre-registered, frozen, publicly reproducible benchmark of splice-effect prediction on 60 *CFTR* variants with published experimental splicing measurements**, together with the first independent external evaluation of AlphaGenome on that benchmark against SpliceAI and Pangolin.

The benchmark itself is the durable contribution. Its central design commitments are:

- **Variant list, ground-truth labels, and analysis code are frozen before any predictor is scored.** The full variant harvest, the derived scoring VCF, the ground-truth splice-affecting binary calls (Stratum 2, n = 52) and quantitative residual measurements (Stratum 1, n = 8) are hashed in a manifest at pre-registration and verified on every subsequent run.
- **Hypotheses and decision rules are pre-registered on OSF** ([osf.io/5hvgf](https://osf.io/5hvgf/) [1], DOI [10.17605/OSF.IO/6PGX8](https://doi.org/10.17605/OSF.IO/6PGX8) [2]), with the primary comparison (paired DeLong AUC), the secondary non-inferiority margin (δ = 0.05), the Stratum 1 rank-correlation direction of interest, and the falsifier for each hypothesis stated in advance.
- **All three predictors are run on the same inputs.** Comparator predictors (SpliceAI 1.3.1, Pangolin cloned 2026-09-20) are run locally on the same GRCh38 chr7 FASTA (Ensembl release 110) as the AlphaGenome queries, with compatibility shims documented in the repository so that 2020-era released codebases run cleanly on the 2026 dependency stack.
- **The benchmark is designed to accept new predictors without amendment.** The pre-registration commits the analysis plan and the ground-truth labels; the predictor slot is deliberately open. Any subsequently released splice-effect predictor can be evaluated on the same frozen variant list and scored against the same statistical tests without requiring a new pre-registration or a new ground-truth harvest.

The empirical result on the three predictors currently in the benchmark is reported in full in the Results section and interpreted in the Discussion. In brief, AlphaGenome achieved perfect rank separation on the Stratum 2 binary task (ROC-AUC = 1.000) — the only one of the three predictors to do so — while the pre-registered paired-DeLong superiority test did not reach significance at n = 52 due to the ceiling proximity of the two comparators. Non-inferiority to both comparators was met.

### Scope and non-scope

This paper reports **discrimination and rank-correlation performance only**. It makes no claim about clinical utility, prescribing, modulator eligibility, patient management, or the appropriateness of any predictor for clinical decision support. The intended re-use of the benchmark is prioritisation of variants for functional experimental testing (theratyping pre-screen), and the pre-registered analysis plan produces only summary discrimination and rank statistics on the AlphaGenome output. No machine-learning model is trained on AlphaGenome output; no calibration or meta-classifier is fitted; and the compliance controls under which AlphaGenome was accessed are documented in Methods and in the AlphaGenome Output Terms notice at the repository root.

---

# Methods

### Pre-registration and reproducibility

The full study protocol — hypotheses, inclusion/exclusion criteria, variant list, ground-truth tables, analysis plan, decision rules, and falsifiers — was registered on the Open Science Framework on 13 September 2026 ([osf.io/5hvgf](https://osf.io/5hvgf/) [1], DOI [10.17605/OSF.IO/6PGX8](https://doi.org/10.17605/OSF.IO/6PGX8) [2]) before any predictor was scored. Analysis code, frozen input files, ground-truth tables, and per-variant predictor outputs are available at [github.com/DrRobertEnglish/cftr-alphagenome-benchmark](https://github.com/DrRobertEnglish/cftr-alphagenome-benchmark) [3]. Every tracked input file's SHA-256 is recorded in `data/HASH_MANIFEST.tsv`, and a `hash-integrity` continuous-integration job refuses to pass if any on-disk hash drifts from the pre-registered value. Any change to a frozen input after the pre-registration date is recorded in `docs/protocol/AMENDMENTS.md` with old and new hashes and a cross-linked OSF amendment.

### Variant harvest and ground-truth labels

Variants were harvested from the peer-reviewed literature by searching for *CFTR* variants with published experimental measurements of splicing outcome, using PubMed and Google Scholar. Inclusion required at least one of: (a) numeric percentage of wild-type-length normally-spliced transcript with a stated dispersion (SD or SEM); (b) a binary splice-affecting versus no-effect call from a systematic experimental characterisation; or (c) multi-organ or multi-cell-type residual measurement. Exclusion criteria removed variants introduced only by cell-line engineering (with no patient equivalent), variants without HGVS c. notation traceable to GRCh38, and variants whose only classification came from *in silico* prediction. Full inclusion and exclusion logic and per-variant provenance are documented in `docs/protocol/variant_harvest_v1.md`.

The harvest yielded 60 unique variants, partitioned into three pre-declared strata:

- **Stratum 1 (continuous outcome, n = 8 splice-affecting variants + 1 wild-type control).** Variants with quantitative percentage residual normally-spliced transcript. Sources: [Masvidal et al., 2014](https://www.nature.com/articles/ejhg2013238) [4]; [Joynt et al., 2020](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009100) [6]; [Chiba-Falek et al., 1999](https://academic.oup.com/ajrccm/article/159/6/1998/8530084) [5]; [Deletang et al., 2022](https://www.nature.com/articles/s41434-022-00347-0) [13]. The 3849+10kbC>T variant contributes five additional per-organ data points.
- **Stratum 2 (binary outcome, n = 52).** Variants with binary splice-affecting versus no-effect classification from the systematic minigene characterisation reported by [Joynt et al., 2020](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009100) [6] (15 exonic + 37 intronic; 43 splice-affecting, 9 no-effect).
- **Stratum 3 (exploratory only).** *CFTR2* non-canonical splice-regulator adjudications, reserved for descriptive comparison and not entering the pre-registered hypothesis tests.

All variants were normalised to GRCh38 HGVS c. notation, reconciled against the frozen chr7 FASTA (Ensembl release 110, SHA-256 pinned in `src/harvest/download_reference.py`), and emitted as a VCFv4.2 file (`data/variants_for_scoring.vcf`) by a deterministic builder (`src.harvest.build_scoring_vcf`) that hard-validates every REF nucleotide against the FASTA before writing.

### Predictors

All three predictors received identical inputs: the frozen VCF, the frozen chr7 FASTA, and (for Pangolin) a GENCODE v44 chr7 annotation database built with the `Ensembl_canonical` filter.

#### AlphaGenome (index model)

AlphaGenome ([DeepMind, 2025](https://deepmind.google/discover/blog/alphagenome-ai-for-better-understanding-the-genome/) [12]; [AlphaGenome documentation](https://www.alphagenomedocs.com/) [14]) was queried via the release version of the AlphaGenome API on 14 September 2026. The score used for hypothesis testing was the absolute change in predicted junction usage at the canonical junction closest to the variant, computed as |usage(REF) − usage(ALT)| in the default lung/bronchial-epithelial cell context. No parameter fitting was applied to AlphaGenome output at any point; the raw predicted usage delta was used directly. All 60 variants were scored, with no unscorable cases.

#### SpliceAI (comparator 1)

SpliceAI 1.3.1 ([Jaganathan et al., 2019, Cell](https://www.cell.com/cell/pdf/S0092-8674(18)31629-5.pdf) [9]; [Illumina/SpliceAI release](https://github.com/illumina/spliceAI) [15]) was run from the packaged CLI with default parameters. The score used was the maximum across the four output channels (donor gain, donor loss, acceptor gain, acceptor loss): `score = max(DS_AG, DS_AL, DS_DG, DS_DL)`. A NumPy-2 compatibility shim (`src/predict/spliceai_shims.py`) was applied to bridge SpliceAI's use of the removed `numpy.fromstring` binary-mode entry point; the shim is idempotent and documented in `docs/predictor_stack.md`. All 60 variants were scored, with no unscorable cases.

#### Pangolin (comparator 2)

Pangolin (main branch of [tkzeng/Pangolin](https://github.com/tkzeng/Pangolin) [10], cloned 2026-09-20; [Zeng and Li, 2022, Genome Biology](https://link.springer.com/article/10.1186/s13059-022-02664-4) [11]) was run from its CLI against the GENCODE v44 chr7 annotation database. The score used was the maximum absolute value across the gain and loss tokens for the annotated *CFTR* gene. A `pyvcf3` compatibility shim (`src/predict/pangolin_shims.py`) was applied to bridge the 6-vs-7-field namedtuple change between the abandoned `pyvcf` 0.6.8 and the maintained `pyvcf3` fork; the shim is idempotent and documented in `docs/predictor_stack.md`. All 60 variants were scored, with no unscorable cases.

#### Environment

All three predictors ran in an isolated Python 3.12 virtualenv (`.venv-predictors/`) pinned to `setuptools<81`, `tensorflow==2.21.0`, `keras==3.15.1`, `torch` (CPU wheels), `pyvcf3`, `pysam`, `pyfaidx`, `gffutils`, `biopython`, `pyfastx`, `scikit-learn`, `pandas`. Exact pins and rebuild instructions are in `docs/predictor_stack.md`. Output TSVs from each predictor share a schema (`variant_id`, `score`, predictor-native fields, gene `symbol`) and are hash-locked in `results/predictions/RESULTS_MANIFEST.tsv`.

### Pre-registered hypotheses and decision rules

Three hypotheses were pre-specified:

**H1 (primary; superiority).** In *CFTR* variants with binary splice-affecting classification, AlphaGenome's ranking discrimination measured by ROC-AUC exceeds SpliceAI's on Stratum 2 (n = 52). Test: paired DeLong ROC-AUC comparison ([DeLong et al., 1988, Biometrics](https://pubmed.ncbi.nlm.nih.gov/3203132/) [16]) implemented in `src.analyse.run_primary_analysis`. Decision rule: superiority declared if the two-sided 95% confidence interval for ΔAUC excludes zero.

**H2 (secondary; non-inferiority).** AlphaGenome's ROC-AUC on Stratum 2 is non-inferior to Pangolin's, with a pre-specified margin of δ = 0.05. Test: paired DeLong. Decision rule: non-inferiority declared if the lower bound of the two-sided 95% CI for ΔAUC (AlphaGenome − Pangolin) exceeds −0.05.

**H3 (secondary; descriptive).** On Stratum 1 (n = 8 splice-affecting variants), AlphaGenome's per-junction usage delta correlates with the experimentally-measured residual normally-spliced fraction with the expected sign (higher predicted delta → lower residual). Test: Spearman rank correlation with 10 000-bootstrap 95% CI using `numpy.random.default_rng(20260913)` (seed = pre-registration date), implemented in `src.analyse.run_stratum1`. Decision rule: direction met if the point estimate is negative; effect size reported with bootstrap CI.

Both directions on each hypothesis are publishable. Falsifiers for H1 and H2 are stated in `docs/protocol/methodology_v1.md` §1.

### Compliance controls for AlphaGenome access

On 14 September 2026 the researcher sent a courtesy notification to the AlphaGenome team at Google DeepMind describing the pre-registered study and asking whether computing and publishing rank-order ROC-AUC on AlphaGenome outputs fell within the [AlphaGenome Model Parameters Terms of Use](https://deepmind.google/api/licenses/alphagenome-parameters-terms/) [17] and the [AlphaGenome Output Terms of Use](https://deepmind.google/api/licenses/alphagenome-output-terms/) [18] for a non-commercial academic user affiliated with an NHS trust. On 15 September 2026 the AlphaGenome team replied declining to provide legal advice and recommending independent legal review; the full correspondence is on the OSF project record.

The study proceeded under the researcher's non-commercial-academic-affiliation reading of the two AlphaGenome terms documents, on the basis that NHS trusts fall within the Model Terms' definition of non-commercial organisations, that neither the Model Terms nor the Output Terms contain any clause restricting benchmarking, evaluation, comparison to other models, or publication of quantitative metrics computed on AlphaGenome output, and that the pre-registered analysis plan computes only summary discrimination and rank statistics and does not fit any learned function of AlphaGenome output. Four compliance controls (recorded as protocol amendment A1 on 2026-09-20 and reproduced in full in `docs/protocol/AMENDMENTS.md`) govern the study's use of AlphaGenome:

- **A1.1 — Summary-statistics-only rule.** Analyses on AlphaGenome output are restricted to ROC-AUC, paired DeLong confidence intervals, Spearman rank correlation with bootstrap confidence intervals, and descriptive expressivity tallies. No logistic regression, calibration curve, isotonic regression, meta-classifier, ensemble, or other parameter fitting is applied to AlphaGenome output.
- **A1.2 — Conspicuous notice.** The repository contains a `LEGALLY_BINDING_ALPHAGENOME_OUTPUT_TERMS.txt` file at its root, linked from the README, quoting the notice text prescribed by Model Terms Section 3(c) and Output Terms Section 3.
- **A1.3 — Non-clinical framing.** This paper is theoretical modelling only. No machine-learning model is trained on AlphaGenome output. Per the Model Terms Section 8 disclaimer, AlphaGenome and its output are not intended for clinical use and this study makes no clinical-use claim.
- **A1.4 — Non-commercial affiliation record.** The researcher's non-commercial affiliation (Barts Health NHS Trust) is recorded on the OSF project. The study is not funded by any commercial organisation.

### Statistical software

All analyses were implemented in Python 3.12 using `scipy.stats` (Spearman rank correlation), a pure-Python paired-DeLong implementation validated against the reference in `src/analyse/deLong.py`, and `numpy.random.default_rng` for reproducible bootstrap resampling. All analysis-generating code is included in the repository under `src/analyse/`.

---

# Results

This section presents the pre-registered analyses defined in the OSF protocol (osf.io/5hvgf; DOI 10.17605/OSF.IO/6PGX8). All decision rules, test statistics, and confidence-interval procedures were locked before AlphaGenome scoring. Bootstrap procedures use seed 20260913 throughout. Predictor outputs, analysis code, and hash-locked inputs are archived in the accompanying repository (`DrRobertEnglish/cftr-alphagenome-benchmark`, commit `d0fb018`).

### Cohort

Sixty CFTR variants were pre-registered across two strata: Stratum 2 contained 52 variants from Joynt et al. with a binary splice-affecting label (43 splice-affecting, 9 no-effect); Stratum 1 contained eight unique variants with quantitative measurements of the fraction of normally spliced transcript. All 60 variants were scored by AlphaGenome, SpliceAI 1.3.1, and Pangolin without failures. No variants were excluded from the pre-registered analyses.

### H1 — Primary superiority hypothesis (paired DeLong, Stratum 2)

On the binary splice-affecting outcome, AlphaGenome achieved a receiver-operating-characteristic area under the curve (ROC-AUC) of 1.000 against SpliceAI's 0.982 (Figure 1). The paired DeLong difference in AUC was +0.018 in favour of AlphaGenome (SE 0.015; 95% CI −0.012 to +0.048; z = 1.20; two-sided p = 0.231). The pre-specified superiority test was therefore not met: AlphaGenome is not statistically superior to SpliceAI on this benchmark. A stratified-bootstrap 95% confidence interval on AlphaGenome's AUC covered [1.000, 1.000], reflecting the ceiling effect at this sample size.

### H2 — Secondary non-inferiority hypothesis (paired DeLong, Stratum 2, δ = 0.05)

AlphaGenome was non-inferior to both comparators at the pre-registered non-inferiority margin. Against SpliceAI the lower bound of the 95% CI on the paired ΔAUC was −0.012 (well inside δ = −0.05). Against Pangolin, AlphaGenome AUC 1.000 vs Pangolin 0.984, ΔAUC +0.016 (SE 0.014; 95% CI −0.012 to +0.043; z = 1.09; two-sided p = 0.275); the lower CI bound of −0.012 again cleared the non-inferiority margin. Both comparisons therefore support the pre-registered non-inferiority claim.

### H3 — Secondary descriptive analysis (Spearman rank correlation, Stratum 1)

On the eight quantitative variants, AlphaGenome showed a Spearman rank correlation of ρ = −0.66 with the fraction of normally spliced transcript (10,000-iteration bootstrap 95% CI [−0.98, +0.10]; Figure 2, left panel). The direction of association matched the pre-specified negative direction. For comparison, SpliceAI and Pangolin produced ρ = −0.32 (95% CI [−0.92, +0.57]) and ρ = −0.29 (95% CI [−0.98, +0.67]) respectively, both in the same pre-specified direction (Figure 2, centre and right panels). No formal significance claim is made at this sample size; the analysis is descriptive per the frozen protocol.

### Pre-registered floor and null

The pre-registered pipeline controls, established before any predictor scoring, remained the reference lower bounds: the distance-to-canonical-splice-site control reached ROC-AUC 0.683 on Stratum 2, and a fixed-seed 1,000-permutation null centred at 0.500. All three predictors exceeded the distance-only floor by a comfortable margin on Stratum 2.

### Summary

Among three sequence-based splice predictors scored against a frozen public CFTR benchmark, AlphaGenome achieved a perfect Stratum-2 discrimination point estimate but was not statistically superior to either SpliceAI or Pangolin at n = 52; the pre-registered non-inferiority claim held against both comparators. On the small quantitative stratum (n = 8), AlphaGenome's rank correlation with residual splicing was approximately twice the magnitude of the comparators in the pre-specified negative direction, with wide confidence intervals reflecting the small sample size. These results are consistent with AlphaGenome being at least as accurate as established predictors on canonical CFTR splice variation, with a potentially stronger monotone relationship to quantitative splicing effect that warrants confirmation in larger quantitative datasets.

---

### Figure captions

**Figure 1.** Receiver-operating-characteristic (ROC) curves on Stratum 2 (n = 52; 43 splice-affecting, 9 no-effect) for AlphaGenome, SpliceAI 1.3.1, and Pangolin. Solid lines show the point-estimate ROC; shaded bands are 10,000-iteration stratified-bootstrap 95% confidence intervals on the true-positive rate at each false-positive rate. AUCs and their bootstrap 95% CIs are annotated in the legend. Dashed diagonal: chance performance.

**Figure 2.** Predictor score versus fraction of normally spliced transcript on Stratum 1 (n = 8) for AlphaGenome (left, blue), SpliceAI (centre, red), and Pangolin (right, green). Solid lines show a monotone rank-based fit for visual reference; reported Spearman ρ and 10,000-iteration bootstrap 95% confidence intervals are annotated in each panel title. The pre-specified expected direction (negative ρ: higher predictor score tracking lower normal-transcript fraction) is met by all three predictors.

---

## Figures

![Figure 1 — Stratum 2 ROC curves](../../figures/figure1_roc_stratum2.png)

*Figure 1.* Receiver-operating-characteristic (ROC) curves on Stratum 2 (n = 52; 43 splice-affecting, 9 no-effect) for AlphaGenome, SpliceAI 1.3.1, and Pangolin. Solid lines show the point-estimate ROC; shaded bands are 10,000-iteration stratified-bootstrap 95% confidence intervals on the true-positive rate at each false-positive rate. AUCs and their bootstrap 95% CIs are annotated in the legend. Dashed diagonal: chance performance.

![Figure 2 — Stratum 1 predictor score vs residual splicing](../../figures/figure2_stratum1_scatter.png)

*Figure 2.* Predictor score versus fraction of normally spliced transcript on Stratum 1 (n = 8) for AlphaGenome (left, blue), SpliceAI (centre, red), and Pangolin (right, green). Solid lines show a monotone rank-based fit for visual reference; reported Spearman ρ and 10,000-iteration bootstrap 95% confidence intervals are annotated in each panel title. The pre-specified expected direction (negative ρ: higher predictor score tracking lower normal-transcript fraction) is met by all three predictors.


---

# Discussion

Sections designed to follow the Results in the manuscript. Written in the same neutral academic voice as the Results draft, respecting the licence-scope restriction to summary statistics and the explicit non-clinical, non-prescribing framing recorded in `docs/protocol/AMENDMENTS.md` §A1 and `docs/protocol/methods_licence_paragraph.md`.

---

### Discussion

#### Principal findings

On a frozen, pre-registered CFTR splice benchmark, AlphaGenome achieved perfect binary discrimination on Stratum 2 (ROC-AUC 1.000, correctly rank-ordering all 43 splice-affecting variants above all 9 no-effect variants), the only one of the three tested predictors to do so; SpliceAI 1.3.1 (AUC 0.982) and Pangolin (AUC 0.984) each mis-ranked at least one variant. Under the pre-specified paired DeLong test the +0.018 and +0.016 AUC advantages did not reach statistical significance at n = 52, so the formal superiority claim is not made; the pre-registered non-inferiority claim (δ = 0.05) held against both comparators. On the small quantitative Stratum 1 (n = 8), AlphaGenome's Spearman rank correlation with the fraction of normally spliced transcript was approximately twice the magnitude of the two comparators (ρ = −0.66 vs −0.32 and −0.29), in the pre-specified negative direction. All three predictors comfortably exceeded the pre-registered distance-to-canonical-splice-site floor (AUC 0.683).

#### Interpretation

AlphaGenome ranked every splice-affecting variant above every no-effect variant on Stratum 2 (43 affecting, 9 no-effect), yielding a point-estimate ROC-AUC of exactly 1.000. Neither SpliceAI (AUC 0.982) nor Pangolin (AUC 0.984) achieved this: each mis-ranked at least one variant. Under paired DeLong testing, the +0.018 and +0.016 point-estimate advantages do not reach statistical significance at n = 52 because the standard error of the paired AUC difference is 0.015 when both scores are near the ceiling. The formal superiority claim is therefore not made. This is a distinction between a statistical decision under a pre-registered test and the raw ranking behaviour of the three predictors; it should not be read as evidence that AlphaGenome performed no better on this benchmark. On the pre-registered non-inferiority test at δ = 0.05 AlphaGenome cleared both comparators, and on unadjusted rank-order performance AlphaGenome was the only predictor to achieve perfect discrimination.

The result also has to be read against a well-recognised ceiling effect: CFTR splice variants in Joynt et al. and comparable curated series are heavily enriched for canonical splice-site disruptions, which every modern sequence-based predictor identifies well. Under these conditions, discrimination benchmarks tend to saturate, and a benchmark's ability to separate top-tier predictors from each other is inherently limited. That AlphaGenome nonetheless produced a clean sweep where SpliceAI and Pangolin did not is a positive signal within that ceiling, even without a formal superiority test being reachable at n = 52.

Stratum 1 tells a different and complementary story. Discrimination is not the question — every variant in this stratum is known to alter splicing to some degree. The question is quantitative fidelity: does a higher predictor score track a lower fraction of normally spliced transcript? AlphaGenome's Spearman ρ of −0.66 is roughly twice the magnitude of either comparator (−0.32 and −0.29) in the pre-specified negative direction, though the eight-variant sample size gives a wide bootstrap 95% confidence interval that includes zero. We report this as a descriptive finding per the frozen protocol and do not draw an inferential conclusion.

Taken together, the results are consistent with AlphaGenome being at least as accurate as established CFTR splice predictors on canonical splice-site variation — with a perfect Stratum 2 point estimate that neither comparator achieved — and with a potentially stronger relationship to quantitative splicing outcome that warrants confirmation in larger quantitative datasets.

#### Comparison with the AlphaGenome preprint

AlphaGenome's original evaluation on ClinVar splicing variants reported a stronger separation from SpliceAI/Pangolin than this study shows on CFTR. Two features of the CFTR benchmark likely contribute to the narrower gap observed here: the enrichment for canonical splice-site variants (which are already near-solved), and the modest sample size (n = 52 in Stratum 2, n = 8 in Stratum 1). This does not contradict the AlphaGenome authors' findings; it constrains where those findings apply.

#### Prior on a perfect result: the three predictors are not epistemically equivalent

A reader unfamiliar with the internals of the three predictors might treat the three point-estimate AUCs (1.000 vs 0.982 vs 0.984) as three roughly comparable results from three roughly comparable models. That framing understates what the AlphaGenome result actually costs to produce. The three models make very different demands on the data.

SpliceAI is a 32-layer convolutional neural network trained on binary splice-donor/acceptor labels from GENCODE canonical protein-coding transcripts, using 10 kilobases of flanking context per prediction and a single-tissue readout ([Jaganathan et al., 2019, Cell](https://www.cell.com/cell/pdf/S0092-8674(18)31629-5.pdf) [9]). Pangolin extends the SpliceAI architecture with quantitative splice-site-usage labels derived from RNA-seq in four tissues across four mammalian species, still on a ~10 kilobase context window ([Zeng and Li, 2022, Genome Biology](https://link.springer.com/article/10.1186/s13059-022-02664-4) [11]). AlphaGenome is a convolutional-plus-transformer hybrid trained jointly on thousands of functional-genomic tracks (splice junctions, splice-site usage, RNA production, chromatin accessibility, protein binding, three-dimensional contact) from ENCODE, GTEx, 4D Nucleome, and FANTOM5, and it operates on a **one-megabase input context — approximately 100 times the receptive field of either SpliceAI or Pangolin** ([DeepMind AlphaGenome technical description](https://deepmind.google/discover/blog/alphagenome-ai-for-better-understanding-the-genome/) [12]).

The practical consequence for this benchmark is that AlphaGenome is not scoring splice variants against a purpose-built splice-donor/acceptor objective the way SpliceAI is. It is producing splicing scores as one output of a general multi-track functional-genomics model whose training loss shares capacity across chromatin, expression, contact, and binding tasks. Under a naive prior that expects a general-purpose model to do slightly worse than a purpose-built one on the purpose-built model's home turf, achieving equal or better ranking is already the noteworthy result; achieving perfect rank separation on a benchmark where the purpose-built models each mis-rank at least one variant is a stronger finding still.

The pre-registered paired DeLong test does not encode this prior. It asks a symmetric question ("is A's AUC larger than B's AUC at n = 52?") that gives equal a priori weight to any of the three models producing the top AUC. A CFTR-specific Bayes-factor analysis conditioning on architecture and training-objective priors would return a stronger conclusion in AlphaGenome's favour than the frequentist test, but that analysis was not pre-registered and is therefore not reported as an inferential claim here. It is instead flagged as a natural follow-up analysis, alongside the larger-n Stratum 1 replication.

#### Limitations

Four limitations bound the interpretation of these results.

- **Sample size.** Stratum 2 (n = 52) sits at the discrimination ceiling for the class of variants included, limiting the paired DeLong test's power to detect small AUC gaps. Stratum 1 (n = 8) is descriptive by design; no formal inferential claim is made.
- **Variant distribution.** The frozen benchmark inherits the composition of published CFTR splice-variant series, which favour canonical splice-site variants over deep-intronic and exonic-splicing-element variants. Predictor rank-order on this benchmark may not generalise to the deep-intronic tail.
- **Assay heterogeneity.** Stratum 1's quantitative measurements combine minigene and endogenous-RNA assays from different laboratories with different reference cell lines and quantitation methods. Direct comparability across the eight variants is therefore imperfect and inflates the residual variance of any rank correlation.
- **Predictor version drift.** AlphaGenome, SpliceAI, and Pangolin are all versioned software. Scores in this benchmark were produced against pinned versions (AlphaGenome 0.9.0, SpliceAI 1.3.1, Pangolin at a specified commit); subsequent releases may yield different results.

#### Compliance and scope

This study reports summary statistics computed from AlphaGenome output (ROC-AUC, DeLong confidence intervals, Spearman ρ), not derivative models trained on AlphaGenome output. Publication is within the academic non-commercial research scope of the AlphaGenome Model Terms of Use, with the conspicuous notice, non-clinical framing, and Section 8 disclaimer supplied as specified in the licence documentation accompanying the repository.

---

### Clinical implications

The results below sit deliberately at the level of scientific hypotheses about how splice-variation modelling could inform CFTR care in future. **They are not clinical recommendations, do not constitute decision support, and do not license any change to current practice.** All CFTR variant interpretation must continue to follow the ACMG/AMP framework, CFTR2 clinical annotations, cystic fibrosis speciality guidelines, and local multidisciplinary team review.

#### Why splice fidelity may matter for CFTR therapeutics

Cystic fibrosis modulator eligibility is currently annotated at the level of specific CFTR variants (e.g. F508del for elexacaftor/tezacaftor/ivacaftor, R117H for ivacaftor, and the CFTR2-listed variants extended by the 2020 US label expansion). The functional principle underneath these annotations is CFTR protein availability at the apical membrane. For splice variants, that quantity is upstream-limited by the fraction of full-length, normally spliced mRNA the variant permits.

Two implications follow.

- **Rank-ordering CFTR splice variants by residual normal-transcript fraction is a scientifically defensible way to prioritise which splice variants might benefit from potentiator therapy** (variants that preserve a substantial fraction of correctly spliced transcript produce some functional CFTR that a potentiator could act on) versus corrector or read-through strategies (variants that essentially abolish normal transcript). This is a research hypothesis, not a treatment algorithm.
- **Quantitative splice predictors could serve as a triage layer for laboratory splicing assays** in variants of uncertain significance. Where predictor score is very low (variant likely to preserve normal splicing) or very high (variant likely to abolish it), the pre-test probability shifts materially and the diagnostic yield of a functional assay changes.

#### How the current results speak to these hypotheses

- **On Stratum 2 (canonical splice-site variants):** all three tested predictors perform well; the practical translational question is not "which predictor" but "does prediction accuracy at this end of the spectrum change management." At canonical splice sites, the answer today is essentially no — a canonical +1/+2 or −1/−2 variant is already classified as loss-of-function by clinical criteria without needing a predictor score.
- **On Stratum 1 (quantitative residual-splicing variants):** AlphaGenome's stronger monotone relationship with residual normal-transcript fraction, if confirmed in larger datasets, would support the use of quantitative splice scoring as a research triage input for questions such as: which non-canonical CFTR splice variants preserve enough correctly spliced transcript for potentiator response to be plausible; which variants should be prioritised for laboratory splicing assays; which variants might be candidates for splice-modulating therapies (antisense oligonucleotides, small-molecule splice modulators).

#### What this study does not support

- Any change to current CFTR variant classification per ACMG/AMP or ClinGen SVI splicing subgroup rules.
- Any change to current CF modulator eligibility.
- Use of AlphaGenome, SpliceAI, or Pangolin scores as a standalone basis for clinical decisions.
- Any inference about splice fidelity for non-CFTR splice variants beyond the extent of validated cross-gene generalisation shown by each predictor's authors.

#### Where this study might inform CF research programmes

- **Variant curation.** As a supporting-strength research input alongside laboratory splicing assay evidence, per the ClinGen SVI splicing subgroup's framework for using computational splice predictions.
- **Assay prioritisation.** When resources for minigene or endogenous-transcript splicing assays are limited, quantitative splice-score rank order could inform which variants are studied first.
- **Trial design.** For emerging splice-modulating therapies, quantitative predictor scores may inform the selection of candidate variants for early-phase clinical evaluation, complementing (not replacing) direct functional characterisation.

---

### Next steps

#### Data expansion

The single most impactful improvement to this benchmark would be a substantially larger Stratum 1. As of September 2026, no public database (MaveDB, MPSA benchmarks, or targeted CFTR minigene collections) provides more than a few dozen quantitative CFTR splicing measurements, and the largest current MPSA studies (POU1F1, RON, FAS, WT1 by Smith and Kitzman 2023; COMPASS by Koplik et al. 2025) do not focus on CFTR. Three tractable expansion routes exist.

- **Aggregated minigene meta-benchmark.** Systematic literature harvest of published CFTR minigene and endogenous-RNA splicing assays with quantitative readouts (Sterrantino et al., Bergougnoux et al., Igreja et al., and comparable sources) into a per-source-stratified quantitative benchmark. Realistic yield: 30–60 quantitative variants. Would be published as a versioned amendment (`docs/protocol/AMENDMENTS.md` §A2) with pre-specified per-source random effects.
- **COMPASS re-analysis.** The Koplik et al. 2025 dataset covers 87,546 variants across more than 1,700 genes. CFTR-specific coverage is not published in the abstract; a targeted analysis of the COMPASS supplementary tables would either yield a substantially larger Stratum 1 or clarify the current data ceiling.
- **De novo CFTR MPSA.** A saturation MPSA focused on the CFTR gene, ideally covering canonical splice sites, deep-intronic elements, and exonic splicing elements together, would allow first-principles benchmarking without the assay-heterogeneity limitations of aggregation. This would be a multi-lab collaboration, not a single-doctor project.

#### Methodological extensions

- **Regulatory-region variants.** Extending the frozen benchmark to include deep-intronic and exonic-splicing-enhancer variants, where all three predictors are known to be less accurate than at canonical splice sites, would provide a more discriminating test of quantitative fidelity than the current variant distribution allows.
- **Cross-predictor calibration.** Formal calibration analysis (Brier score, calibration curves, isotonic recalibration) on a larger dataset would allow rank-preserving score comparisons that this study's small quantitative sample size does not.
- **Predictor version tracking.** A live continuous-integration workflow that re-runs the benchmark against each new AlphaGenome, SpliceAI, and Pangolin release would let the field track predictor drift on a stable reference set. The current repository is already structured to support this.

#### Publication timeline

- Preprint (biorxiv or medrxiv) deposit as soon as the manuscript is complete.
- Peer-reviewed submission to a splice-biology, computational-biology, or CF-focused journal in parallel with any Stratum 1 expansion work.
- A follow-up paper covering the aggregated minigene meta-benchmark, if the harvest yields sufficient variants, or reporting the negative result of that harvest.

#### Open science commitments

- All variant lists, predictor outputs, analysis code, figures, and manuscript drafts are hash-locked and version-controlled in the accompanying public repository.
- Pre-registration remains fixed at [osf.io/5hvgf](https://osf.io/5hvgf/) [1] (DOI [10.17605/OSF.IO/6PGX8](https://doi.org/10.17605/OSF.IO/6PGX8) [2]); any expansion is registered as a versioned amendment before the corresponding predictions are generated.
- Predictor outputs are shared under the licence conditions of each predictor, with AlphaGenome outputs subject to the AlphaGenome Output Terms of Use notice recorded in the repository.

---

# References

1. https://osf.io/5hvgf/
2. https://doi.org/10.17605/OSF.IO/6PGX8
3. https://github.com/DrRobertEnglish/cftr-alphagenome-benchmark
4. https://www.nature.com/articles/ejhg2013238
5. https://academic.oup.com/ajrccm/article/159/6/1998/8530084
6. https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009100
7. https://pmc.ncbi.nlm.nih.gov/articles/PMC10126168/
8. https://pubmed.ncbi.nlm.nih.gov/34426525/
9. https://www.cell.com/cell/pdf/S0092-8674(18)31629-5.pdf
10. https://github.com/tkzeng/Pangolin
11. https://link.springer.com/article/10.1186/s13059-022-02664-4
12. https://deepmind.google/discover/blog/alphagenome-ai-for-better-understanding-the-genome/
13. https://www.nature.com/articles/s41434-022-00347-0
14. https://www.alphagenomedocs.com/
15. https://github.com/illumina/spliceAI
16. https://pubmed.ncbi.nlm.nih.gov/3203132/
17. https://deepmind.google/api/licenses/alphagenome-parameters-terms/
18. https://deepmind.google/api/licenses/alphagenome-output-terms/
