# Methods

## Pre-registration and reproducibility

The full study protocol — hypotheses, inclusion/exclusion criteria, variant list, ground-truth tables, analysis plan, decision rules, and falsifiers — was registered on the Open Science Framework on 13 September 2026 ([osf.io/5hvgf](https://osf.io/5hvgf/), DOI [10.17605/OSF.IO/6PGX8](https://doi.org/10.17605/OSF.IO/6PGX8)) before any predictor was scored. Analysis code, frozen input files, ground-truth tables, and per-variant predictor outputs are available at [github.com/DrRobertEnglish/cftr-alphagenome-benchmark](https://github.com/DrRobertEnglish/cftr-alphagenome-benchmark). Every tracked input file's SHA-256 is recorded in `data/HASH_MANIFEST.tsv`, and a `hash-integrity` continuous-integration job refuses to pass if any on-disk hash drifts from the pre-registered value. Any change to a frozen input after the pre-registration date is recorded in `docs/protocol/AMENDMENTS.md` with old and new hashes and a cross-linked OSF amendment.

## Variant harvest and ground-truth labels

Variants were harvested from the peer-reviewed literature by searching for *CFTR* variants with published experimental measurements of splicing outcome, using PubMed and Google Scholar. Inclusion required at least one of: (a) numeric percentage of wild-type-length normally-spliced transcript with a stated dispersion (SD or SEM); (b) a binary splice-affecting versus no-effect call from a systematic experimental characterisation; or (c) multi-organ or multi-cell-type residual measurement. Exclusion criteria removed variants introduced only by cell-line engineering (with no patient equivalent), variants without HGVS c. notation traceable to GRCh38, and variants whose only classification came from *in silico* prediction. Full inclusion and exclusion logic and per-variant provenance are documented in `docs/protocol/variant_harvest_v1.md`.

The harvest yielded 60 unique variants, partitioned into three pre-declared strata:

- **Stratum 1 (continuous outcome, n = 8 splice-affecting variants + 1 wild-type control).** Variants with quantitative percentage residual normally-spliced transcript. Sources: [Masvidal et al., 2014](https://www.nature.com/articles/ejhg2013238); [Joynt et al., 2020](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009100); [Chiba-Falek et al., 1999](https://academic.oup.com/ajrccm/article/159/6/1998/8530084); [Deletang et al., 2022](https://www.nature.com/articles/s41434-022-00347-0). The 3849+10kbC>T variant contributes five additional per-organ data points.
- **Stratum 2 (binary outcome, n = 52).** Variants with binary splice-affecting versus no-effect classification from the systematic minigene characterisation reported by [Joynt et al., 2020](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009100) (15 exonic + 37 intronic; 43 splice-affecting, 9 no-effect).
- **Stratum 3 (exploratory only).** *CFTR2* non-canonical splice-regulator adjudications, reserved for descriptive comparison and not entering the pre-registered hypothesis tests.

All variants were normalised to GRCh38 HGVS c. notation, reconciled against the frozen chr7 FASTA (Ensembl release 110, SHA-256 pinned in `src/harvest/download_reference.py`), and emitted as a VCFv4.2 file (`data/variants_for_scoring.vcf`) by a deterministic builder (`src.harvest.build_scoring_vcf`) that hard-validates every REF nucleotide against the FASTA before writing.

## Predictors

All three predictors received identical inputs: the frozen VCF, the frozen chr7 FASTA, and (for Pangolin) a GENCODE v44 chr7 annotation database built with the `Ensembl_canonical` filter.

### AlphaGenome (index model)

AlphaGenome ([DeepMind, 2025](https://deepmind.google/discover/blog/alphagenome-ai-for-better-understanding-the-genome/); [AlphaGenome documentation](https://www.alphagenomedocs.com/)) was queried via the release version of the AlphaGenome API on 14 September 2026. The score used for hypothesis testing was the absolute change in predicted junction usage at the canonical junction closest to the variant, computed as |usage(REF) − usage(ALT)| in the default lung/bronchial-epithelial cell context. No parameter fitting was applied to AlphaGenome output at any point; the raw predicted usage delta was used directly. All 60 variants were scored, with no unscorable cases.

### SpliceAI (comparator 1)

SpliceAI 1.3.1 ([Jaganathan et al., 2019, Cell](https://www.cell.com/cell/pdf/S0092-8674(18)31629-5.pdf); [Illumina/SpliceAI release](https://github.com/illumina/spliceAI)) was run from the packaged CLI with default parameters. The score used was the maximum across the four output channels (donor gain, donor loss, acceptor gain, acceptor loss): `score = max(DS_AG, DS_AL, DS_DG, DS_DL)`. A NumPy-2 compatibility shim (`src/predict/spliceai_shims.py`) was applied to bridge SpliceAI's use of the removed `numpy.fromstring` binary-mode entry point; the shim is idempotent and documented in `docs/predictor_stack.md`. All 60 variants were scored, with no unscorable cases.

### Pangolin (comparator 2)

Pangolin (main branch of [tkzeng/Pangolin](https://github.com/tkzeng/Pangolin), cloned 2026-09-20; [Zeng and Li, 2022, Genome Biology](https://link.springer.com/article/10.1186/s13059-022-02664-4)) was run from its CLI against the GENCODE v44 chr7 annotation database. The score used was the maximum absolute value across the gain and loss tokens for the annotated *CFTR* gene. A `pyvcf3` compatibility shim (`src/predict/pangolin_shims.py`) was applied to bridge the 6-vs-7-field namedtuple change between the abandoned `pyvcf` 0.6.8 and the maintained `pyvcf3` fork; the shim is idempotent and documented in `docs/predictor_stack.md`. All 60 variants were scored, with no unscorable cases.

### Environment

All three predictors ran in an isolated Python 3.12 virtualenv (`.venv-predictors/`) pinned to `setuptools<81`, `tensorflow==2.21.0`, `keras==3.15.1`, `torch` (CPU wheels), `pyvcf3`, `pysam`, `pyfaidx`, `gffutils`, `biopython`, `pyfastx`, `scikit-learn`, `pandas`. Exact pins and rebuild instructions are in `docs/predictor_stack.md`. Output TSVs from each predictor share a schema (`variant_id`, `score`, predictor-native fields, gene `symbol`) and are hash-locked in `results/predictions/RESULTS_MANIFEST.tsv`.

## Pre-registered hypotheses and decision rules

Three hypotheses were pre-specified:

**H1 (primary; superiority).** In *CFTR* variants with binary splice-affecting classification, AlphaGenome's ranking discrimination measured by ROC-AUC exceeds SpliceAI's on Stratum 2 (n = 52). Test: paired DeLong ROC-AUC comparison ([DeLong et al., 1988, Biometrics](https://pubmed.ncbi.nlm.nih.gov/3203132/)) implemented in `src.analyse.run_primary_analysis`. Decision rule: superiority declared if the two-sided 95% confidence interval for ΔAUC excludes zero.

**H2 (secondary; non-inferiority).** AlphaGenome's ROC-AUC on Stratum 2 is non-inferior to Pangolin's, with a pre-specified margin of δ = 0.05. Test: paired DeLong. Decision rule: non-inferiority declared if the lower bound of the two-sided 95% CI for ΔAUC (AlphaGenome − Pangolin) exceeds −0.05.

**H3 (secondary; descriptive).** On Stratum 1 (n = 8 splice-affecting variants), AlphaGenome's per-junction usage delta correlates with the experimentally-measured residual normally-spliced fraction with the expected sign (higher predicted delta → lower residual). Test: Spearman rank correlation with 10 000-bootstrap 95% CI using `numpy.random.default_rng(20260913)` (seed = pre-registration date), implemented in `src.analyse.run_stratum1`. Decision rule: direction met if the point estimate is negative; effect size reported with bootstrap CI.

Both directions on each hypothesis are publishable. Falsifiers for H1 and H2 are stated in `docs/protocol/methodology_v1.md` §1.

## Compliance controls for AlphaGenome access

On 14 September 2026 the researcher sent a courtesy notification to the AlphaGenome team at Google DeepMind describing the pre-registered study and asking whether computing and publishing rank-order ROC-AUC on AlphaGenome outputs fell within the [AlphaGenome Model Parameters Terms of Use](https://deepmind.google/api/licenses/alphagenome-parameters-terms/) and the [AlphaGenome Output Terms of Use](https://deepmind.google/api/licenses/alphagenome-output-terms/) for a non-commercial academic user affiliated with an NHS trust. On 15 September 2026 the AlphaGenome team replied declining to provide legal advice and recommending independent legal review; the full correspondence is on the OSF project record.

The study proceeded under the researcher's non-commercial-academic-affiliation reading of the two AlphaGenome terms documents, on the basis that NHS trusts fall within the Model Terms' definition of non-commercial organisations, that neither the Model Terms nor the Output Terms contain any clause restricting benchmarking, evaluation, comparison to other models, or publication of quantitative metrics computed on AlphaGenome output, and that the pre-registered analysis plan computes only summary discrimination and rank statistics and does not fit any learned function of AlphaGenome output. Four compliance controls (recorded as protocol amendment A1 on 2026-09-20 and reproduced in full in `docs/protocol/AMENDMENTS.md`) govern the study's use of AlphaGenome:

- **A1.1 — Summary-statistics-only rule.** Analyses on AlphaGenome output are restricted to ROC-AUC, paired DeLong confidence intervals, Spearman rank correlation with bootstrap confidence intervals, and descriptive expressivity tallies. No logistic regression, calibration curve, isotonic regression, meta-classifier, ensemble, or other parameter fitting is applied to AlphaGenome output.
- **A1.2 — Conspicuous notice.** The repository contains a `LEGALLY_BINDING_ALPHAGENOME_OUTPUT_TERMS.txt` file at its root, linked from the README, quoting the notice text prescribed by Model Terms Section 3(c) and Output Terms Section 3.
- **A1.3 — Non-clinical framing.** This paper is theoretical modelling only. No machine-learning model is trained on AlphaGenome output. Per the Model Terms Section 8 disclaimer, AlphaGenome and its output are not intended for clinical use and this study makes no clinical-use claim.
- **A1.4 — Non-commercial affiliation record.** The researcher's non-commercial affiliation (Barts Health NHS Trust) is recorded on the OSF project. The study is not funded by any commercial organisation.

## Statistical software

All analyses were implemented in Python 3.12 using `scipy.stats` (Spearman rank correlation), a pure-Python paired-DeLong implementation validated against the reference in `src/analyse/deLong.py`, and `numpy.random.default_rng` for reproducible bootstrap resampling. All analysis-generating code is included in the repository under `src/analyse/`.
