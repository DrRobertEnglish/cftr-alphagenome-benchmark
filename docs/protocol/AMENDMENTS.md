# Protocol amendments

Any change to the frozen inputs (variant list, reference annotation, ground-truth tables) or to the pre-registered analysis plan after the OSF registration date (2026-09-13) must be recorded here **before** any further predictor query.

The `hash-integrity` CI job refuses to pass unless every tracked file's on-disk SHA-256 matches the recorded hash in `data/HASH_MANIFEST.tsv`. Regenerating that manifest is the mechanical part of an amendment — the substantive part is a versioned OSF amendment cross-linking to the specific change.

## Format

Each amendment gets a section below with the following fields:

- **Amendment number:** monotonically increasing integer.
- **Date:** YYYY-MM-DD (UTC).
- **OSF amendment DOI or link:** required.
- **Files changed:** repo-relative paths.
- **Old SHA-256 / New SHA-256:** for each changed file.
- **Rationale:** one paragraph.
- **Effect on hypotheses:** whether H1/H2/H3 decision rules are affected. If yes, the amendment must state whether the original hypotheses are being replaced (rare, requires justification) or supplemented (preferred; adds an exploratory hypothesis to Other Planned Analyses).

## Amendments

### Note on derived-frozen inputs added after registration

The hash manifest additionally tracks files that are **algorithmically derived** from an OSF-registered input by a deterministic, auditable transformation:

| Derived file | Derived from | Transformation |
|---|---|---|
| `data/variants_for_scoring.vcf` | `data/variant_list_frozen.tsv` | `src.harvest.build_scoring_vcf` — sort by POS ascending, hard-validate every REF nucleotide against the GRCh38 chr7 FASTA (hash pinned in `src.harvest.download_reference`), emit VCFv4.2. |

These are not protocol amendments — they add no new degrees of freedom over what is already OSF-locked. They exist so that CI can catch accidental corruption of derived scoring inputs (e.g. a re-run that produces a different VCF from the same TSV means the builder or the FASTA has drifted).

### Substantive amendments

#### 2026-09-20 (later that day) — Amendment #2: report both DeepMind-recommended composite AND pre-registered per-junction rule as parallel primary analyses

**Amendment number:** 2
**Date:** 2026-09-20
**OSF amendment DOI or link:** pending upload (draft in `docs/protocol/AMENDMENTS.md`)
**Files changed:**
- `docs/protocol/AMENDMENTS.md` (this entry)
- `scripts/colab_alphagenome_scoring_v2.ipynb` (new; adds pre-registered-rule cells alongside composite-rule cells)
- `scripts/prereg_scoring_cells.py` (new; source for the new cells with full design docstring)
- `scripts/build_updated_colab_notebook.py` (new; deterministically builds the v2 notebook)
- Following the Colab run: `results/predictions/alphagenome_prereg.tsv` (new), `data/HASH_MANIFEST.tsv` (updated), `docs/manuscript/methods_v0_2.md` (rewritten scoring section).

**Rationale.** During pre-submission critical review of preprint v0.1, we discovered that the AlphaGenome scoring rule executed in the Colab notebook did not correspond to the rule written into the OSF pre-registration.

The pre-registration specified: `AlphaGenome per-junction usage delta = |usage(REF) − usage(ALT)| at the canonical junction closest to the variant, default lung/bronchial epithelial cell context` ([OSF pre-registration v1.0](../protocol/osf_preregistration_v1.md), Predictors §1; [methodology v1.0](../protocol/methodology_v1.md) §4.1).

The Colab notebook implemented: `max(splice_sites) + max(splice_site_usage) + max(splice_junctions) / 5`, aggregated across all tissues and genes by taking the maximum absolute score across all tracks and genes ([scripts/colab_alphagenome_scoring.ipynb](../../scripts/colab_alphagenome_scoring.ipynb), cell 8).

The composite score was not chosen post-hoc: it is the DeepMind-documented recommended splicing composite specified verbatim in the [AlphaGenome FAQ](https://www.alphagenomedocs.com/faqs.html) ("How to score splicing variants") and in the [AlphaGenome splicing tutorial](https://www.alphagenomedocs.com/colabs/splicing_variant_scoring.html). Its provenance predates this study. However, the pre-registration did not specify this composite: it specified a narrower per-junction rule that only became partially implementable when the AlphaGenome SDK's actual output structure was inspected in detail (`SPLICE_SITE_USAGE` returns per-track magnitudes aggregated within a gene mask, not per-junction values; per-junction resolution requires the separate `SpliceJunctionScorer`, whose output does expose `junction_Start` and `junction_End` columns).

We therefore adopt the following resolution, chosen after the researcher and the analysis pipeline both reviewed the situation:

1. The DeepMind-recommended composite score is reported as **primary analysis (Rule A)**. It is defensible on the basis of DeepMind's own documentation and gives a reproducible, deterministic scoring procedure that any independent group can replicate.
2. The OSF-preregistered per-junction rule is operationalised as **primary analysis (Rule B)** using the SpliceJunctionScorer with a pre-declared cell-context ontology set (see `scripts/prereg_scoring_cells.py` docstring) and a distance-to-variant tie-break. The exact ontology set — `{UBERON:0002048, UBERON:0002185, CL:0002145, CL:1000271, CL:0002632}` — is frozen in this amendment and cannot be adjusted post-hoc.
3. **Both rules are reported side-by-side as parallel primary analyses.** Neither is designated "the" primary. The paper's Results section carries columns for each rule and the Discussion section addresses agreement and disagreement between them.
4. Any variant for which Rule B yields no in-context tracks (i.e. no tracks in the pre-declared ontology set have per-junction data for the variant's interval) is recorded as `NA_no_cell_context` and excluded from Rule B analyses only; it remains in Rule A analyses. This exclusion criterion is frozen in this amendment.
5. Multiplicity is addressed by pre-registering that inference from either rule (Rule A OR Rule B) is treated as one hypothesis test, i.e. the H1/H2/H3 alpha threshold applies to whichever rule the paper's inferential claim rests on for that hypothesis. When both rules give the same conclusion at nominal alpha, the finding is described as robust across scoring rules. When they differ, the disagreement is reported prominently and no unqualified inferential claim is made.

**Effect on hypotheses.** H1, H2 and H3 (as pre-registered) are executed twice: once with Rule A scores as the AlphaGenome predictor, once with Rule B scores. The decision rules themselves (paired DeLong for H1 and H2, Spearman for H3) are unchanged. This is a **supplementation** rather than a replacement in the sense of the format above: the pre-registered rule is retained (as Rule B) and an additional pre-declared rule (Rule A) is added.

**Compliance with A1.1 (summary-statistics-only rule).** The `sj / 5` weighting in Rule A is not a fitted parameter — it is a fixed weight taken from DeepMind's documentation. Neither rule involves training or fitting a model on AlphaGenome output. Both rules produce a single scalar per variant that then enters standard rank-order statistics. Rule A therefore does not tension A1.1.

**Time-order note.** This amendment is filed **after** the discovery of the deviation but **before** any pre-registered-rule scores exist. The pre-registered-rule (Rule B) TSV is generated only after this amendment is committed. Rule A scores already exist in `results/predictions/alphagenome.tsv` (SHA-256 recorded 2026-09-15). The paper's analysis code will be re-run once Rule B scores are in the repository.

---

**Amendment number:** 1
**Date:** 2026-09-20

Context. On 2026-09-14 the researcher sent a courtesy notification to the AlphaGenome team at Google DeepMind describing the pre-registered study and asking whether computing and publishing rank-order ROC-AUC on AlphaGenome outputs fell within the AlphaGenome Model Parameters Terms of Use and the AlphaGenome Output Terms of Use for a non-commercial academic user affiliated with an NHS trust. On 2026-09-15 the AlphaGenome team replied declining to provide legal advice and recommending independent legal review (see the OSF project record for the full correspondence). No further correspondence has been received.

Decision. The study proceeds under the researcher's non-commercial-academic-affiliation reading of the two AlphaGenome terms documents, on the basis that (a) NHS trusts fall within the Model Terms' definition of non-commercial organisations ("universities, non-profit organizations and research institutes, educational, journalism and government bodies"), (b) neither the Model Terms nor the Output Terms contain any clause restricting benchmarking, evaluation, comparison to other models, or publication of quantitative metrics computed on AlphaGenome output, and (c) the pre-registered analysis plan computes only summary discrimination and rank statistics and does not fit any learned function of AlphaGenome output.

Compliance controls. The following controls are recorded on this date and applied to all subsequent work in this repository, all preprint and journal submissions, and all supplementary data releases:

- **A1.1 Summary-statistics-only rule.** Analyses computed on AlphaGenome output are restricted to ROC-AUC, paired DeLong confidence intervals, Spearman rank correlation with bootstrap confidence intervals, and descriptive expressivity tallies. No logistic regression, calibration, isotonic regression, meta-classifier, ensemble, or other parameter fitting is applied to AlphaGenome output. This operationalises firewall F1.
- **A1.2 Conspicuous notice.** The repository contains a `LEGALLY_BINDING_ALPHAGENOME_OUTPUT_TERMS.txt` file at its root, linked from the README, quoting the notice text prescribed by Model Terms Section 3(c) and Output Terms Section 3 and referencing the Output Terms URL. Any supplementary file containing per-variant AlphaGenome scores or AlphaGenome-derived statistics is accompanied by a reference to this notice.
- **A1.3 Non-clinical framing.** The Methods section of the manuscript will quote the Model Terms Section 8 disclaimer verbatim, will state that the paper is theoretical modelling only, will state that no machine learning model is trained on AlphaGenome output, and will reference the Output Terms URL. The pre-registered banned-word list continues to be enforced across the full manuscript.
- **A1.4 Non-commercial affiliation record.** The researcher's non-commercial affiliation (Barts Health NHS Trust) is recorded on the OSF project. The study is not funded by any commercial organisation and produces no output for the benefit of any commercial organisation.

Degrees of freedom added. None. A1 does not change the frozen variant list, the pre-registered hypotheses, the pre-registered outcome definitions, or the pre-registered analysis plan. It records the operational compliance controls under which the pre-registered plan is executed.

Pre-registered fallback status. The fallback path in Section 9 ("Publish expressivity-only report if DeepMind licence clarification denies rank-statistic use") is not invoked because the AlphaGenome reply did not deny rank-statistic use — it declined to make any legal determination. The fallback remains available at any point up to submission if, before then, either the AlphaGenome team communicates an explicit restriction or an independent legal opinion advises against publication of rank-order metrics.


## Predictor scoring log (informational, not amendments)

### 2026-09-20 — SpliceAI 1.3.1 and Pangolin (tkzeng, main branch, cloned 2026-09-20) both scored

Both comparator predictors were run against `data/variants_for_scoring.vcf` on the frozen chr7 GRCh38 FASTA (Ensembl release 110). SpliceAI was run via the packaged CLI with a NumPy-2 compatibility shim (`src/predict/spliceai_shims.py`); Pangolin was run via a pyvcf3 compatibility shim (`src/predict/pangolin_shims.py`) against a GENCODE v44 chr7 annotation database built with `Ensembl_canonical` filter.

Outputs are recorded in `results/predictions/RESULTS_MANIFEST.tsv` with their SHA-256 hashes:

- `results/predictions/spliceai.tsv` — 60/60 variants scored, no unscorable variants.
- `results/predictions/pangolin.tsv` — 60/60 variants scored, no unscorable variants.

These are results, not inputs, so they are recorded in a separate results manifest rather than `data/HASH_MANIFEST.tsv`. The primary analysis has not yet been run against the AlphaGenome scores — those remain gated on the AlphaGenome licence decision (cron `ad481e61`, fires 2026-09-27).

A smoke run using the SpliceAI and Pangolin scores as if they were the two comparators (mechanics test only, not a pre-registered analysis) confirmed the pipeline runs end-to-end and produces sensible AUCs. Those smoke outputs live under `results/smoke/` (gitignored).
