# Results (draft v0.2)

This section presents the pre-registered analyses defined in the OSF protocol ([osf.io/5hvgf](https://osf.io/5hvgf/); DOI [10.17605/OSF.IO/6PGX8](https://doi.org/10.17605/OSF.IO/6PGX8)). All decision rules, test statistics, and confidence-interval procedures were locked before AlphaGenome scoring. Bootstrap procedures use seed 20260913 (the pre-registration date). Each pre-registered hypothesis is reported under both AlphaGenome scoring rules — the DeepMind-recommended merged composite (Rule A) and the OSF-preregistered per-junction score in a lung/bronchial cell context (Rule B) — as parallel primary analyses (protocol Amendment #2). Predictor outputs, analysis code, and hash-locked inputs are archived in the accompanying repository ([DrRobertEnglish/cftr-alphagenome-benchmark](https://github.com/DrRobertEnglish/cftr-alphagenome-benchmark), commit tip at time of writing `6399063`).

## Cohort

Sixty CFTR variants were pre-registered across two analytical strata: Stratum 2 contained 52 variants from Joynt et al. with a binary splice-affecting label (43 splice-affecting, 9 no-effect); Stratum 1 contained eight unique splice-affecting variants with quantitative measurements of the fraction of normally spliced transcript, plus one wild-type control. All 60 variants were scored by AlphaGenome under both scoring rules, by SpliceAI 1.3.1, and by Pangolin without failures. No variants were excluded from the pre-registered analyses.

## H1 — Primary superiority hypothesis (paired DeLong, Stratum 2)

On the binary splice-affecting outcome, both AlphaGenome scoring rules achieved high discrimination but neither exceeded SpliceAI at nominal α (Figure 1; Table 1). Under Rule A (composite), AlphaGenome achieved a receiver-operating-characteristic area under the curve (ROC-AUC) of 1.0000 against SpliceAI's 0.9819. The paired DeLong difference in AUC was +0.0181 in favour of AlphaGenome (SE 0.0151; 95% CI −0.0115 to +0.0477; z = 1.20; two-sided p = 0.231). Under Rule B (per-junction, lung/bronchial context), AlphaGenome achieved AUC 0.9922 against SpliceAI's 0.9819 (ΔAUC +0.0103; SE 0.0142; 95% CI −0.0176 to +0.0382; z = 0.726; two-sided p = 0.468). The pre-specified superiority test was not met under either rule.

## H2 — Secondary non-inferiority hypothesis (paired DeLong, Stratum 2, δ = 0.05)

AlphaGenome was non-inferior to both comparators at the pre-registered non-inferiority margin under both scoring rules (Table 1). Against SpliceAI, the lower bound of the 95% CI on the paired ΔAUC was −0.0115 under Rule A and −0.0176 under Rule B — both well inside the pre-specified δ = −0.05. Against Pangolin, AlphaGenome AUC was 1.0000 (Rule A) or 0.9922 (Rule B) versus Pangolin's 0.9845; the ΔAUC lower CI bound was −0.0113 under Rule A and comparable under Rule B, both clearing the non-inferiority margin. All four comparisons therefore support the pre-registered non-inferiority claim, and the finding is robust across scoring rules.

## H3 — Descriptive analysis of quantitative splicing (Spearman rank correlation, Stratum 1)

On the eight quantitative variants, the two scoring rules gave point estimates of the same sign but materially different magnitudes (Figure 2; Table 1). Rule A produced a Spearman rank correlation of ρ = −0.659 with the fraction of normally spliced transcript (10,000-iteration bootstrap 95% CI [−0.98, +0.10]). Rule B produced ρ = −0.287 (95% CI [−0.85, +0.61]). Both estimates matched the pre-specified negative direction. For reference, SpliceAI and Pangolin produced ρ = −0.32 (95% CI [−0.92, +0.57]) and ρ = −0.29 (95% CI [−0.98, +0.67]) respectively (Figure 2, centre and right panels). No formal significance claim is made at this sample size (n = 8); H3 is descriptive per the frozen protocol. Under Rule A, AlphaGenome's monotone tracking of residual splicing is approximately twice the magnitude of the comparators; under Rule B, it is comparable to the comparators. The bootstrap confidence intervals cross zero under both rules, and the point estimates are not statistically distinguishable from one another at this sample size.

## Table 1 — Parallel-rule summary

| Hypothesis | Metric | Rule A (composite) | Rule B (per-junction, cell-context) | Comparator |
|---|---|---|---|---|
| H1 vs SpliceAI | AUC (95% CI) | 1.0000 | 0.9922 | SpliceAI 0.9819 |
| H1 vs SpliceAI | ΔAUC (95% CI); p | +0.0181 (−0.0115 to +0.0477); 0.231 | +0.0103 (−0.0176 to +0.0382); 0.468 | — |
| H1 decision | Superiority | not met | not met | — |
| H2 vs SpliceAI | Non-inferiority (δ = 0.05) | met | met | — |
| H2 vs Pangolin | AUC | 1.0000 | 0.9922 | Pangolin 0.9845 |
| H2 vs Pangolin | ΔAUC; p | +0.0155; ≈ 0.24 | +0.0078; 0.328 | — |
| H2 decision | Non-inferiority (δ = 0.05) | met | met | — |
| H3 Stratum 1 (n = 8) | Spearman ρ (bootstrap 95% CI) | −0.659 (−0.98 to +0.10) | −0.287 (−0.85 to +0.61) | SpliceAI −0.32; Pangolin −0.29 |
| H3 direction | Negative point estimate | met | met | met |

## Divergence between scoring rules on Stratum 2 and Stratum 1

Rule A achieves perfect rank separation on Stratum 2; Rule B does not. The AUC loss under Rule B (0.9922 versus 1.0000) is arithmetically explained by a single splice-affecting variant (Joynt_E02, Rule B score 0.221) ranking below three no-effect variants (Joynt_E03/E04/E05, Rule B scores 0.346/0.427/0.502) — three wrong pairs out of 43 × 9 = 387, giving an AUC loss of 3/387 = 0.0078 and a resulting AUC of 0.9922. The three no-effect variants share a distinguishing feature under Rule B: the closest lung/bronchial-ontology canonical junction lies 883, 884, and 2 892 base pairs from the variant respectively. The pre-registered per-junction rule scores those variants against a distant junction whose per-junction absolute score is not near zero, producing a score higher than the score assigned to Joynt_E02 (whose in-context junction is 215 bp away). Under Rule A the composite maximum across all tracks and genes correctly assigns Joynt_E02 a higher score than any of the three no-effect variants.

The same mechanism accounts for the Stratum 1 divergence. Leave-one-out analysis on the eight Stratum 1 splice-affecting variants (plus the wild-type control) shows that two variants drive the difference in ρ between rules. Dropping S1_01_Masvidal (40% normally spliced, Rule B score 6.28, closest in-context junction 164 bp away) raises Rule B's ρ magnitude from −0.29 to −0.61, matching Rule A. Dropping the wild-type control S1_06_Joynt_WT (100% normally spliced, Rule B score 0.502, closest in-context junction 2 892 bp away) reduces Rule B's ρ magnitude from −0.29 to −0.04. Under Rule A the wild-type control receives a score of 0.19 (appropriately near zero); under Rule B it receives a score of the same order of magnitude as several moderately splice-affecting variants.

The unified pattern is that when the pre-declared cell-context ontology set contains no canonical junction within a few hundred base pairs of the variant, the closest-junction fallback selects a distant junction whose per-junction absolute score is not zero and is not a reliable measure of the variant's true splicing impact. Rule A's aggregation across all tracks and genes launders this out because the composite is a maximum over the full track/gene set; a spurious high in one distant junction does not swamp the true signal in the near-variant region. A per-variant table listing chosen-junction distances, per-rule scores, and outcome labels is provided as Supplementary Table S1 (`results/ruleB/DIVERGENCE_DIAGNOSIS.md` in the repository).

Because the two rules were pre-declared as parallel primary analyses (protocol Amendment #2), no unqualified claim of perfect discrimination is made: perfect rank separation is Rule-A-specific. The claim that AlphaGenome is non-inferior to SpliceAI and Pangolin on this benchmark is robust across scoring rules; the claim that AlphaGenome discriminates perfectly on Stratum 2, and the claim that its monotone tracking of quantitative residual splicing is roughly twice the magnitude of the comparators, are Rule A findings and are not reproduced under Rule B.

## Pre-registered floor and null

The pre-registered pipeline controls, established before any predictor scoring, remained the reference lower bounds: the distance-to-canonical-splice-site control reached ROC-AUC 0.683 on Stratum 2, and a fixed-seed 1 000-permutation null centred at 0.500. All three predictors exceeded the distance-only floor by a comfortable margin on Stratum 2 under either AlphaGenome scoring rule.

## Summary

Among three sequence-based splice predictors scored against a frozen public CFTR benchmark, AlphaGenome was non-inferior to SpliceAI and Pangolin at n = 52 under both pre-specified scoring rules; the non-inferiority claim held robustly across scoring rules. Perfect rank separation was achieved only under the DeepMind-recommended composite score (Rule A) and not under the pre-registered per-junction score in a lung/bronchial cell context (Rule B), where a single splice-affecting variant misranked below three no-effect variants whose closest in-context junctions were more than 880 bp from the variant. On the small quantitative stratum (n = 8), Rule A's monotone tracking of residual splicing was approximately twice the magnitude of the comparators; Rule B's was comparable, and the two rule-specific point estimates were not statistically distinguishable at this sample size. The systematic divergence between the two rules on this dataset identifies a substantive methodological property of per-junction scoring in a restricted cell-context ontology set: the closest-junction fallback becomes unreliable when there is no annotated in-context junction near the variant. These results are consistent with AlphaGenome being at least as accurate as established predictors on canonical CFTR splice variation and identify a scoring-rule choice that materially affects reported discrimination on this benchmark, warranting confirmation in larger and more junction-dense variant sets.

---

## Figure captions

**Figure 1.** Receiver-operating-characteristic (ROC) curves on Stratum 2 (n = 52; 43 splice-affecting, 9 no-effect) for AlphaGenome (Rule A composite, solid; Rule B per-junction, dashed), SpliceAI 1.3.1, and Pangolin. Shaded bands are 10 000-iteration stratified-bootstrap 95% confidence intervals on the true-positive rate at each false-positive rate. AUCs and their bootstrap 95% CIs are annotated in the legend. Dashed diagonal: chance performance.

**Figure 2.** Predictor score versus fraction of normally spliced transcript on Stratum 1 (n = 8) for AlphaGenome under Rule A (top-left, blue) and Rule B (top-right, orange), SpliceAI (bottom-left, red), and Pangolin (bottom-right, green). Solid lines show a monotone rank-based fit for visual reference; reported Spearman ρ and 10 000-iteration bootstrap 95% confidence intervals are annotated in each panel title. The pre-specified expected direction (negative ρ: higher predictor score tracking lower normal-transcript fraction) is met by both AlphaGenome scoring rules and by both comparators.
