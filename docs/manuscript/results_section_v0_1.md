# Results (draft v0.1)

This section presents the pre-registered analyses defined in the OSF protocol (osf.io/5hvgf; DOI 10.17605/OSF.IO/6PGX8). All decision rules, test statistics, and confidence-interval procedures were locked before AlphaGenome scoring. Bootstrap procedures use seed 20260913 throughout. Predictor outputs, analysis code, and hash-locked inputs are archived in the accompanying repository (`DrRobertEnglish/cftr-alphagenome-benchmark`, commit `d0fb018`).

## Cohort

Sixty CFTR variants were pre-registered across two strata: Stratum 2 contained 52 variants from Joynt et al. with a binary splice-affecting label (43 splice-affecting, 9 no-effect); Stratum 1 contained eight unique variants with quantitative measurements of the fraction of normally spliced transcript. All 60 variants were scored by AlphaGenome, SpliceAI 1.3.1, and Pangolin without failures. No variants were excluded from the pre-registered analyses.

## H1 — Primary superiority hypothesis (paired DeLong, Stratum 2)

On the binary splice-affecting outcome, AlphaGenome achieved a receiver-operating-characteristic area under the curve (ROC-AUC) of 1.000 against SpliceAI's 0.982 (Figure 1). The paired DeLong difference in AUC was +0.018 in favour of AlphaGenome (SE 0.015; 95% CI −0.012 to +0.048; z = 1.20; two-sided p = 0.231). The pre-specified superiority test was therefore not met: AlphaGenome is not statistically superior to SpliceAI on this benchmark. A stratified-bootstrap 95% confidence interval on AlphaGenome's AUC covered [1.000, 1.000], reflecting the ceiling effect at this sample size.

## H2 — Secondary non-inferiority hypothesis (paired DeLong, Stratum 2, δ = 0.05)

AlphaGenome was non-inferior to both comparators at the pre-registered non-inferiority margin. Against SpliceAI the lower bound of the 95% CI on the paired ΔAUC was −0.012 (well inside δ = −0.05). Against Pangolin, AlphaGenome AUC 1.000 vs Pangolin 0.984, ΔAUC +0.016 (SE 0.014; 95% CI −0.012 to +0.043; z = 1.09; two-sided p = 0.275); the lower CI bound of −0.012 again cleared the non-inferiority margin. Both comparisons therefore support the pre-registered non-inferiority claim.

## H3 — Secondary descriptive analysis (Spearman rank correlation, Stratum 1)

On the eight quantitative variants, AlphaGenome showed a Spearman rank correlation of ρ = −0.66 with the fraction of normally spliced transcript (10,000-iteration bootstrap 95% CI [−0.98, +0.10]; Figure 2, left panel). The direction of association matched the pre-specified negative direction. For comparison, SpliceAI and Pangolin produced ρ = −0.32 (95% CI [−0.92, +0.57]) and ρ = −0.29 (95% CI [−0.98, +0.67]) respectively, both in the same pre-specified direction (Figure 2, centre and right panels). No formal significance claim is made at this sample size; the analysis is descriptive per the frozen protocol.

## Pre-registered floor and null

The pre-registered pipeline controls, established before any predictor scoring, remained the reference lower bounds: the distance-to-canonical-splice-site control reached ROC-AUC 0.683 on Stratum 2, and a fixed-seed 1,000-permutation null centred at 0.500. All three predictors exceeded the distance-only floor by a comfortable margin on Stratum 2.

## Summary

Among three sequence-based splice predictors scored against a frozen public CFTR benchmark, AlphaGenome achieved a perfect Stratum-2 discrimination point estimate but was not statistically superior to either SpliceAI or Pangolin at n = 52; the pre-registered non-inferiority claim held against both comparators. On the small quantitative stratum (n = 8), AlphaGenome's rank correlation with residual splicing was approximately twice the magnitude of the comparators in the pre-specified negative direction, with wide confidence intervals reflecting the small sample size. These results are consistent with AlphaGenome being at least as accurate as established predictors on canonical CFTR splice variation, with a potentially stronger monotone relationship to quantitative splicing effect that warrants confirmation in larger quantitative datasets.

---

## Figure captions

**Figure 1.** Receiver-operating-characteristic (ROC) curves on Stratum 2 (n = 52; 43 splice-affecting, 9 no-effect) for AlphaGenome, SpliceAI 1.3.1, and Pangolin. Solid lines show the point-estimate ROC; shaded bands are 10,000-iteration stratified-bootstrap 95% confidence intervals on the true-positive rate at each false-positive rate. AUCs and their bootstrap 95% CIs are annotated in the legend. Dashed diagonal: chance performance.

**Figure 2.** Predictor score versus fraction of normally spliced transcript on Stratum 1 (n = 8) for AlphaGenome (left, blue), SpliceAI (centre, red), and Pangolin (right, green). Solid lines show a monotone rank-based fit for visual reference; reported Spearman ρ and 10,000-iteration bootstrap 95% confidence intervals are annotated in each panel title. The pre-specified expected direction (negative ρ: higher predictor score tracking lower normal-transcript fraction) is met by all three predictors.
