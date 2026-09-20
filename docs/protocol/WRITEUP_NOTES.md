# Write-up notes — decisions to record in the manuscript

Running list of choices made after OSF registration (2026-09-13) that are consistent with the frozen protocol but should be described plainly in the paper's Methods and (where relevant) Limitations sections. These are not protocol amendments; they are execution decisions the reader has a right to see.

Add new items to the bottom with the date; do not rewrite history.

---

## 2026-09-13 — Masvidal complex allele c.[2657+5G>A;2562T>G] excluded from the machine-readable variant list

**Where to describe:** Methods → Variant harvest and Methods → Predictor scoring.

**What happened.** The Masvidal 2014 harvest contained one complex allele — two SNVs in cis — with a minigene measurement of 23% ± 0.4% normally spliced transcript. When we materialised the OSF-registered narrative harvest into `data/variant_list_frozen.tsv` we excluded this row from the analysable set.

**Why.** All three sequence predictors compared in the primary analysis (SpliceAI, Pangolin, AlphaGenome) score a single-variant VCF-style input (one chromosome, one position, one reference allele, one alternate allele). None accepts a two-variant cis-linked complex allele as one input, and there is no principled way to combine two separate predictor scores into one composite score without changing the analysis. Scoring the two constituent SNVs independently would double-count them in Stratum 2 (c.2657+5G>A is already present as a Joynt Stratum 2 row) and would not correspond to any measurement in the ground truth.

**Effect on the pre-registered hypotheses.** None. The Stratum 1 rank-correlation analysis (H3) is the only place this variant could have entered; excluding it changes the Stratum 1 effective sample from n=9 to n=8. The pre-registered H3 threshold (n≥9 variants with numeric %) is not met with n=8 alone, but H3 is descriptive and will be reported as such with a bootstrap CI, not as a hypothesis test. The primary H1 test is unaffected because Stratum 2 (n=52 Joynt) does not include this variant.

**How to phrase it.** "One complex allele from the Masvidal 2014 harvest, comprising two cis-linked SNVs, was excluded from the analysable variant list because the single-variant predictors compared in this study cannot score two linked substitutions as one input. This exclusion was documented in the code and results are unchanged if the constituent SNVs are analysed separately as they are already present in Stratum 2."

---

## 2026-09-13 — Stratum 2 class balance is 43:9 splice-affecting to no-effect

**Where to describe:** Methods → Statistical analysis, and Limitations.

**What happened.** The 52 variants in Joynt Table 1 and Table 2 partition into 43 splice-affecting and 9 no-splice-effect on the RNA-effect column. This is inherited from the source paper's variant selection — Joynt's design deliberately over-sampled canonical splice sites where missplicing is essentially guaranteed.

**Why it matters for interpretation.** The 43:9 skew widens the null distribution for ROC-AUC on this dataset. In the frozen negative-control run (1000 label permutations, seed 20260913), the null 95% CI for AUC is [0.296, 0.698] with SD 0.106. A predictor scoring 0.68 on Stratum 2 is not clearly distinguishable from chance on this dataset.

**Effect on the pre-registered hypotheses.** None on the H1 decision rule. The pre-registered primary test is paired DeLong of AlphaGenome vs SpliceAI on the identical 52-variant set, which conditions on the shared labels and so does not depend on class-balance-driven null width. But the interpretive language must be careful: an AlphaGenome AUC that is high in absolute terms is only clinically meaningful if it beats SpliceAI by a margin the paired test can detect at α=0.05.

**Planned robustness add-on (exploratory, in Other Planned Analyses).** Report each predictor's Stratum 2 AUC alongside a stratified-bootstrap 95% CI in addition to the pre-registered paired-DeLong CI. This is exploratory and does not change the frozen H1 decision rule. It is a robustness display, not a test.

**How to phrase it.** "The Joynt 2020 variant set was enriched for canonical splice sites and therefore carries a splice-affecting to no-effect ratio of 43:9. This skew widens the null distribution for any single-predictor AUC on this dataset (empirical 95% null CI [0.296, 0.698] from 1000 label permutations with a frozen seed), so single-predictor AUCs are reported with stratified-bootstrap CIs alongside the pre-registered paired-DeLong test that compares predictors on the same 52 variants."
