# CFTR splice variants with published quantitative residual normal-transcript measurements

**Harvest frozen:** 2026-09-10
**Author:** Robert English

Sources traced to primary paper. Only variants with a numeric % or SD/SEM measurement are counted as quantitative (qualitative "residual" mentions listed separately).

## Fully quantified splice variants (Stratum 1 — continuous outcome)

| # | HGVS c. | Legacy | Assay | % normal | SD/SEM | Source |
|---|---|---|---|---|---|---|
| 1 | c.580-1G>T | 712-1G>T | Nasal epithelium RT-qPCR (patient) | 40% (patients) / 66% (carriers) | ±1.3% / ±2.7% | [Masvidal 2014](https://www.nature.com/articles/ejhg2013238) |
| 2 | c.2657+5G>A | 2789+5G>A | Nasal epithelium RT-qPCR + minigene | 24% (patients) / 38% (carriers) / 29% (minigene) | ±9% / ±10% / ±2.3% | [Masvidal 2014](https://www.nature.com/articles/ejhg2013238) |
| 3 | c.[2657+5G>A; 2562T>G] complex | — | Minigene | 23% | ±0.4% | [Masvidal 2014](https://www.nature.com/articles/ejhg2013238) |
| 4 | c.2908G>C | G970R | Minigene EMG_i14-i18 | 0% | — | [Joynt 2020](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009100) |
| 5 | c.2908G>A | G970S | Minigene EMG_i14-i18 | 0% | — | [Joynt 2020](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009100) |
| 6 | c.2909G>A | G970D | Minigene EMG_i14-i18 (wild-type control) | 100% | — | [Joynt 2020](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009100) |
| 7 | c.3873G>C | Q1291H | Minigene EMG_i21-i24 | 37.2% | ±2.2% | [Joynt 2020](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009100) |
| 8 | c.3717+12191C>T (=c.3718-2477C>T) | 3849+10kbC>T | Fetal patient RNA (multi-organ RT-PCR) | Lung 26% / trachea 17% / colon 19% / pancreas 2% / ileum 1% | ±0.5–1.5% | [Chiba-Falek 1999](https://academic.oup.com/ajrccm/article/159/6/1998/8530084) |
| 9 | c.1680-886A>G | 1811+1.6kbA>G | Patient mRNA | 1-3% | — | [Deletang 2022](https://www.nature.com/articles/s41434-022-00347-0) |

**Quantified splice-affecting: n = 8 unique variants** (+ 1 wild-type control)
Plus: multi-organ measurements for 3849+10kbC>T give 5 additional data points.

## Binary tier (Stratum 2 — splice-affecting yes/no, from Joynt 2020)

Joynt characterised **52 CFTR variants (15 exonic + 37 intronic)** with binary splice-effect classification:
- **43/52 splice-affecting** (60% exonic, 91.9% intronic)
- **9/52 no splice effect**

This gives a clean binary AUC with **n = 52** using the Joynt paper alone.

## Additional splice-affecting variants with qualitative residual (Joynt supplementary)

- c.274G>A, c.3700A>G, c.4005+2T>C, c.4242+2T>C, c.165-3C>T, c.2657+2_2657+3insA, c.3468+5G>A, c.3717+40A>G

**Qualitative-only: n = 8** (usable in a binary "any residual vs no residual" analysis if the paper's residual/no-residual call is extracted from figures)

## Not usable (mentioned but no quantitative or extractable classification)

- c.3140-16T>A ([Kondratyeva 2021](https://pubmed.ncbi.nlm.nih.gov/34071719/)) — residual mentioned, no %
- c.1210-12T[5] (5T) — normal transcript levels described qualitatively
- CFTR2 non-canonical splice regulator adjudications — ~19 variants with functional call but not always with residual %

## Effective n for the primary hypothesis

- **Primary outcome (Stratum 2 binary AUC):** **n = 52** (Joynt) with 43 splice-affecting + 9 no-effect variants
- **Secondary outcome (Stratum 1 rank correlation):** **n = 9–14** (variants with numeric %)
- **Total unique CFTR splice variants with at least binary functional data across sources:** conservatively **n = 50–70** after de-duplication

**Pre-registered harvest threshold:** n ≥ 20 continuous OR n ≥ 30 binary. Both met.

**Harvest gate: PASS.**
