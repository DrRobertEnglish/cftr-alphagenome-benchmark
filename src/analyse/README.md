# Analysis pipeline

Frozen implementation of the pre-registered decision rules for H1, H2, H3.
Every runner reads the frozen variant list plus one or two predictor score
TSVs and produces derived outputs only.

## Modules

- `io.py` — shared readers (`load_variants`, `load_predictor_scores`) and
  frame builders (`build_stratum2_frame`, `build_stratum1_frame`). Enforces
  the pre-registered inclusion rules in one place.
- `paired_delong.py` — Sun & Xu (2014) fast paired-DeLong implementation, the
  Python cross-check for the R `pROC::roc.test(method="delong", paired=TRUE)`
  primary implementation. `superiority_check` applies the H1 rule;
  `non_inferiority_check` applies the H2 rule with margin δ.
- `run_primary_analysis.py` — end-to-end runner for H1 (superiority) and H2
  (non-inferiority) on Stratum 2, plus the exploratory stratified-bootstrap
  AUC CIs documented in `docs/protocol/WRITEUP_NOTES.md`.
- `run_stratum1.py` — end-to-end runner for H3 (descriptive Spearman ρ with
  10 000-iter bootstrap CI) on Stratum 1.

## Predictor score TSV schema

Both runners expect a TSV with at minimum:

| column | type | notes |
|---|---|---|
| `variant_id` | string | must match the `variant_id` column in `data/variant_list_frozen.tsv` |
| `<score column>` | float | the predictor's continuous score for the variant |

The score column name is passed with `--score-column` (default `score`).
Missing values are rejected. Rows for variants not in the frozen list are
ignored. Rows for variants that are in the frozen list but not scored by a
predictor are reported in the coverage block and dropped from the analysis.

## Running the primary analysis (H1 + H2)

    python -m src.analyse.run_primary_analysis \
        --variants data/variant_list_frozen.tsv \
        --scores-a results/predictions/alphagenome.tsv --label-a AlphaGenome \
        --scores-b results/predictions/spliceai.tsv    --label-b SpliceAI \
        --score-column score \
        --output results/primary_alphagenome_vs_spliceai.tsv

Emits three files at the output path:
- `<name>.tsv` — flat machine-readable metrics.
- `<name>.json` — nested JSON with full coverage, DeLong result, H1 and H2
  decisions, bootstrap CIs, and provenance.
- `<name>.md` — plain-language decision block for pasting into the paper.

## Running the H3 analysis

    python -m src.analyse.run_stratum1 \
        --variants data/variant_list_frozen.tsv \
        --scores results/predictions/alphagenome.tsv \
        --label AlphaGenome --score-column score \
        --output results/stratum1_alphagenome.tsv

Same output pattern: `.tsv`, `.json`, `.md`.

## Smoke test

Run `bash scripts/smoke_analysis_pipeline.sh` to exercise the whole pipeline
end-to-end against the splice-site-distance floor benchmark (mechanics test
only, not a scientific comparison). Outputs go to `results/smoke/` which is
gitignored.

## What is NOT in these runners by design

- No API calls to AlphaGenome, SpliceAI, or Pangolin. Predictor scores are
  read from static TSV files. This keeps the analysis reproducible offline
  and independent of any external service's availability.
- No implicit inclusion of variants outside the frozen list. Add-variant
  amendments would require a version bump of the OSF pre-registration.
- No re-computation of the pre-registered decision rules from the data. H1
  is superiority with p < 0.05; H2 is non-inferiority with δ = 0.05.
