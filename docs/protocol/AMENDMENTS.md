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

*None yet — protocol is at v1.0 as registered on 2026-09-13.*

## Predictor scoring log (informational, not amendments)

### 2026-09-20 — SpliceAI 1.3.1 and Pangolin (tkzeng, main branch, cloned 2026-09-20) both scored

Both comparator predictors were run against `data/variants_for_scoring.vcf` on the frozen chr7 GRCh38 FASTA (Ensembl release 110). SpliceAI was run via the packaged CLI with a NumPy-2 compatibility shim (`src/predict/spliceai_shims.py`); Pangolin was run via a pyvcf3 compatibility shim (`src/predict/pangolin_shims.py`) against a GENCODE v44 chr7 annotation database built with `Ensembl_canonical` filter.

Outputs are recorded in `results/predictions/RESULTS_MANIFEST.tsv` with their SHA-256 hashes:

- `results/predictions/spliceai.tsv` — 60/60 variants scored, no unscorable variants.
- `results/predictions/pangolin.tsv` — 60/60 variants scored, no unscorable variants.

These are results, not inputs, so they are recorded in a separate results manifest rather than `data/HASH_MANIFEST.tsv`. The primary analysis has not yet been run against the AlphaGenome scores — those remain gated on the AlphaGenome licence decision (cron `ad481e61`, fires 2026-09-27).

A smoke run using the SpliceAI and Pangolin scores as if they were the two comparators (mechanics test only, not a pre-registered analysis) confirmed the pipeline runs end-to-end and produces sensible AUCs. Those smoke outputs live under `results/smoke/` (gitignored).
