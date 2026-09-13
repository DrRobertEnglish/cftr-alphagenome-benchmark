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

*None yet — protocol is at v1.0 as registered on 2026-09-13.*
