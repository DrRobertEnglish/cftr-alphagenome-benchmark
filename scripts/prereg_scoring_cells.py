"""Pre-registration-faithful AlphaGenome scoring — additional cells for the Colab notebook.

This module contains the code for two additional Colab cells that operationalise the
OSF-preregistered scoring rule ("per-junction usage delta at the canonical junction
closest to the variant, default lung/bronchial epithelial cell context"). The cells
run alongside the DeepMind-recommended composite scoring already in the notebook,
producing a second TSV (`alphagenome_prereg.tsv`) with per-variant pre-registered
scores. This lets us report both rules as parallel primary analyses per the
2026-09-20 amendment (see docs/protocol/AMENDMENTS.md).

Design notes (derived from the AlphaGenome SDK documentation as of 2026-09-20):

1. The SDK's SPLICE_SITE_USAGE scorer returns per-track magnitudes aggregated over
   sites within a gene mask. It does not expose per-junction, per-position outputs
   that would let us implement "|usage(REF) - usage(ALT)| at the canonical junction
   closest to variant" using SPLICE_SITE_USAGE alone.

2. The SpliceJunctionScorer (SPLICE_JUNCTIONS output type) does expose per-junction,
   per-track data via the `junction_Start` and `junction_End` columns in the tidy
   DataFrame. This is the AlphaGenome SDK object that most closely maps to the
   pre-registered rule's "at the canonical junction closest to variant" clause.

3. "Default lung / bronchial epithelial cell context" is operationalised as the
   union of tracks whose ontology CURIE is in the pre-declared set:
     - UBERON:0002048 (lung)
     - UBERON:0002185 (bronchus)
     - CL:0002145 (ciliated columnar cell of tracheobronchial tree)
     - CL:1000271 (lung ciliated cell)
     - CL:0002632 (epithelial cell of lower respiratory tract)
   If no tracks match this set for a variant, the pre-reg-faithful score is recorded
   as NA and the variant is excluded from pre-registered-rule analyses. This
   fallback behaviour is documented in the amendment.

4. The per-variant score is: for the junction closest to the variant position
   (measured as min(|variant_pos - junction_Start|, |variant_pos - junction_End|),
   summed over the two junction endpoints — i.e. distance from variant to the
   nearest endpoint of the junction), take max(|raw_score|) across the pre-declared
   ontology-filtered tracks. This is a magnitude (non-signed), matching the
   pre-registered "absolute change" wording.

5. Deliberate limitations:
   - "Canonical" is operationalised as "closest to the variant" because AlphaGenome
     does not carry a canonical-junction flag; this is the same approximation used
     by SpliceAI's Δ-max operationalisation of "closest splice site".
   - The junction-distance tie-break rule is: if two junctions are equidistant,
     use the lower junction_Start.
   - Cell-context ontology set is fixed by this module and pre-declared before
     any scoring; it is not adjusted post-hoc.

Cell 8b (code): score every variant with SpliceJunctionScorer, keeping the tidy
DataFrame per variant so we can select the junction after the fact.
Cell 8c (code): compute the pre-registered rule score per variant and write TSV.

These cells DO NOT require re-scoring with SPLICE_SITES or SPLICE_SITE_USAGE — the
main scoring cell (cell 8) already covers those. They only need one additional
API call per variant, using the SpliceJunctionScorer specifically.
"""

# ==============================================================================
# CELL 8b — score every variant with the per-junction SPLICE_JUNCTIONS scorer
# ==============================================================================

CELL_8B_CODE = '''
# 4b. Pre-registration-faithful scoring — per-junction SPLICE_JUNCTIONS output
#
# The OSF-preregistered rule scores each variant as |usage(REF) - usage(ALT)|
# at the canonical junction closest to the variant, in a default lung/bronchial
# epithelial cell context. The AlphaGenome SDK's SplicerJunctionScorer exposes
# per-junction, per-track scores that let us implement this. We score every
# variant with only the SpliceJunctionScorer here (the other two scorers are
# already covered by cell 8) and keep the tidy DataFrame per variant.
#
# NOTE: This cell adds ~1 API call per variant and takes ~2-5 minutes total.

import time, pandas as pd

# Pre-declared cell-context ontology set (see prereg_scoring_cells.py docstring
# for provenance). This must be frozen before any scoring; do not change it
# post-hoc.
PREREG_CELL_CONTEXT_ONTOLOGIES = frozenset({
    'UBERON:0002048',  # lung
    'UBERON:0002185',  # bronchus
    'CL:0002145',      # ciliated columnar cell of tracheobronchial tree
    'CL:1000271',      # lung ciliated cell
    'CL:0002632',      # epithelial cell of lower respiratory tract
})

junction_scorer = vs.RECOMMENDED_VARIANT_SCORERS['SPLICE_JUNCTIONS']

per_variant_junction_dfs = {}   # variant_id -> tidy DataFrame or None
junction_scoring_errors = []

for i, row in enumerate(variants, 1):
    t0 = time.time()
    vid = row['variant_id']
    try:
        variant = genome.Variant(
            chromosome=f"chr{row['chrom']}",
            position=int(row['grch38_pos']),
            reference_bases=row['ref'],
            alternate_bases=row['alt'],
        )
        interval = variant.reference_interval.resize(SEQ_LEN)
        scores = model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=[junction_scorer],
        )
        df = vs.tidy_scores([scores])
        per_variant_junction_dfs[vid] = df
        dt = time.time() - t0
        n_rows = len(df) if df is not None else 0
        print(f'[{i:>2}/60] {vid}: OK ({n_rows} junction-track rows) in {dt:.1f}s')
    except Exception as e:
        per_variant_junction_dfs[vid] = None
        junction_scoring_errors.append((vid, f'{type(e).__name__}: {str(e)[:200]}'))
        print(f'[{i:>2}/60] {vid}: FAILED ({type(e).__name__})')

print()
print(f'Scored: {sum(1 for v in per_variant_junction_dfs.values() if v is not None)}/60')
print(f'Failed: {len(junction_scoring_errors)}/60')
if junction_scoring_errors:
    print('Failures:')
    for vid, err in junction_scoring_errors:
        print(f'  {vid}: {err}')
'''

# ==============================================================================
# CELL 8c — compute the pre-registered rule score per variant, write TSV
# ==============================================================================

CELL_8C_CODE = '''
# 4c. Compute pre-registration-faithful score per variant
#
# For each variant with a successful SpliceJunctionScorer output:
#   1. Filter tidy DataFrame to rows whose ontology_curie is in
#      PREREG_CELL_CONTEXT_ONTOLOGIES.
#   2. If no rows match, mark this variant as NA_no_cell_context.
#   3. Otherwise, for each remaining row compute
#         d = min(|variant_pos - junction_Start|, |variant_pos - junction_End|)
#      and pick the row(s) with the minimum d (canonical-junction proxy).
#   4. Among those rows, take max(|raw_score|) across tracks. This is the
#      pre-reg-faithful score.
#
# The output TSV alphagenome_prereg.tsv has one row per variant with columns
# for the score, the chosen junction coordinates, the number of matching
# in-context tracks, and the closest-junction distance in bp.

import csv, hashlib, math

def prereg_score_for_variant(vid, variant_pos, df):
    """Return dict of pre-reg-faithful score fields, or NA fields if not scorable."""
    out = {
        'variant_id': vid,
        'prereg_score': '',
        'chosen_junction_start': '',
        'chosen_junction_end': '',
        'chosen_junction_distance_bp': '',
        'n_in_context_tracks': '',
        'note': '',
    }
    if df is None or len(df) == 0:
        out['note'] = 'NA_no_junction_scoring'
        return out
    if 'ontology_curie' not in df.columns:
        out['note'] = 'NA_missing_ontology_column'
        return out
    if 'junction_Start' not in df.columns or 'junction_End' not in df.columns:
        out['note'] = 'NA_missing_junction_columns'
        return out
    in_context = df[df['ontology_curie'].isin(PREREG_CELL_CONTEXT_ONTOLOGIES)].copy()
    n_in_context_tracks = int(in_context['track_name'].nunique()) if 'track_name' in in_context.columns else len(in_context)
    if len(in_context) == 0:
        out['note'] = 'NA_no_cell_context'
        out['n_in_context_tracks'] = 0
        return out
    # Distance from variant to nearest endpoint of each junction.
    in_context['distance'] = in_context.apply(
        lambda r: min(abs(variant_pos - int(r['junction_Start'])),
                      abs(variant_pos - int(r['junction_End']))),
        axis=1,
    )
    dmin = in_context['distance'].min()
    closest = in_context[in_context['distance'] == dmin]
    # Tie-break: lowest junction_Start
    closest = closest.sort_values(['junction_Start', 'junction_End'])
    j_start = int(closest.iloc[0]['junction_Start'])
    j_end = int(closest.iloc[0]['junction_End'])
    # Score: max(|raw_score|) across tracks matching that junction, in-context.
    matched = in_context[
        (in_context['junction_Start'] == j_start) &
        (in_context['junction_End'] == j_end)
    ]
    score = float(matched['raw_score'].abs().max())
    out['prereg_score'] = f'{score:.6f}'
    out['chosen_junction_start'] = str(j_start)
    out['chosen_junction_end'] = str(j_end)
    out['chosen_junction_distance_bp'] = str(int(dmin))
    out['n_in_context_tracks'] = str(n_in_context_tracks)
    out['note'] = ''
    return out


prereg_results = []
for row in variants:
    vid = row['variant_id']
    variant_pos = int(row['grch38_pos'])
    df = per_variant_junction_dfs.get(vid)
    r = prereg_score_for_variant(vid, variant_pos, df)
    r['chrom'] = row['chrom']
    r['pos'] = row['grch38_pos']
    r['ref'] = row['ref']
    r['alt'] = row['alt']
    prereg_results.append(r)

OUT_PATH_PREREG = 'alphagenome_prereg.tsv'
fieldnames_prereg = [
    'variant_id', 'chrom', 'pos', 'ref', 'alt',
    'prereg_score',
    'chosen_junction_start', 'chosen_junction_end', 'chosen_junction_distance_bp',
    'n_in_context_tracks',
    'note',
]
with open(OUT_PATH_PREREG, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fieldnames_prereg, delimiter='\\t', lineterminator='\\n')
    w.writeheader()
    for r in prereg_results:
        w.writerow({k: r.get(k, '') for k in fieldnames_prereg})

with open(OUT_PATH_PREREG, 'rb') as f:
    sha_prereg = hashlib.sha256(f.read()).hexdigest()
print(f'Wrote {OUT_PATH_PREREG} - SHA-256 {sha_prereg}')

n_scored = sum(1 for r in prereg_results if r['prereg_score'])
n_no_context = sum(1 for r in prereg_results if r['note'] == 'NA_no_cell_context')
n_other_na = sum(1 for r in prereg_results if r['note'] and r['note'] != 'NA_no_cell_context')
print(f'Pre-reg scored: {n_scored}/60')
print(f'NA (no in-context tracks): {n_no_context}/60')
print(f'NA (other): {n_other_na}/60')
if n_other_na:
    for r in prereg_results:
        if r['note'] and r['note'] != 'NA_no_cell_context':
            print(f'  {r["variant_id"]}: {r["note"]}')
'''


# ==============================================================================
# CELL 11b — download the pre-reg TSV alongside the composite TSV
# ==============================================================================

CELL_11B_CODE = '''
from google.colab import files
files.download(OUT_PATH_PREREG)
'''


if __name__ == '__main__':
    # Print the cells so they can be pasted into the notebook (or fed to a builder).
    print('# ============= CELL 8b =============')
    print(CELL_8B_CODE)
    print('# ============= CELL 8c =============')
    print(CELL_8C_CODE)
    print('# ============= CELL 11b =============')
    print(CELL_11B_CODE)
