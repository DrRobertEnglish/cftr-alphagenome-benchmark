#!/usr/bin/env python3
"""Build an updated Colab notebook that scores the frozen variants using BOTH
the DeepMind-recommended composite (as before) AND the OSF-preregistered rule
(new cells from prereg_scoring_cells.py).

Reads `scripts/colab_alphagenome_scoring.ipynb`, inserts three new cells
(8b: SPLICE_JUNCTIONS scoring, 8c: pre-reg rule computation, 11b: download the
pre-reg TSV), and writes `scripts/colab_alphagenome_scoring_v2.ipynb`.

Run: python scripts/build_updated_colab_notebook.py
"""

from __future__ import annotations
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
IN_NB = REPO / 'scripts' / 'colab_alphagenome_scoring.ipynb'
OUT_NB = REPO / 'scripts' / 'colab_alphagenome_scoring_v2.ipynb'
CELLS_PY = REPO / 'scripts' / 'prereg_scoring_cells.py'


def load_cell_sources_from_module() -> tuple[str, str, str]:
    """Import cell code strings from prereg_scoring_cells.py without executing the module."""
    src = CELLS_PY.read_text()
    # Parse out the three assigned string literals. Use ast to be robust to future edits.
    import ast
    tree = ast.parse(src)
    cells: dict[str, str] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
            if name in ('CELL_8B_CODE', 'CELL_8C_CODE', 'CELL_11B_CODE') and isinstance(node.value, ast.Constant):
                cells[name] = node.value.value
    if set(cells) != {'CELL_8B_CODE', 'CELL_8C_CODE', 'CELL_11B_CODE'}:
        raise RuntimeError(f'Failed to load all three cells: got {sorted(cells)}')
    return cells['CELL_8B_CODE'], cells['CELL_8C_CODE'], cells['CELL_11B_CODE']


def make_code_cell(source: str) -> dict:
    """Colab-compatible code cell."""
    return {
        'cell_type': 'code',
        'metadata': {},
        'execution_count': None,
        'outputs': [],
        'source': source.lstrip('\n').splitlines(keepends=True),
    }


def make_markdown_cell(source: str) -> dict:
    return {
        'cell_type': 'markdown',
        'metadata': {},
        'source': source.lstrip('\n').splitlines(keepends=True),
    }


def main() -> None:
    nb = json.loads(IN_NB.read_text())
    cell_8b, cell_8c, cell_11b = load_cell_sources_from_module()

    # Find the last SPLICE-scoring code cell (cell 8 in v1) so we can insert after it.
    insert_after_idx = None
    for i, c in enumerate(nb['cells']):
        src = ''.join(c.get('source', []))
        if c['cell_type'] == 'code' and 'model.score_variant' in src and 'splicing_scorers' in src:
            insert_after_idx = i
    if insert_after_idx is None:
        raise RuntimeError('Could not locate the main scoring cell (looking for splicing_scorers + score_variant).')

    # Build the insertion block: markdown header + code cell + markdown header + code cell.
    new_cells_after_scoring = [
        make_markdown_cell(
            '## 4b. Pre-registration-faithful scoring (per-junction, cell-context filtered)\n\n'
            'This block implements the scoring rule that was pre-registered on the OSF: '
            '`|usage(REF) - usage(ALT)|` at the canonical junction closest to the variant, '
            'in a default lung / bronchial epithelial cell context. It runs alongside the '
            'DeepMind-recommended composite already scored above; the paper reports both as '
            'parallel primary analyses (see `docs/protocol/AMENDMENTS.md` amendment #2, '
            '2026-09-20).'
        ),
        make_code_cell(cell_8b),
        make_markdown_cell(
            '## 4c. Compute pre-registered score per variant and write TSV\n\n'
            'Selects, for each variant, the junction closest to the variant position within '
            'the pre-declared cell-context ontology set, and takes `max(|raw_score|)` across '
            'matching tracks at that junction. Writes `alphagenome_prereg.tsv` and reports '
            'the SHA-256.'
        ),
        make_code_cell(cell_8c),
    ]

    # Find the download cell (last code cell in v1) so we add the pre-reg download after it.
    download_idx = None
    for i, c in enumerate(nb['cells']):
        src = ''.join(c.get('source', []))
        if c['cell_type'] == 'code' and 'files.download' in src:
            download_idx = i
    if download_idx is None:
        raise RuntimeError('Could not locate the files.download cell in v1 notebook.')

    # Splice. Insert 4b/4c right after the scoring cell.
    new_cells = (
        nb['cells'][: insert_after_idx + 1]
        + new_cells_after_scoring
        + nb['cells'][insert_after_idx + 1 :]
    )
    # Find the (new) index of the download cell after insertion.
    for i, c in enumerate(new_cells):
        src = ''.join(c.get('source', []))
        if c['cell_type'] == 'code' and 'files.download' in src:
            new_download_idx = i
            break
    # Append pre-reg download markdown + code cells right after the composite download.
    new_cells = (
        new_cells[: new_download_idx + 1]
        + [
            make_markdown_cell(
                '## 6b. Download `alphagenome_prereg.tsv`\n\n'
                'This downloads the pre-registered-rule TSV to your local machine. Send it '
                'back to the chat session alongside `alphagenome.tsv`; the two files together '
                'let the paper report both scoring rules.'
            ),
            make_code_cell(cell_11b),
        ]
        + new_cells[new_download_idx + 1 :]
    )
    nb['cells'] = new_cells

    # Update the intro markdown to mention the second scoring rule.
    for i, c in enumerate(nb['cells']):
        src = ''.join(c.get('source', []))
        if c['cell_type'] == 'markdown' and src.startswith('# AlphaGenome scoring'):
            new_intro = (
                '# AlphaGenome scoring - 60 frozen CFTR variants\n\n'
                'This notebook scores the pre-registered variant list from the '
                '[cftr-alphagenome-benchmark](https://github.com/DrRobertEnglish/cftr-alphagenome-benchmark) '
                'repository using the official AlphaGenome SDK, computes **two** scoring rules:\n\n'
                '1. The DeepMind-recommended **merged splicing score** '
                '(`max(splice_sites) + max(splice_site_usage) + max(splice_junctions)/5`), '
                'aggregated across all tracks and genes.\n'
                '2. The **OSF-preregistered rule**: `|usage(REF) - usage(ALT)|` at the '
                'canonical junction closest to the variant, in a default lung / bronchial '
                'epithelial cell context.\n\n'
                'It writes two hash-verified TSVs (`alphagenome.tsv` and '
                '`alphagenome_prereg.tsv`) that you drop back into the local repo for the '
                'H1/H2/H3 analyses. Both rules are reported as parallel primary analyses in '
                'the paper (see `docs/protocol/AMENDMENTS.md` amendment #2).\n\n'
                'Run each cell in order. Total runtime is roughly 10-25 minutes depending on '
                'Colab load. The API key never leaves your Colab runtime.\n\n'
                '**Pre-registration:** [osf.io/5hvgf](https://osf.io/5hvgf/) - DOI '
                '[10.17605/OSF.IO/6PGX8](https://doi.org/10.17605/OSF.IO/6PGX8)\n\n'
                '**Licence:** AlphaGenome output is subject to the [AlphaGenome Output Terms '
                'of Use](https://www.alphagenomedocs.com/licenses/). Computed statistics '
                '(ROC-AUC, DeLong CIs, Spearman rho) are summary statistics, not derivative '
                'models. No machine learning model is trained on AlphaGenome output.\n'
            )
            nb['cells'][i] = make_markdown_cell(new_intro)
            break

    OUT_NB.parent.mkdir(parents=True, exist_ok=True)
    OUT_NB.write_text(json.dumps(nb, indent=1) + '\n')
    print(f'Wrote {OUT_NB}')
    print(f'Cells: {len(nb["cells"])} (was {len(json.loads(IN_NB.read_text())["cells"])})')


if __name__ == '__main__':
    main()
