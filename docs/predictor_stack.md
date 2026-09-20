# Predictor stack

The primary hypothesis (H1) compares AlphaGenome against SpliceAI. The secondary non-inferiority hypothesis (H2) compares AlphaGenome against Pangolin. Both comparator predictors are run locally so that the same VCF, FASTA, and hash-manifested inputs are seen by every model, and so that the run can be reproduced from a clean checkout.

This document records the exact environment those comparators run in, and the compatibility work needed to keep 2020-era released predictor codebases compiling on 2026 dependency stacks.

## Isolated Python environment

Both predictors need pins that conflict with the sandbox default. They live in a dedicated `.venv-predictors/` virtualenv at repo root (gitignored).

- Python 3.12 (SpliceAI's TensorFlow wheels don't cover 3.14 yet).
- `setuptools<81` (SpliceAI 1.3.1 still imports the deprecated `pkg_resources`; setuptools 84 no longer ships it).
- `tensorflow==2.21.0`, `numpy>=2`, `scipy>=1.13`, `keras==3.15.1` for SpliceAI.
- `torch` (CPU wheels from the PyTorch index) for Pangolin.
- `pyvcf3` (maintained fork of the abandoned `pyvcf`; Pangolin's imports say `import vcf` and pyvcf3 provides that).
- `pysam`, `pyfaidx`, `gffutils`, `biopython`, `pyfastx`, `scikit-learn`, `pandas`.

Bring the env up in one shot:

```
python3.12 -m venv .venv-predictors
source .venv-predictors/bin/activate
pip install "setuptools<81"
pip install spliceai==1.3.1 tensorflow==2.21.0 pysam pyfaidx pandas numpy scipy scikit-learn
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install pyvcf3 gffutils biopython pyfastx
git clone --depth 1 https://github.com/tkzeng/Pangolin.git /tmp/Pangolin
pip install /tmp/Pangolin
```

## Reference data

Both predictors are run against the same reference FASTA:

- `data/reference/GRCh38_chr7.fa` — chr7 from Ensembl release 110 (SHA-256 pinned in `src/harvest/download_reference.py`). Chromosome names are bare (`7`, not `chr7`), matching SpliceAI's example convention.

Pangolin additionally needs a `gffutils` DB derived from GENCODE:

- `data/reference/gencode/gencode.v44.chr7.ensembl.gtf.gz` — chr7-only slice of GENCODE v44 primary annotation, with `chr` prefixes stripped so contigs align with the FASTA.
- `data/reference/gencode/gencode.v44.chr7.ensembl.db` — gffutils DB built with `--filter Ensembl_canonical`, matching Pangolin's default annotation strategy.

Both are rebuildable with:

```
python -m src.harvest.download_reference
python -m src.harvest.download_gencode
```

## Compatibility shims

Two shims live under `src/predict/`. Both are idempotent and safe to apply repeatedly. Both are enforced by our runners as a subprocess preamble so the shim fires before the predictor's module is imported.

### `spliceai_shims.py`

SpliceAI 1.3.1's `one_hot_encode` calls `numpy.fromstring(seq, numpy.int8)`. NumPy 2 keeps that entry point as a stub that raises `ValueError` in binary mode. The shim probes for that failure and, if present, replaces `numpy.fromstring` with a `numpy.frombuffer` wrapper that also handles SpliceAI's `str`+latin-1 encoding path.

### `pangolin_shims.py`

Pangolin targets the unmaintained `pyvcf` 0.6.8, whose `vcf.parser._Info` namedtuple has 6 fields. The maintained `pyvcf3` fork added a 7th `type_code` field. The shim wraps `_Info` so that six-arg calls (Pangolin's) still work by appending `None`.

## Output schema

Both `run_spliceai.py` and `run_pangolin.py` write to `results/predictions/<predictor>.tsv` with a shared schema:

```
variant_id    score    <predictor-native fields ...>    symbol
```

- `variant_id` matches `data/variant_list_frozen.tsv`.
- `score` is in `[0, 1]`, higher = more likely splice-affecting. For SpliceAI this is `max(DS_AG, DS_AL, DS_DG, DS_DL)`. For Pangolin it is the max absolute value across gain and loss tokens for the annotated gene(s).
- `symbol` is the gene symbol used for pre-registered filtering (`CFTR`).

The analysis pipeline (`src.analyse.run_primary_analysis`, `src.analyse.run_stratum1`) reads these files via `--score-column score`.

The resulting scored TSVs are hash-locked in `results/predictions/RESULTS_MANIFEST.tsv`. Any change to the pipeline that produces different scores must be recorded there.

## AlphaGenome

AlphaGenome scoring is deliberately deferred. The AlphaGenome licence request to Google DeepMind is currently in a courtesy-notice window; scoring will begin once either a permission response arrives, an explicit restriction triggers the pre-registered expressivity-only fallback, or the 14-day window closes without response and the pre-registered interpretation (non-commercial academic use permitted) takes effect. The cron `ad481e61` fires on 2026-09-27 to make that decision.

## Reproducibility caveats

- Model weights: SpliceAI ships its four weight files inside the `spliceai` pypi package (frozen with the version pin). Pangolin ships weights in the git repository (frozen at the commit cloned above). Any change in either invalidates the recorded output hashes.
- The gffutils DB's binary encoding can vary across gffutils versions. `src/harvest/download_gencode.py` records a hash for the DB used at scoring time but treats a mismatch as a warning rather than a hard failure. Predictor outputs remain the ground truth for reproducibility auditing.
- Neither predictor uses random sampling at scoring time, so re-running the pipeline on the same inputs and weights produces byte-identical output TSVs.
