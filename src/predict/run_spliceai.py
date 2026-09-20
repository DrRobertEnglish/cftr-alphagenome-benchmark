"""Run SpliceAI 1.3.1 against the frozen scoring VCF and emit per-variant
scores in the schema src.analyse.run_primary_analysis expects.

Handles two upstream-SpliceAI issues at import time (see spliceai_shims.py):
  1. NumPy 2 removed np.fromstring — shimmed to np.frombuffer.
  2. pkg_resources availability requires setuptools < 81 in the venv.

Output schema (results/predictions/spliceai.tsv):
    variant_id \\t score \\t ds_ag \\t ds_al \\t ds_dg \\t ds_dl \\t symbol

`score` is the maximum of the four SpliceAI delta scores (DS_AG, DS_AL, DS_DG,
DS_DL), which is the standard summary score used for pathogenicity ranking in
the SpliceAI paper (Jaganathan et al., Cell 2019). Downstream ROC-AUC on this
score is the pre-registered predictor value used in H1 and H2.

Score assumptions the analyse pipeline relies on:
  - Higher score = more likely to affect splicing (positive direction of
    association with binary_splice_affecting=1).
  - Values are in [0, 1].
  - Every frozen variant scored; unscored variants (e.g. deletions longer than
    2*D bp) are written with score=NaN so run_primary_analysis explicitly
    reports them under `n_missing_a` / `n_missing_b` rather than silently
    dropping them.

Usage:

    python -m src.predict.run_spliceai \\
        --vcf data/variants_for_scoring.vcf \\
        --fasta data/reference/GRCh38_chr7.fa \\
        --output results/predictions/spliceai.tsv

Runs offline. No API calls. TensorFlow CPU backend (~1-2 min per variant).
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

# Apply NumPy-2 compatibility shim BEFORE any spliceai import happens
from src.predict.spliceai_shims import apply_shims

apply_shims()

import pandas as pd


def run_spliceai_cli(
    input_vcf: Path,
    output_vcf: Path,
    fasta: Path,
    distance: int = 50,
) -> None:
    """Invoke the SpliceAI CLI as a subprocess, with the shim applied in the
    child via -c preamble."""
    # The shim must be applied in the child before spliceai imports.
    # Cleanest: run the CLI via a preamble that applies the shim.
    preamble = (
        f"import sys; sys.path.insert(0, {str(Path.cwd())!r}); "
        "from src.predict.spliceai_shims import apply_shims; apply_shims(); "
        "from spliceai.__main__ import main; "
        f"sys.argv = ['spliceai', '-I', {str(input_vcf)!r}, '-O', {str(output_vcf)!r}, "
        f"'-R', {str(fasta)!r}, '-A', 'grch38', '-D', {str(distance)!r}]; "
        "sys.exit(main())"
    )
    proc = subprocess.run(
        [sys.executable, "-c", preamble],
        check=False,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise SystemExit(
            f"SpliceAI failed (rc={proc.returncode}):\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    print(proc.stderr, file=sys.stderr)


def parse_spliceai_vcf(vcf_path: Path) -> pd.DataFrame:
    """Extract per-variant SpliceAI scores from the annotated VCF.

    SpliceAI's INFO adds SpliceAI=<ALLELE>|<SYMBOL>|<DS_AG>|<DS_AL>|<DS_DG>|
    <DS_DL>|<DP_AG>|<DP_AL>|<DP_DG>|<DP_DL>. A variant with no gene at the
    position gets no SpliceAI= entry; write those out with NaN scores.
    """
    rows: list[dict] = []
    with vcf_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            variant_id = fields[2]
            info = fields[7]
            spliceai_field = None
            for kv in info.split(";"):
                if kv.startswith("SpliceAI="):
                    spliceai_field = kv[len("SpliceAI="):]
                    break
            if spliceai_field is None or spliceai_field in (".", ""):
                rows.append({
                    "variant_id": variant_id,
                    "score": float("nan"),
                    "ds_ag": float("nan"), "ds_al": float("nan"),
                    "ds_dg": float("nan"), "ds_dl": float("nan"),
                    "symbol": "",
                })
                continue
            # If SpliceAI returns multiple alts (comma-separated), take the
            # first — our VCF is single-alt per row by construction.
            first = spliceai_field.split(",")[0]
            parts = first.split("|")
            if len(parts) < 6:
                rows.append({
                    "variant_id": variant_id, "score": float("nan"),
                    "ds_ag": float("nan"), "ds_al": float("nan"),
                    "ds_dg": float("nan"), "ds_dl": float("nan"),
                    "symbol": "",
                })
                continue
            _allele, symbol, ds_ag_s, ds_al_s, ds_dg_s, ds_dl_s = parts[:6]
            try:
                ds_ag = float(ds_ag_s); ds_al = float(ds_al_s)
                ds_dg = float(ds_dg_s); ds_dl = float(ds_dl_s)
                score = max(ds_ag, ds_al, ds_dg, ds_dl)
            except ValueError:
                ds_ag = ds_al = ds_dg = ds_dl = float("nan")
                score = float("nan")
            rows.append({
                "variant_id": variant_id, "score": score,
                "ds_ag": ds_ag, "ds_al": ds_al, "ds_dg": ds_dg, "ds_dl": ds_dl,
                "symbol": symbol,
            })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path,
                        default=Path("data/variants_for_scoring.vcf"))
    parser.add_argument("--fasta", type=Path,
                        default=Path("data/reference/GRCh38_chr7.fa"))
    parser.add_argument("--output", type=Path,
                        default=Path("results/predictions/spliceai.tsv"))
    parser.add_argument("--distance", type=int, default=50,
                        help="SpliceAI -D parameter (default 50)")
    parser.add_argument("--intermediate-vcf", type=Path,
                        default=Path("results/predictions/spliceai_raw.vcf"))
    args = parser.parse_args()

    if not args.fasta.exists():
        raise SystemExit(
            f"FASTA {args.fasta} not found. Run: "
            f"python -m src.harvest.download_reference"
        )
    if not args.vcf.exists():
        raise SystemExit(
            f"VCF {args.vcf} not found. Run: "
            f"python -m src.harvest.build_scoring_vcf"
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.intermediate_vcf.parent.mkdir(parents=True, exist_ok=True)

    print(f"Running SpliceAI on {args.vcf}")
    run_spliceai_cli(args.vcf, args.intermediate_vcf, args.fasta,
                     distance=args.distance)

    df = parse_spliceai_vcf(args.intermediate_vcf)
    df.to_csv(args.output, sep="\t", index=False)
    n_scored = int(df["score"].notna().sum())
    n_missing = len(df) - n_scored
    print(f"Wrote {args.output}")
    print(f"  Variants: {len(df)} ({n_scored} scored, {n_missing} unscorable)")


if __name__ == "__main__":
    main()
