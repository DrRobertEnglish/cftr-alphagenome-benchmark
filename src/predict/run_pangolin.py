"""Run Pangolin on the frozen scoring VCF.

Applies a pyvcf3 compatibility shim, invokes the pangolin CLI as a subprocess
against a chr7-scoped GENCODE gffutils DB, then parses the produced VCF into
our shared predictor output schema:

    variant_id \t score \t pangolin_gene \t pangolin_raw \t symbol

Where `score` = maximum absolute Pangolin gain/loss across the annotated CFTR
gene, in [0, 1], higher = more likely splice-affecting (matches the ranking
convention used by run_spliceai.py).

Pre-registration: OSF DOI 10.17605/OSF.IO/6PGX8.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_VCF = REPO_ROOT / "data" / "variants_for_scoring.vcf"
DEFAULT_FASTA = REPO_ROOT / "data" / "reference" / "GRCh38_chr7.fa"
DEFAULT_DB = REPO_ROOT / "data" / "reference" / "gencode" / "gencode.v44.chr7.ensembl.db"
DEFAULT_OUTPUT = REPO_ROOT / "results" / "predictions" / "pangolin.tsv"


PREAMBLE = r"""
import sys
from src.predict.pangolin_shims import apply_shims
apply_shims()
from pangolin.pangolin import main
sys.argv[0] = 'pangolin'
sys.exit(main())
"""


def run_pangolin_cli(
    vcf_path: Path,
    fasta_path: Path,
    db_path: Path,
    intermediate_prefix: Path,
) -> None:
    """Invoke Pangolin CLI with shims applied via subprocess preamble."""
    argv = [
        sys.executable,
        "-c",
        PREAMBLE,
        str(vcf_path),
        str(fasta_path),
        str(db_path),
        str(intermediate_prefix),
        "-d",
        "50",
        "-m",
        "True",
    ]
    subprocess.run(argv, check=True, cwd=str(REPO_ROOT))


# Pangolin INFO format:
#   gene_ensembl_id|pos:score_gain|pos:score_loss|warnings,...
# Multiple genes are comma-separated. Each per-gene block has a gain (positive
# delta at the top-scoring position) and a loss (negative delta at the
# top-scoring position). We take the max absolute value across genes as the
# ranking score.
_TOKEN_RE = re.compile(r"(-?\d+):(-?[\d.]+)")


def extract_max_abs_score(pangolin_info: str) -> tuple[float, str]:
    """Return (max_abs_score, comma-separated list of gene ids)."""
    if not pangolin_info or pangolin_info == ".":
        return 0.0, ""
    genes = []
    max_abs = 0.0
    for gene_block in pangolin_info.split(","):
        parts = gene_block.split("|")
        if not parts:
            continue
        gene_id = parts[0]
        genes.append(gene_id)
        for token in parts[1:]:
            m = _TOKEN_RE.fullmatch(token)
            if m:
                v = abs(float(m.group(2)))
                max_abs = max(max_abs, v)
    return max_abs, ",".join(genes)


def parse_pangolin_vcf(vcf_path: Path) -> list[dict]:
    """Parse the Pangolin-annotated VCF into row dicts."""
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
            pangolin = "."
            for kv in info.split(";"):
                if kv.startswith("Pangolin="):
                    pangolin = kv[len("Pangolin=") :]
                    break
            score, genes = extract_max_abs_score(pangolin)
            rows.append(
                {
                    "variant_id": variant_id,
                    "score": score,
                    "pangolin_gene": genes,
                    "pangolin_raw": pangolin,
                    "symbol": "CFTR",
                }
            )
    return rows


def write_tsv(rows: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    columns = ["variant_id", "score", "pangolin_gene", "pangolin_raw", "symbol"]
    with out_path.open("w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write(
                "\t".join(
                    [
                        str(row["variant_id"]),
                        f"{row['score']:.6g}",
                        row["pangolin_gene"],
                        row["pangolin_raw"],
                        row["symbol"],
                    ]
                )
                + "\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, default=DEFAULT_VCF)
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--intermediate-vcf",
        type=Path,
        default=None,
        help="Path prefix for the Pangolin CLI output VCF. Default: sibling of --output",
    )
    args = parser.parse_args()

    if args.intermediate_vcf is None:
        args.intermediate_vcf = args.output.with_suffix("").with_name(
            args.output.stem + "_raw"
        )

    args.intermediate_vcf.parent.mkdir(parents=True, exist_ok=True)

    print(f"Running Pangolin on {args.vcf}", flush=True)
    run_pangolin_cli(
        vcf_path=args.vcf,
        fasta_path=args.fasta,
        db_path=args.db,
        intermediate_prefix=args.intermediate_vcf,
    )

    produced_vcf = Path(str(args.intermediate_vcf) + ".vcf")
    if not produced_vcf.exists():
        raise SystemExit(f"Pangolin did not produce expected VCF at {produced_vcf}")

    rows = parse_pangolin_vcf(produced_vcf)
    scored = sum(1 for r in rows if r["score"] > 0 or r["pangolin_gene"])
    unscorable = len(rows) - scored
    write_tsv(rows, args.output)
    print(f"Wrote {args.output}")
    print(f"  Variants: {len(rows)} ({scored} scored, {unscorable} unscorable)")


if __name__ == "__main__":
    main()
