"""Build a scoring-ready VCF from the frozen variant list.

Reads data/variant_list_frozen.tsv, looks up each variant's REF nucleotide in
the local GRCh38 chr7 FASTA (data/reference/GRCh38_chr7.fa), validates that
the REF in the frozen TSV matches the reference base at that position, and
writes data/variants_for_scoring.vcf in VCFv4.2 format compatible with
SpliceAI (-A grch38) and Pangolin.

Hard-fails on any REF mismatch. If a mismatch appears, the frozen TSV has a
data-integrity problem (strand slip, coordinate conversion error, or VariantValidator
inconsistency) and must be investigated before scoring.

VCF conventions:
  - Chromosome column: "7" (no chr prefix), matches Ensembl FASTA naming.
  - VCFv4.2 header only; no sample columns.
  - INFO column carries variant_id and stratum as key=value pairs for
    downstream joining back to the frozen TSV.
  - Sorted by POS ascending.

Idempotent. Overwrites the output VCF on every run.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd
from pyfaidx import Fasta

VCF_HEADER_TEMPLATE = """##fileformat=VCFv4.2
##fileDate={date}
##source=cftr-alphagenome-benchmark/src.harvest.build_scoring_vcf
##reference=GRCh38 (Ensembl release 110, chr7 only)
##INFO=<ID=variant_id,Number=1,Type=String,Description="Variant ID from data/variant_list_frozen.tsv">
##INFO=<ID=stratum,Number=1,Type=String,Description="Study stratum (1=quantitative, 2=binary)">
##contig=<ID=7,length=159345973>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_and_build(
    variants_tsv: Path,
    fasta_path: Path,
    output_vcf: Path,
) -> dict:
    """Return dict of {n_variants, n_snv, n_indel, sha256}."""
    df = pd.read_csv(variants_tsv, sep="\t")
    required = {"variant_id", "chrom", "grch38_pos", "ref", "alt", "stratum"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"variant TSV missing required columns: {sorted(missing)}")

    fa = Fasta(str(fasta_path))
    if "7" not in fa:
        raise ValueError(f"FASTA {fasta_path} does not contain chromosome '7'")

    # Validate every REF
    mismatches: list[str] = []
    for _, row in df.iterrows():
        chrom = str(row["chrom"])
        if chrom != "7":
            mismatches.append(
                f"{row['variant_id']}: expected chrom=7, got {chrom!r}"
            )
            continue
        pos = int(row["grch38_pos"])
        ref = str(row["ref"]).upper()
        # pyfaidx is 0-based half-open; VCF POS is 1-based
        ref_len = len(ref)
        actual = fa["7"][pos - 1 : pos - 1 + ref_len].seq.upper()
        if actual != ref:
            mismatches.append(
                f"{row['variant_id']} at chr{chrom}:{pos}: expected REF={ref}, "
                f"FASTA has {actual}"
            )
    if mismatches:
        msg = "REF mismatches between frozen TSV and GRCh38 chr7 FASTA:\n  - " + \
              "\n  - ".join(mismatches)
        raise SystemExit(msg + f"\n\n{len(mismatches)} mismatch(es). Refusing to write VCF.")

    # Sort by POS ascending for VCF convention
    df_sorted = df.sort_values("grch38_pos", kind="mergesort").reset_index(drop=True)

    # Emit VCF
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    # Header is emitted with a deterministic fileDate (the OSF registration
    # date) so the VCF hash is reproducible run-to-run rather than drifting
    # with wall-clock time.
    lines: list[str] = [VCF_HEADER_TEMPLATE.format(date="20260913")]
    for _, row in df_sorted.iterrows():
        info = f"variant_id={row['variant_id']};stratum={row['stratum']}"
        lines.append(
            f"7\t{int(row['grch38_pos'])}\t{row['variant_id']}\t"
            f"{str(row['ref']).upper()}\t{str(row['alt']).upper()}\t.\t.\t{info}"
        )
    output_vcf.write_text("\n".join(lines) + "\n")

    n_snv = int(
        ((df_sorted["ref"].str.len() == 1) & (df_sorted["alt"].str.len() == 1)).sum()
    )
    n_indel = len(df_sorted) - n_snv
    return {
        "n_variants": len(df_sorted),
        "n_snv": n_snv,
        "n_indel": n_indel,
        "sha256": sha256_of(output_vcf),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=Path,
                        default=Path("data/variant_list_frozen.tsv"))
    parser.add_argument("--fasta", type=Path,
                        default=Path("data/reference/GRCh38_chr7.fa"))
    parser.add_argument("--output", type=Path,
                        default=Path("data/variants_for_scoring.vcf"))
    args = parser.parse_args()

    if not args.fasta.exists():
        raise SystemExit(
            f"FASTA {args.fasta} not found. Run: "
            f"python -m src.harvest.download_reference"
        )

    summary = validate_and_build(args.variants, args.fasta, args.output)
    print(f"Wrote {args.output}")
    print(f"  Variants: {summary['n_variants']} ({summary['n_snv']} SNV, {summary['n_indel']} indel)")
    print(f"  SHA-256: {summary['sha256']}")
    print("All REF nucleotides validated against GRCh38 chr7 FASTA.")


if __name__ == "__main__":
    main()
