"""Download the GRCh38 chr7 primary-assembly FASTA from Ensembl release 110.

Every variant in data/variant_list_frozen.tsv is on chr7, so we cache only
that chromosome (~46 MB gzipped, ~155 MB unzipped) rather than the full 3 GB
primary assembly. The predictor runners (src/predict/run_spliceai.py,
src/predict/run_pangolin.py) point at this file for REF lookup and flanking
context.

Idempotent. If the file already exists AND its SHA-256 matches the recorded
hash in data/HASH_MANIFEST.tsv, this script exits cleanly without
re-downloading.

Source URL (frozen at Ensembl release 110):
  https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/
    Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz

Chromosome naming: `7` (no `chr` prefix). SpliceAI and Pangolin both accept
this naming; the VCF builder writes matching `7` in the CHROM column.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import shutil
import sys
import urllib.request
from pathlib import Path

ENSEMBL_URL = (
    "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/"
    "Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz"
)
EXPECTED_SHA256_FA = "71b9f38d87ecb8db7067498602e5ae44bc09a6080e276c11a24a3d2933b1ed93"
EXPECTED_SHA256_GZ = "d20985eabdf27638b81067c94189acc3e330730ba23f77c8ffa15576a590ad4a"


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(dest_dir: Path, force: bool = False) -> None:
    dest_dir.mkdir(parents=True, exist_ok=True)
    gz_path = dest_dir / "GRCh38_chr7.fa.gz"
    fa_path = dest_dir / "GRCh38_chr7.fa"

    if fa_path.exists() and not force:
        got = sha256_of(fa_path)
        if got == EXPECTED_SHA256_FA:
            print(f"OK: {fa_path} already present and matches expected SHA-256.")
            return
        print(f"WARNING: {fa_path} exists but hash mismatch (got {got})")
        print("Re-downloading.")

    print(f"Downloading {ENSEMBL_URL} -> {gz_path}")
    with urllib.request.urlopen(ENSEMBL_URL) as resp, gz_path.open("wb") as out:
        shutil.copyfileobj(resp, out)

    got_gz = sha256_of(gz_path)
    if got_gz != EXPECTED_SHA256_GZ:
        raise SystemExit(
            f"Downloaded gz SHA-256 mismatch. Expected {EXPECTED_SHA256_GZ}, got {got_gz}. "
            "Ensembl may have republished release 110; investigate before proceeding."
        )

    print(f"Decompressing -> {fa_path}")
    with gzip.open(gz_path, "rb") as src, fa_path.open("wb") as dst:
        shutil.copyfileobj(src, dst)

    got_fa = sha256_of(fa_path)
    if got_fa != EXPECTED_SHA256_FA:
        raise SystemExit(
            f"Decompressed FASTA SHA-256 mismatch. Expected {EXPECTED_SHA256_FA}, got {got_fa}."
        )
    print(f"OK: {fa_path} downloaded and hash-verified ({got_fa[:16]}...)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dest-dir", type=Path,
        default=Path("data/reference"),
        help="Where to place GRCh38_chr7.fa[.gz] (default: data/reference)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-download even if the file already exists with a matching hash",
    )
    args = parser.parse_args()
    try:
        download(args.dest_dir, force=args.force)
    except (urllib.error.URLError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
