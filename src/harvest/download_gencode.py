"""Download GENCODE v44 primary annotation and derive a chr7 gffutils DB.

Pangolin needs a gffutils database keyed on the same chromosome naming as the
reference FASTA. Our FASTA uses bare Ensembl-style contig names (`7`), while
GENCODE emits UCSC-style names (`chr7`), so we strip the prefix before building
the DB. The GENCODE release (v44) matches Ensembl 110, which matches the
release of the chr7 FASTA fetched by download_reference.py — see the CFTR
methodology doc for the frozen release choice.

Usage:
    python -m src.harvest.download_gencode
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = REPO_ROOT / "data" / "reference" / "gencode"

GENCODE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_44/gencode.v44.primary_assembly.annotation.gtf.gz"
)
# SHA-256 of the chr7-only, prefix-stripped GTF and the resulting gffutils DB.
# These pin the exact annotation Pangolin was scored against.
SHA256_CHR7_ENSEMBL_GTF = "08e10fe8ac713cb71dcf045012b44d94ae9bc7f9d4f02d0f81892e4b105014f0"
SHA256_CHR7_ENSEMBL_DB = "f482950ffe019f2a53947602f317b56fa6d08c362104a8765ef3d4f8e26c324e"


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    print(f"Downloading {url}", flush=True)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as resp, dest.open("wb") as out:
        shutil.copyfileobj(resp, out)


def slice_and_reprefix(full_gz: Path, out_gz: Path) -> None:
    """Extract chr7 rows and strip the `chr` prefix so contigs match FASTA."""
    print(f"Slicing to chr7 and reprefixing to Ensembl naming: {out_gz}", flush=True)
    with gzip.open(full_gz, "rt") as fin, gzip.open(out_gz, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.split("\t", 1)
            if fields[0] == "chr7":
                fout.write("7\t" + fields[1])


def build_db(gtf_gz: Path, db_path: Path) -> None:
    """Invoke tkzeng/Pangolin's create_db.py script.

    We defer to their script rather than reimplementing so the DB is
    schema-identical to what the pangolin CLI expects.
    """
    print(f"Building gffutils DB: {db_path}", flush=True)
    # tkzeng/Pangolin ships scripts/create_db.py; it may live in the checkout
    # left in /tmp, or on the user's PATH. Fall back to a Python search.
    candidates = [
        Path("/tmp/Pangolin/scripts/create_db.py"),
        REPO_ROOT / "third_party" / "Pangolin" / "scripts" / "create_db.py",
    ]
    script = next((p for p in candidates if p.exists()), None)
    if script is None:
        raise SystemExit(
            "create_db.py not found. Clone tkzeng/Pangolin and retry, "
            "or drop create_db.py into third_party/Pangolin/scripts/."
        )
    subprocess.run(
        [sys.executable, str(script), str(gtf_gz), "--filter", "Ensembl_canonical"],
        check=True,
    )
    # Script writes DB next to input, with .db suffix
    produced = gtf_gz.with_suffix("")  # strip .gz
    produced = produced.with_suffix(".db")  # replace .gtf with .db
    if produced != db_path:
        produced.rename(db_path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--keep-full",
        action="store_true",
        help="Retain the 48 MB full-genome GENCODE GTF (default: delete after slicing).",
    )
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    full_gz = OUT_DIR / "gencode.v44.primary_assembly.annotation.gtf.gz"
    chr7_gz = OUT_DIR / "gencode.v44.chr7.ensembl.gtf.gz"
    db_path = OUT_DIR / "gencode.v44.chr7.ensembl.db"

    if not chr7_gz.exists():
        if not full_gz.exists():
            download(GENCODE_URL, full_gz)
        slice_and_reprefix(full_gz, chr7_gz)
        if not args.keep_full:
            full_gz.unlink()

    got = sha256_of(chr7_gz)
    if got != SHA256_CHR7_ENSEMBL_GTF:
        raise SystemExit(
            f"SHA-256 mismatch for chr7 GTF.\n  expected: {SHA256_CHR7_ENSEMBL_GTF}\n  got:      {got}"
        )
    print(f"OK GTF: {chr7_gz} ({got[:12]}...)")

    if not db_path.exists():
        build_db(chr7_gz, db_path)

    got_db = sha256_of(db_path)
    if got_db != SHA256_CHR7_ENSEMBL_DB:
        # gffutils DBs can be non-deterministic across gffutils versions;
        # warn rather than fail
        print(
            f"WARN DB hash differs from recorded pin.\n"
            f"  expected: {SHA256_CHR7_ENSEMBL_DB}\n"
            f"  got:      {got_db}\n"
            f"  (gffutils DBs can vary across library versions; verify predictor outputs.)",
            file=sys.stderr,
        )
    else:
        print(f"OK DB:  {db_path} ({got_db[:12]}...)")


if __name__ == "__main__":
    main()
