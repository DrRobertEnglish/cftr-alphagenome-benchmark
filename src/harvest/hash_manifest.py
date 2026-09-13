"""
Hash-manifest generator and verifier for frozen inputs.

The frozen variant list, the reference exon annotation, and any other file
under `data/` that a predictor batch reads must have a SHA-256 hash recorded
in `data/HASH_MANIFEST.tsv`. Any change to those files after the OSF
pre-registration date requires a versioned OSF amendment before further
predictor queries — the CI workflow refuses to pass if the recorded hashes
and on-disk hashes disagree.

Usage:
    # Regenerate the manifest (only allowed when explicitly amending):
    python -m src.harvest.hash_manifest --write

    # Verify (CI mode — non-zero exit on mismatch):
    python -m src.harvest.hash_manifest --check
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
MANIFEST_PATH = REPO_ROOT / "data" / "HASH_MANIFEST.tsv"

# Files whose hashes are tracked. Order matters for reproducibility of the
# manifest file itself. Paths are relative to REPO_ROOT. Add entries here as
# the study progresses; the CI check enforces them all.
TRACKED_FILES: list[str] = [
    "data/reference/CFTR_ENST00000003084_11_exons_GRCh38.tsv",
    # "data/variant_list_frozen.tsv",   # add on freeze
]


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def read_manifest() -> dict[str, str]:
    if not MANIFEST_PATH.exists():
        return {}
    m: dict[str, str] = {}
    with MANIFEST_PATH.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise ValueError(f"malformed manifest line: {line!r}")
            if parts == ["path", "sha256"]:
                continue  # header row
            m[parts[0]] = parts[1]
    return m


def write_manifest(entries: dict[str, str]) -> None:
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST_PATH.open("w") as f:
        f.write("# SHA-256 manifest of frozen inputs.\n")
        f.write("# Rebuilding this file after the OSF registration date requires\n")
        f.write("# a documented, dated amendment in docs/protocol/AMENDMENTS.md.\n")
        f.write("path\tsha256\n")
        for path in sorted(entries):
            f.write(f"{path}\t{entries[path]}\n")


def compute_all() -> dict[str, str]:
    out: dict[str, str] = {}
    for rel in TRACKED_FILES:
        p = REPO_ROOT / rel
        if not p.exists():
            raise FileNotFoundError(f"tracked file missing: {rel}")
        out[rel] = sha256_of(p)
    return out


def cmd_write() -> int:
    entries = compute_all()
    write_manifest(entries)
    for path, h in sorted(entries.items()):
        print(f"{path}\t{h}")
    print(f"\nwrote {len(entries)} entries to {MANIFEST_PATH.relative_to(REPO_ROOT)}")
    return 0


def cmd_check() -> int:
    recorded = read_manifest()
    live = compute_all()

    missing_from_manifest = sorted(set(live) - set(recorded))
    missing_from_disk = sorted(set(recorded) - set(live))
    mismatched = sorted(
        p for p in set(live) & set(recorded) if recorded[p] != live[p]
    )

    if not (missing_from_manifest or missing_from_disk or mismatched):
        print(f"OK: {len(live)} tracked file(s) match recorded hashes")
        return 0

    if missing_from_manifest:
        print("MISSING FROM MANIFEST:")
        for p in missing_from_manifest:
            print(f"  {p} (on disk: {live[p]})")
    if missing_from_disk:
        print("MISSING FROM DISK:")
        for p in missing_from_disk:
            print(f"  {p} (recorded: {recorded[p]})")
    if mismatched:
        print("HASH MISMATCH:")
        for p in mismatched:
            print(f"  {p}\n    recorded: {recorded[p]}\n    on disk:  {live[p]}")
    return 2


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument("--write", action="store_true",
                   help="Regenerate HASH_MANIFEST.tsv from live file hashes.")
    g.add_argument("--check", action="store_true",
                   help="Verify manifest matches live hashes; non-zero exit on mismatch.")
    args = parser.parse_args()
    return cmd_write() if args.write else cmd_check()


if __name__ == "__main__":
    sys.exit(main())
