"""Score the 60 frozen CFTR variants with AlphaGenome via the official SDK.

Reads the frozen variant list from `data/variant_list_frozen.tsv`, calls
`dna_client.score_variant` for each variant with the three splicing scorers
(splice sites, splice site usage, splice junctions), and computes the paper's
recommended merged splicing score:

    alphagenome_splicing = max(splice_sites) + max(splice_site_usage)
                                             + max(splice_junctions) / 5.0

For each scorer, the aggregation across genes and tracks is `max(abs(raw))`, per
the AlphaGenome docs. The output TSV `results/predictions/alphagenome.tsv`
carries one row per variant with columns:

    variant_id, chrom, pos, ref, alt,
    splice_sites, splice_site_usage, splice_junctions, score, error

`score` is the merged splicing score (higher = larger predicted splicing effect,
sign-blind). `error` is empty on success and holds the exception string on
failure. Failures are logged but do not abort the run; the pre-registered
analysis will drop unscorable variants at the loading step.

The API key is read from `ALPHA_GENOME_API_KEY` (loaded from `.env` if present,
which is gitignored). This runner writes the results TSV atomically, computes
its SHA-256 hash, and updates `results/predictions/RESULTS_MANIFEST.tsv`.
"""

from __future__ import annotations

import csv
import hashlib
import os
import sys
import time
from pathlib import Path
from typing import Any


def _load_env(dotenv_path: Path) -> None:
    if not dotenv_path.exists():
        return
    for line in dotenv_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            continue
        key, val = line.split("=", 1)
        key = key.strip()
        val = val.strip().strip('"').strip("'")
        if key and key not in os.environ:
            os.environ[key] = val


REPO_ROOT = Path(__file__).resolve().parents[2]
VARIANTS_TSV = REPO_ROOT / "data" / "variant_list_frozen.tsv"
OUT_TSV = REPO_ROOT / "results" / "predictions" / "alphagenome.tsv"
MANIFEST_TSV = REPO_ROOT / "results" / "predictions" / "RESULTS_MANIFEST.tsv"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def _update_manifest(target: Path) -> None:
    """Update or insert the target file's SHA-256 in RESULTS_MANIFEST.tsv."""
    manifest_lines: list[list[str]] = []
    header = ["path", "sha256", "bytes"]
    if MANIFEST_TSV.exists():
        with MANIFEST_TSV.open() as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if not row:
                    continue
                if row[0] == "path":
                    header = row
                    continue
                manifest_lines.append(row)
    rel = str(target.relative_to(REPO_ROOT))
    digest = _sha256(target)
    size = target.stat().st_size
    manifest_lines = [row for row in manifest_lines if row and row[0] != rel]
    manifest_lines.append([rel, digest, str(size)])
    manifest_lines.sort(key=lambda r: r[0])
    with MANIFEST_TSV.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(manifest_lines)


def _score_row(model: Any, tidy_scores_fn: Any, splicing_scorers: list, variant_row: dict, sequence_length: int) -> dict:
    from alphagenome.data import genome

    variant = genome.Variant(
        chromosome=f"chr{variant_row['chrom']}",
        position=int(variant_row["grch38_pos"]),
        reference_bases=variant_row["ref"],
        alternate_bases=variant_row["alt"],
    )
    interval = variant.reference_interval.resize(sequence_length)
    scores = model.score_variant(
        interval=interval,
        variant=variant,
        variant_scorers=splicing_scorers,
    )
    df = tidy_scores_fn([scores])
    row: dict[str, Any] = {
        "variant_id": variant_row["variant_id"],
        "chrom": variant_row["chrom"],
        "pos": variant_row["grch38_pos"],
        "ref": variant_row["ref"],
        "alt": variant_row["alt"],
        "splice_sites": "",
        "splice_site_usage": "",
        "splice_junctions": "",
        "score": "",
        "error": "",
    }
    if df is None or len(df) == 0:
        row["error"] = "empty scorer output"
        return row
    # Two column-name paths: output_type is documented on tidy_scores;
    # variant_scorer may contain the scorer identity. Try exact match on
    # each column, then fall back to case-insensitive equality on the
    # stringified value. aggregate by max(abs(raw_score)) per AlphaGenome docs.
    import re

    def _max_abs(output_type_name: str) -> float | None:
        pat = re.compile(rf"(^|[^A-Z_]){re.escape(output_type_name)}([^A-Z_]|$)")
        for col in ("output_type", "variant_scorer"):
            if col not in df.columns:
                continue
            col_str = df[col].astype(str)
            # Exact-uppercase match first
            sub = df[col_str.str.upper() == output_type_name]
            if len(sub) == 0:
                # Word-boundary regex match (avoids SPLICE_SITES matching SPLICE_SITE_USAGE)
                sub = df[col_str.str.upper().str.contains(pat, regex=True, na=False)]
            if len(sub) > 0 and "raw_score" in sub.columns:
                return float(sub["raw_score"].abs().max())
        return None

    ss = _max_abs("SPLICE_SITES")
    ssu = _max_abs("SPLICE_SITE_USAGE")
    sj = _max_abs("SPLICE_JUNCTIONS")
    row["splice_sites"] = f"{ss:.6f}" if ss is not None else ""
    row["splice_site_usage"] = f"{ssu:.6f}" if ssu is not None else ""
    row["splice_junctions"] = f"{sj:.6f}" if sj is not None else ""
    if ss is None or ssu is None or sj is None:
        row["error"] = f"missing scorer output (ss={ss}, ssu={ssu}, sj={sj})"
        return row
    merged = ss + ssu + sj / 5.0
    row["score"] = f"{merged:.6f}"
    return row


def main() -> int:
    _load_env(REPO_ROOT / ".env")
    api_key = os.environ.get("ALPHA_GENOME_API_KEY", "").strip()
    if not api_key:
        print(
            "ERROR: ALPHA_GENOME_API_KEY not set. Populate .env at the repo root.",
            file=sys.stderr,
        )
        return 2

    # Imports are deferred so a missing SDK yields a clear error.
    from alphagenome.models import dna_client, variant_scorers as vs_module

    splicing_scorers = [
        vs_module.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"],
        vs_module.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITE_USAGE"],
        vs_module.RECOMMENDED_VARIANT_SCORERS["SPLICE_JUNCTIONS"],
    ]

    print(f"[{time.strftime('%H:%M:%S')}] Connecting to AlphaGenome...")
    model = dna_client.create(api_key, timeout=60.0)
    print(f"[{time.strftime('%H:%M:%S')}] Connected. Loading frozen variants from {VARIANTS_TSV}")

    with VARIANTS_TSV.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
    print(f"[{time.strftime('%H:%M:%S')}] {len(rows)} variants to score")

    sequence_length = dna_client.SEQUENCE_LENGTH_1MB
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    tmp_out = OUT_TSV.with_suffix(".tsv.tmp")

    fieldnames = [
        "variant_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "splice_sites",
        "splice_site_usage",
        "splice_junctions",
        "score",
        "error",
    ]
    with tmp_out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for i, variant_row in enumerate(rows, 1):
            t0 = time.time()
            try:
                out = _score_row(model, vs_module.tidy_scores, splicing_scorers, variant_row, sequence_length)
            except Exception as e:
                out = {k: "" for k in fieldnames}
                out.update(
                    {
                        "variant_id": variant_row["variant_id"],
                        "chrom": variant_row["chrom"],
                        "pos": variant_row["grch38_pos"],
                        "ref": variant_row["ref"],
                        "alt": variant_row["alt"],
                        "error": f"{type(e).__name__}: {str(e)[:200]}",
                    }
                )
            writer.writerow(out)
            f.flush()
            dt = time.time() - t0
            status = "OK" if not out["error"] else f"ERR {out['error'][:80]}"
            print(f"[{time.strftime('%H:%M:%S')}] {i:3d}/{len(rows)} {variant_row['variant_id']:12s} {dt:5.1f}s {status}")

    tmp_out.replace(OUT_TSV)
    _update_manifest(OUT_TSV)
    digest = _sha256(OUT_TSV)
    print(f"\n[{time.strftime('%H:%M:%S')}] Wrote {OUT_TSV} ({OUT_TSV.stat().st_size} bytes)")
    print(f"[{time.strftime('%H:%M:%S')}] SHA-256: {digest}")
    print(f"[{time.strftime('%H:%M:%S')}] Manifest updated at {MANIFEST_TSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
