"""Build the frozen machine-readable CFTR variant TSV from the OSF-registered
variant harvest (docs/protocol/variant_harvest_v1.md).

Sources of truth:

- Joynt et al. 2020, PLOS Genetics 16(10):e1009100 — Tables 1 (exonic, n=15)
  and 2 (intronic, n=37). Full 52-variant list constitutes Stratum 2.
- Masvidal et al. 2014, EJHG — quantitative % normally spliced for c.580-1G>T,
  c.2657+5G>A, c.[2657+5G>A; 2562T>G] complex.
- Chiba-Falek et al. 1999, AJRCCM — quantitative % normally spliced for
  c.3717+12191C>T (3849+10kbC>T) across multiple tissues.
- Deletang et al. 2022, Gene Therapy — quantitative % for c.1680-886A>G.

GRCh38 coordinates are resolved via the VariantValidator public REST API
(https://rest.variantvalidator.org/) which anchors HGVS c. notation on
the CFTR MANE Select transcript NM_000492.4 (== ENST00000003084.11).

Output: data/variant_list_frozen.tsv with columns:
    variant_id | source | stratum | hgvs_c | legacy_name |
    chrom | grch38_pos | ref | alt |
    binary_splice_affecting | pct_normally_spliced | pct_sd | notes
"""
from __future__ import annotations

import csv
import json
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
OUT_TSV = REPO_ROOT / "data" / "variant_list_frozen.tsv"

VV_TRANSCRIPT = "NM_000492.4"  # CFTR MANE Select, matches ENST00000003084.11
VV_BASE = "https://rest.variantvalidator.org/VariantValidator/variantvalidator/GRCh38"


@dataclass
class Variant:
    variant_id: str  # short unique key
    source: str  # primary paper
    stratum: str  # "1" (quantitative) or "2" (binary from Joynt)
    hgvs_c: str
    legacy_name: str
    binary_splice_affecting: int  # 1 == splices differently vs WT, 0 == no splice effect
    pct_normally_spliced: float | None = None  # Stratum 1 only
    pct_sd: float | None = None  # standard deviation or SEM as reported
    notes: str = ""
    chrom: str = ""
    grch38_pos: int = 0
    ref: str = ""
    alt: str = ""
    resolved_hgvs_g: str = ""

    def as_row(self) -> dict:
        return {
            "variant_id": self.variant_id,
            "source": self.source,
            "stratum": self.stratum,
            "hgvs_c": self.hgvs_c,
            "legacy_name": self.legacy_name,
            "chrom": self.chrom,
            "grch38_pos": self.grch38_pos,
            "ref": self.ref,
            "alt": self.alt,
            "binary_splice_affecting": self.binary_splice_affecting,
            "pct_normally_spliced": ""
            if self.pct_normally_spliced is None
            else f"{self.pct_normally_spliced:.4f}",
            "pct_sd": "" if self.pct_sd is None else f"{self.pct_sd:.4f}",
            "resolved_hgvs_g": self.resolved_hgvs_g,
            "notes": self.notes,
        }


# ---------------------------------------------------------------------------
# JOYNT 2020 — TABLE 1 (exonic, n=15)
# Splice-affecting flag uses the "RNA effect" column verbatim from the paper:
# "Missplices" == 1, "Does not missplice" == 0.
# ---------------------------------------------------------------------------
JOYNT_EXONIC: list[Variant] = [
    Variant("Joynt_E01", "Joynt2020", "2", "c.166G>A",  "E56K",   0),
    Variant("Joynt_E02", "Joynt2020", "2", "c.274G>A",  "E92K",   1, notes="Missplices per RNA-effect column; primary isoform normally spliced"),
    Variant("Joynt_E03", "Joynt2020", "2", "c.580G>A",  "G194R",  0),
    Variant("Joynt_E04", "Joynt2020", "2", "c.581G>T",  "G194V",  0),
    Variant("Joynt_E05", "Joynt2020", "2", "c.2909G>A", "G970D",  0),
    Variant("Joynt_E06", "Joynt2020", "2", "c.3719T>G", "V1240G", 0),
    Variant("Joynt_E07", "Joynt2020", "2", "c.454A>G",  "M152V",  1),
    Variant("Joynt_E08", "Joynt2020", "2", "c.523A>G",  "I175V",  1),
    Variant("Joynt_E09", "Joynt2020", "2", "c.2816A>G", "H939R",  1),
    Variant("Joynt_E10", "Joynt2020", "2", "c.3700A>G", "I1234V", 1),
    Variant("Joynt_E11", "Joynt2020", "2", "c.2908G>C", "G970R",  1),
    Variant("Joynt_E12", "Joynt2020", "2", "c.2908G>A", "G970S",  1),
    Variant("Joynt_E13", "Joynt2020", "2", "c.3717G>C", "R1239S", 1),
    Variant("Joynt_E14", "Joynt2020", "2", "c.3872A>G", "Q1291R", 0),
    Variant("Joynt_E15", "Joynt2020", "2", "c.3873G>C", "Q1291H", 1),
]

# ---------------------------------------------------------------------------
# JOYNT 2020 — TABLE 2 (intronic, n=37)
# Splice-affecting flag from the "RNA effect" column verbatim.
# ---------------------------------------------------------------------------
JOYNT_INTRONIC: list[Variant] = [
    # Canonical splice-site subset (n=19, all Missplice)
    Variant("Joynt_I01", "Joynt2020", "2", "c.164+1G>A", "296+1G>A", 1),
    Variant("Joynt_I02", "Joynt2020", "2", "c.164+2T>C", "296+2T>C", 1),
    Variant("Joynt_I03", "Joynt2020", "2", "c.165-2A>G", "297-2A>G", 1),
    Variant("Joynt_I04", "Joynt2020", "2", "c.273+1G>A", "405+1G>A", 1),
    Variant("Joynt_I05", "Joynt2020", "2", "c.274-1G>A", "406-1G>A", 1),
    Variant("Joynt_I06", "Joynt2020", "2", "c.274-2A>G", "406-2A>G", 1),
    Variant("Joynt_I07", "Joynt2020", "2", "c.489+1G>T", "621+1G>T", 1),
    Variant("Joynt_I08", "Joynt2020", "2", "c.579+1G>T", "711+1G>T", 1),
    Variant("Joynt_I09", "Joynt2020", "2", "c.1584+1G>A", "1716+1G>A", 1),
    Variant("Joynt_I10", "Joynt2020", "2", "c.1585-1G>A", "1717-1G>A", 1),
    Variant("Joynt_I11", "Joynt2020", "2", "c.2658-1G>C", "2790-1G>C", 1),
    Variant("Joynt_I12", "Joynt2020", "2", "c.2658-2A>G", "2790-2A>G", 1),
    Variant("Joynt_I13", "Joynt2020", "2", "c.2988+1G>A", "3120+1G>A", 1),
    Variant("Joynt_I14", "Joynt2020", "2", "c.3469-2A>G", "3601-2A>G", 1),
    Variant("Joynt_I15", "Joynt2020", "2", "c.3717+1G>A", "3849+1G>A", 1),
    Variant("Joynt_I16", "Joynt2020", "2", "c.3718-1G>A", "3850-1G>A", 1),
    Variant("Joynt_I17", "Joynt2020", "2", "c.3873+1G>A", "4005+1G>A", 1),
    Variant("Joynt_I18", "Joynt2020", "2", "c.3873+2T>C", "4005+2T>C", 1, notes="Residual full-length transcript; canonical GT->GC donor"),
    Variant("Joynt_I19", "Joynt2020", "2", "c.4242+2T>C", "4374+2T>C", 1, notes="Residual full-length transcript; canonical GT->GC donor"),

    # Proximal intronic subset (n=15)
    Variant("Joynt_I20", "Joynt2020", "2", "c.164+3_164+4insT", "296+3insT", 0, notes="Not CF-causing per Joynt Table 2 footnote"),
    Variant("Joynt_I21", "Joynt2020", "2", "c.165-3C>T", "297-3C>T", 1),
    Variant("Joynt_I22", "Joynt2020", "2", "c.273+3A>C", "405+3A>C", 1),
    Variant("Joynt_I23", "Joynt2020", "2", "c.489+3A>G", "621+3A>G", 1),
    Variant("Joynt_I24", "Joynt2020", "2", "c.579+3A>G", "711+3A>G", 1),
    Variant("Joynt_I25", "Joynt2020", "2", "c.579+3A>C", "711+3A>C", 1),
    Variant("Joynt_I26", "Joynt2020", "2", "c.579+3A>T", "711+3A>T", 1),
    Variant("Joynt_I27", "Joynt2020", "2", "c.579+5G>A", "711+5G>A", 1),
    Variant("Joynt_I28", "Joynt2020", "2", "c.2657+2_2657+3insA", "2789+2insA", 1),
    Variant("Joynt_I29", "Joynt2020", "2", "c.2657+5G>A", "2789+5G>A", 1),
    Variant("Joynt_I30", "Joynt2020", "2", "c.3468+2_3468+3insT", "3600+2insT", 1),
    Variant("Joynt_I31", "Joynt2020", "2", "c.3468+5G>A", "3600+5G>A", 1),
    Variant("Joynt_I32", "Joynt2020", "2", "c.3717+4A>G", "3849+4A>G", 1),
    Variant("Joynt_I33", "Joynt2020", "2", "c.3717+5G>A", "3849+5G>A", 1),
    Variant("Joynt_I34", "Joynt2020", "2", "c.3718-3T>G", "3850-3T>G", 1),

    # Distal intronic subset (n=3)
    Variant("Joynt_I35", "Joynt2020", "2", "c.164+28A>G", "296+28A>G", 0, notes="Not CF-causing per Joynt Table 2 footnote"),
    Variant("Joynt_I36", "Joynt2020", "2", "c.2620-26A>G", "2752-26A>G", 0, notes="Not CF-causing per Joynt Table 2 footnote"),
    Variant("Joynt_I37", "Joynt2020", "2", "c.3717+40A>G", "3849+40A>G", 1),
]

# ---------------------------------------------------------------------------
# STRATUM 1 (n=9 unique variants + WT control + multi-tissue rows)
# Quantitative % normally spliced from the original papers.
#
# NOTE ON DOUBLE-USE: three of the Stratum 1 variants also appear in Joynt
# Table 2 (c.2657+5G>A, and the Joynt-quantified c.2908G>C, c.2908G>A,
# c.2909G>A, c.3873G>C are already covered in the exonic subset).
# We keep them here with the Stratum 1 quantitative measurement attached so
# they can be recovered by the H3 rank-correlation analysis; the Stratum 2
# analysis will de-duplicate by hgvs_c.
# ---------------------------------------------------------------------------
STRATUM1: list[Variant] = [
    Variant(
        "S1_01_Masvidal", "Masvidal2014", "1", "c.580-1G>T", "712-1G>T",
        1, pct_normally_spliced=40.0, pct_sd=1.3,
        notes="Nasal epithelium RT-qPCR (patient); carrier value 66%+/-2.7% not used here",
    ),
    Variant(
        "S1_02_Masvidal", "Masvidal2014", "1", "c.2657+5G>A", "2789+5G>A",
        1, pct_normally_spliced=24.0, pct_sd=9.0,
        notes="Nasal epithelium RT-qPCR (patient); minigene 29%+/-2.3% not used here",
    ),
    # c.[2657+5G>A;2562T>G] complex allele from Masvidal 2014 minigene
    # (23%+/-0.4% normally spliced) is EXCLUDED from the frozen TSV because
    # single-variant predictors (SpliceAI, Pangolin, AlphaGenome) cannot
    # score two cis-linked SNVs as one input. It is retained in the
    # narrative harvest for reference but not analysed.
    Variant(
        "S1_04_Joynt", "Joynt2020", "1", "c.2908G>C", "G970R",
        1, pct_normally_spliced=0.0, pct_sd=None,
        notes="0% full-length; 78.5%+/-3.2% skipped, 21.5%+/-3.2% 177nt-del",
    ),
    Variant(
        "S1_05_Joynt", "Joynt2020", "1", "c.2908G>A", "G970S",
        1, pct_normally_spliced=0.0, pct_sd=None,
        notes="0% full-length; 31.1%+/-2.7% skipped, 68.3%+/-2.7% 177nt-del",
    ),
    Variant(
        "S1_06_Joynt_WT", "Joynt2020", "1", "c.2909G>A", "G970D",
        0, pct_normally_spliced=100.0, pct_sd=None,
        notes="Wild-type-splicing control per Joynt Table 1",
    ),
    Variant(
        "S1_07_Joynt", "Joynt2020", "1", "c.3873G>C", "Q1291H",
        1, pct_normally_spliced=37.2, pct_sd=2.2,
    ),
    Variant(
        "S1_08_ChibaFalek_lung", "ChibaFalek1999", "1",
        "c.3717+12191C>T", "3849+10kbC>T",
        1, pct_normally_spliced=26.0, pct_sd=1.5,
        notes="Fetal lung RT-PCR",
    ),
    Variant(
        "S1_09_Deletang", "Deletang2022", "1", "c.1680-886A>G", "1811+1.6kbA>G",
        1, pct_normally_spliced=2.0, pct_sd=None,
        notes="Patient mRNA; reported range 1-3%, midpoint used",
    ),
]


def _vv_fetch(hgvs_c: str, retries: int = 3, delay: float = 1.5) -> dict:
    """Call VariantValidator for a single HGVS c. and return the parsed JSON."""
    query = f"{VV_TRANSCRIPT}:{hgvs_c}"
    url = f"{VV_BASE}/{urllib.parse.quote(query, safe=':>')}/{VV_TRANSCRIPT}"
    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"accept": "application/json"})
            with urllib.request.urlopen(req, timeout=45) as fh:
                return json.loads(fh.read().decode("utf-8"))
        except Exception as exc:  # noqa: BLE001 — API returns many failure modes
            last_err = exc
            time.sleep(delay * (attempt + 1))
    raise RuntimeError(f"VariantValidator failed for {hgvs_c}: {last_err}")


def resolve_grch38(variant: Variant) -> None:
    """Populate chrom/pos/ref/alt on the variant via VariantValidator."""
    data = _vv_fetch(variant.hgvs_c)
    key = f"{VV_TRANSCRIPT}:{variant.hgvs_c}"
    # VV sometimes normalises the HGVS string (e.g. spacing) — search keys
    match_key = key if key in data else next(
        (k for k in data if isinstance(data.get(k), dict) and "primary_assembly_loci" in data[k]),
        None,
    )
    if match_key is None:
        raise RuntimeError(f"No primary_assembly_loci for {variant.hgvs_c}: keys={list(data.keys())[:5]}")
    hit = data[match_key]
    grch38 = hit.get("primary_assembly_loci", {}).get("grch38", {})
    vcf = grch38.get("vcf")
    if not vcf:
        raise RuntimeError(f"No GRCh38 VCF block for {variant.hgvs_c}")
    variant.chrom = vcf["chr"]
    variant.grch38_pos = int(vcf["pos"])
    variant.ref = vcf["ref"]
    variant.alt = vcf["alt"]
    variant.resolved_hgvs_g = grch38.get("hgvs_genomic_description", "")


def main() -> None:
    variants: list[Variant] = [*JOYNT_EXONIC, *JOYNT_INTRONIC, *STRATUM1]
    print(f"Resolving GRCh38 coordinates for {len(variants)} variants via VariantValidator...")
    failures: list[tuple[str, str]] = []
    for i, v in enumerate(variants, 1):
        try:
            resolve_grch38(v)
            print(f"  [{i:2d}/{len(variants)}] {v.variant_id:22s} {v.hgvs_c:32s} -> chr{v.chrom}:{v.grch38_pos} {v.ref}>{v.alt}")
        except Exception as exc:  # noqa: BLE001
            failures.append((v.variant_id, str(exc)))
            print(f"  [{i:2d}/{len(variants)}] {v.variant_id:22s} {v.hgvs_c:32s} FAILED: {exc}")
        time.sleep(0.4)  # be polite to the public API

    if failures:
        print(f"\n{len(failures)} variant(s) failed resolution:")
        for vid, msg in failures:
            print(f"  {vid}: {msg}")
        raise SystemExit(1)

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(variants[0].as_row().keys())
    with OUT_TSV.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for v in variants:
            writer.writerow(v.as_row())
    print(f"\nWrote {OUT_TSV.relative_to(REPO_ROOT)} with {len(variants)} rows.")


if __name__ == "__main__":
    main()
