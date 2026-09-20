"""Assemble the CFTR AlphaGenome benchmark preprint.

Concatenates the five section drafts into a single preprint markdown file
with a standard preprint header, ordered sections, embedded figures, and
an aggregated References section built from every unique inline URL
across all sections. Inline links are preserved (readers still see the
markdown link) and are additionally numbered as [n] pointers into the
References list.
"""

from __future__ import annotations

import re
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MS = REPO / "docs" / "manuscript"
OUT = MS / "preprint_v0_1.md"

# Section order for the final preprint. Titles are rewritten to preprint-
# standard H1 headings; each source file's own top-level H1 is stripped
# because the preprint has one title and one H1 per section.
SECTIONS: list[tuple[str, str]] = [
    ("abstract_v0_1.md", "Abstract"),
    ("introduction_v0_1.md", "Introduction"),
    ("methods_v0_1.md", "Methods"),
    ("results_section_v0_1.md", "Results"),
    ("discussion_and_clinical_v0_1.md", "Discussion"),
]

TITLE = (
    "Ranking residual normal splicing in CFTR variants: a pre-registered "
    "benchmark of AlphaGenome against SpliceAI and Pangolin"
)
AUTHOR = "Robert English"
AFFILIATION = "Barts Health NHS Trust, London, United Kingdom"
CORRESPONDENCE = (
    "Correspondence: via the OSF project record "
    "([osf.io/5hvgf](https://osf.io/5hvgf/))."
)
PREPRINT_DATE = date(2026, 9, 20).strftime("%d %B %Y")
OSF_LINE = (
    "**Pre-registration.** OSF [osf.io/5hvgf](https://osf.io/5hvgf/), "
    "DOI [10.17605/OSF.IO/6PGX8](https://doi.org/10.17605/OSF.IO/6PGX8), "
    "registered 13 September 2026."
)
CODE_LINE = (
    "**Code and data.** "
    "[github.com/DrRobertEnglish/cftr-alphagenome-benchmark]"
    "(https://github.com/DrRobertEnglish/cftr-alphagenome-benchmark)."
)


# Markdown link matcher. URLs may contain balanced parentheses (e.g. the
# SpliceAI Cell paper URL contains "S0092-8674(18)31629-5.pdf"), so we
# accept any run of URL-safe characters plus balanced parens, and stop
# only at whitespace or an unmatched closing paren.
_URL_CHAR = r"[^\s()]"
_URL_BAL_PAREN = rf"\({_URL_CHAR}*\)"
MD_LINK_RE = re.compile(
    rf"\[([^\]]+)\]\((https?://(?:{_URL_CHAR}|{_URL_BAL_PAREN})+)\)"
)


def strip_leading_h1(body: str) -> str:
    """Remove a leading top-level heading from a section body."""
    lines = body.lstrip().splitlines()
    if lines and lines[0].startswith("# "):
        lines = lines[1:]
        while lines and not lines[0].strip():
            lines.pop(0)
    return "\n".join(lines).rstrip() + "\n"


def demote_headings(body: str) -> str:
    """Demote all headings by one level (# -> ##, ## -> ###, etc.)."""
    out_lines = []
    for line in body.splitlines():
        m = re.match(r"^(#{1,5})\s+(.*)$", line)
        if m:
            hashes, rest = m.groups()
            out_lines.append("#" + hashes + " " + rest)
        else:
            out_lines.append(line)
    return "\n".join(out_lines)


def collect_urls(body: str, registry: dict[str, int]) -> None:
    """Register every unique URL in the body into a shared numbering."""
    for _, url in MD_LINK_RE.findall(body):
        if url not in registry:
            registry[url] = len(registry) + 1


def annotate_links(body: str, registry: dict[str, int]) -> str:
    """Append [n] pointers after each markdown link (in-place)."""

    def repl(match: re.Match[str]) -> str:
        anchor, url = match.groups()
        n = registry[url]
        return f"[{anchor}]({url}) [{n}]"

    return MD_LINK_RE.sub(repl, body)


def sanitise_first_paragraph_of(section_name: str, body: str) -> str:
    """Small per-section fixes that don't belong in the source files."""
    if section_name == "Results":
        # The source draft has a parenthetical repo/commit note we don't
        # want twice (the preprint has its own Code and data line).
        body = re.sub(
            r"; predictor outputs, analysis code, and hash-locked inputs are archived in the accompanying repository \(`DrRobertEnglish/cftr-alphagenome-benchmark`, commit `[a-f0-9]+`\)\.",
            ".",
            body,
        )
    return body


def build_figures_block() -> str:
    """Return the figures block, inserted between Results and Discussion."""
    return (
        "## Figures\n\n"
        "![Figure 1 — Stratum 2 ROC curves](../../figures/figure1_roc_stratum2.png)\n\n"
        "*Figure 1.* Receiver-operating-characteristic (ROC) curves on "
        "Stratum 2 (n = 52; 43 splice-affecting, 9 no-effect) for "
        "AlphaGenome, SpliceAI 1.3.1, and Pangolin. Solid lines show the "
        "point-estimate ROC; shaded bands are 10,000-iteration "
        "stratified-bootstrap 95% confidence intervals on the true-positive "
        "rate at each false-positive rate. AUCs and their bootstrap 95% CIs "
        "are annotated in the legend. Dashed diagonal: chance performance.\n\n"
        "![Figure 2 — Stratum 1 predictor score vs residual splicing]"
        "(../../figures/figure2_stratum1_scatter.png)\n\n"
        "*Figure 2.* Predictor score versus fraction of normally spliced "
        "transcript on Stratum 1 (n = 8) for AlphaGenome (left, blue), "
        "SpliceAI (centre, red), and Pangolin (right, green). Solid lines "
        "show a monotone rank-based fit for visual reference; reported "
        "Spearman ρ and 10,000-iteration bootstrap 95% confidence intervals "
        "are annotated in each panel title. The pre-specified expected "
        "direction (negative ρ: higher predictor score tracking lower "
        "normal-transcript fraction) is met by all three predictors.\n"
    )


def main() -> None:
    registry: dict[str, int] = {}

    # First pass: collect URLs across every section, in source order.
    prepared: list[tuple[str, str]] = []
    for filename, section_title in SECTIONS:
        raw = (MS / filename).read_text()
        body = strip_leading_h1(raw)
        body = sanitise_first_paragraph_of(section_title, body)
        body = demote_headings(body)
        collect_urls(body, registry)
        prepared.append((section_title, body))

    # Also collect URLs from the constant blocks so they get [n] tags too.
    for constant in (OSF_LINE, CODE_LINE, CORRESPONDENCE):
        collect_urls(constant, registry)

    # Second pass: annotate every link with its [n] pointer.
    annotated_sections: list[tuple[str, str]] = [
        (title, annotate_links(body, registry)) for title, body in prepared
    ]

    parts: list[str] = []
    parts.append(f"# {TITLE}\n")
    parts.append(f"**{AUTHOR}**  \n{AFFILIATION}\n")
    parts.append(f"*Preprint version 0.1 — {PREPRINT_DATE}*\n")
    parts.append(annotate_links(CORRESPONDENCE, registry))
    parts.append(annotate_links(OSF_LINE, registry))
    parts.append(annotate_links(CODE_LINE, registry))
    parts.append("\n---\n")

    figures_inserted = False
    for title, body in annotated_sections:
        parts.append(f"# {title}\n")
        parts.append(body)
        parts.append("\n---\n")
        # Insert the figures block between Results and Discussion so that
        # readers see the figures before the Discussion refers to them.
        if title == "Results" and not figures_inserted:
            parts.append(build_figures_block())
            parts.append("\n---\n")
            figures_inserted = True

    # References
    parts.append("# References\n")
    for url, n in sorted(registry.items(), key=lambda kv: kv[1]):
        parts.append(f"{n}. {url}")
    parts.append("")

    OUT.write_text("\n".join(parts).rstrip() + "\n")
    print(f"Wrote {OUT} ({OUT.stat().st_size:,} bytes, {len(registry)} unique references)")


if __name__ == "__main__":
    main()
