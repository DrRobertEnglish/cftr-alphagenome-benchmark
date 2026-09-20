"""Patch tkzeng/Pangolin to work with pyvcf3.

Pangolin (upstream commit c1e6c33) targets the unmaintained pyvcf 0.6.8 API.
pyvcf3 is the maintained fork; it exposes `import vcf` but its `_Info`
namedtuple has an extra `type_code` field. Rather than editing site-packages
in place, this shim patches the module at import time.
"""

from __future__ import annotations

_APPLIED = False


def apply_shims() -> None:
    """Idempotent monkey-patch for pangolin.pangolin under pyvcf3."""
    global _APPLIED
    if _APPLIED:
        return

    import vcf.parser as _vcf_parser

    _original_Info = _vcf_parser._Info

    def _Info_compat(*args, **kwargs):
        # If caller (Pangolin) passed 6 positional args, append a None
        # for the new type_code field required by pyvcf3.
        if len(args) == 6 and "type_code" not in kwargs:
            args = (*args, None)
        return _original_Info(*args, **kwargs)

    _vcf_parser._Info = _Info_compat  # type: ignore[assignment]
    _APPLIED = True
