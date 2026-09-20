"""Compatibility shims applied at import time so SpliceAI 1.3.1 runs on
current NumPy (>= 2) and Python (>= 3.12).

Upstream SpliceAI has not been updated for NumPy 2 or Python 3.12+. Rather
than pinning the whole scientific stack backwards, we apply two narrow
patches at import time:

1. ``np.fromstring`` — removed in NumPy 2. Re-attached as a small wrapper
   over ``np.frombuffer`` with the same signature SpliceAI needs
   (``np.fromstring(seq, np.int8)`` where ``seq`` is a str).

2. ``pkg_resources`` availability — SpliceAI's ``__init__`` calls
   ``pkg_resources.get_distribution(__name__).version``. If setuptools is
   present but pkg_resources has been split out, we defer to the caller to
   ensure setuptools is installed. Nothing to patch here.

Import ``apply_shims()`` **before** importing anything from ``spliceai``.
"""
from __future__ import annotations

import numpy as np

_SHIMS_APPLIED = False


def apply_shims() -> None:
    """Idempotent. Apply all compatibility shims required for SpliceAI 1.3.1."""
    global _SHIMS_APPLIED
    if _SHIMS_APPLIED:
        return

    # NumPy 2 keeps np.fromstring as a stub that raises ValueError in binary
    # mode. Probe and unconditionally replace when the binary path is broken.
    _needs_shim = False
    try:
        np.fromstring(b"\x00", np.int8)
    except ValueError:
        _needs_shim = True
    except Exception:  # noqa: BLE001
        _needs_shim = True
    if _needs_shim:
        def _fromstring(string, dtype=float, count=-1, sep=""):
            if sep:
                # legacy text-mode: split and parse
                return np.array([dtype(s) for s in string.split(sep)])
            if isinstance(string, str):
                # SpliceAI's actual call site uses str + latin-1 bytes 0x00-0x04
                string = string.encode("latin-1")
            if count == -1:
                return np.frombuffer(string, dtype=dtype)
            return np.frombuffer(string, dtype=dtype, count=count)
        np.fromstring = _fromstring  # type: ignore[attr-defined]

    _SHIMS_APPLIED = True
