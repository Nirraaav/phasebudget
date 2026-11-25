"""ctypes binding to the native dynamic program in `cpp/dp.cpp`.

The library is optional. `make -C cpp` builds `_pbdp.so` next to this file;
without it, or with `PHASEBUDGET_NATIVE=0` in the environment, `available()`
is false and `CoverageDP` runs the pure-Python pass. The native pass answers
to the Python one: `tests/test_native.py` requires equal costs, haplotypes,
assignments and boundary costs, not merely close ones.
"""

from __future__ import annotations

import ctypes
import os
from pathlib import Path
from typing import Optional, Sequence

_LIBRARY_PATH = Path(__file__).with_name("_pbdp.so")
_library: Optional[ctypes.CDLL] = None
_loaded = False

_int32_p = ctypes.POINTER(ctypes.c_int32)
_double_p = ctypes.POINTER(ctypes.c_double)
_int8_p = ctypes.POINTER(ctypes.c_int8)


def _load() -> Optional[ctypes.CDLL]:
    global _library, _loaded
    if _loaded:
        return _library
    _loaded = True
    if os.environ.get("PHASEBUDGET_NATIVE", "1") == "0" or not _LIBRARY_PATH.exists():
        return None
    library = ctypes.CDLL(str(_LIBRARY_PATH))
    library.pbdp_run.restype = ctypes.c_int
    library.pbdp_run.argtypes = [
        ctypes.c_int,
        _int32_p,
        _int32_p,
        _double_p,
        _double_p,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        _double_p,
        _int8_p,
        _int8_p,
        _double_p,
        _double_p,
    ]
    library.pbdp_posterior.restype = ctypes.c_int
    library.pbdp_posterior.argtypes = [
        ctypes.c_int,
        _int32_p,
        _int32_p,
        _double_p,
        _double_p,
        ctypes.c_int,
        ctypes.c_int,
        _double_p,
        _double_p,
    ]
    _library = library
    return _library


def available() -> bool:
    """True when the native library is built and not disabled."""
    return _load() is not None


def run(
    n_variants: int,
    offsets: Sequence[int],
    active: Sequence[int],
    c0: Sequence[float],
    c1: Sequence[float],
    n_fragments: int,
    *,
    want_margins: bool,
    interval: int = 0,
) -> tuple[float, tuple[int, ...], tuple[int, ...], Optional[list[float]], Optional[list[float]]]:
    """One native pass: `(cost, haplotype, assignment, zero, one)`.

    `offsets` delimits each column's slots in `active`, `c0` and `c1`, in the
    order `CoverageDP` assigns bit positions. `zero` and `one` are the
    boundary costs, or None when `want_margins` is false. `interval` is the
    checkpoint spacing; 0 means the square root of the variant count.
    """
    library = _load()
    if library is None:
        raise RuntimeError("the native dynamic program is not built; run `make -C cpp`")
    boundaries = max(n_variants - 1, 0)
    cost = ctypes.c_double(0.0)
    haplotype = (ctypes.c_int8 * max(n_variants, 1))()
    assignment = (ctypes.c_int8 * max(n_fragments, 1))()
    zero = (ctypes.c_double * max(boundaries, 1))() if want_margins else None
    one = (ctypes.c_double * max(boundaries, 1))() if want_margins else None
    status = library.pbdp_run(
        n_variants,
        (ctypes.c_int32 * max(len(offsets), 1))(*offsets),
        (ctypes.c_int32 * max(len(active), 1))(*active),
        (ctypes.c_double * max(len(c0), 1))(*c0),
        (ctypes.c_double * max(len(c1), 1))(*c1),
        n_fragments,
        1 if want_margins else 0,
        interval,
        ctypes.byref(cost),
        haplotype,
        assignment,
        ctypes.cast(zero, _double_p) if zero is not None else None,
        ctypes.cast(one, _double_p) if one is not None else None,
    )
    if status == 1:
        raise RuntimeError("no feasible state; the fragment matrix is malformed")
    if status == 2:
        raise RuntimeError("the native traceback broke")
    if status != 0:
        raise RuntimeError(f"the native dynamic program failed with status {status}")
    return (
        cost.value,
        tuple(haplotype[j] for j in range(n_variants)),
        tuple(assignment[f] for f in range(n_fragments)),
        list(zero[:boundaries]) if zero is not None else None,
        list(one[:boundaries]) if one is not None else None,
    )


def posterior(
    n_variants: int,
    offsets: Sequence[int],
    active: Sequence[int],
    c0: Sequence[float],
    c1: Sequence[float],
    n_fragments: int,
    *,
    interval: int = 0,
) -> list[tuple[float, float]]:
    """`(log Z_0(b), log Z_1(b))` for every boundary, by the native sum-product pass.

    Agrees with `phasebudget.posterior.boundary_log_sums` to rounding;
    `tests/test_native.py` checks it against that and against the oracle.
    """
    library = _load()
    if library is None:
        raise RuntimeError("the native dynamic program is not built; run `make -C cpp`")
    boundaries = max(n_variants - 1, 0)
    if boundaries == 0:
        return []
    zero = (ctypes.c_double * boundaries)()
    one = (ctypes.c_double * boundaries)()
    status = library.pbdp_posterior(
        n_variants,
        (ctypes.c_int32 * max(len(offsets), 1))(*offsets),
        (ctypes.c_int32 * max(len(active), 1))(*active),
        (ctypes.c_double * max(len(c0), 1))(*c0),
        (ctypes.c_double * max(len(c1), 1))(*c1),
        n_fragments,
        interval,
        zero,
        one,
    )
    if status != 0:
        raise RuntimeError(f"the native posterior pass failed with status {status}")
    return list(zip(zero, one))
