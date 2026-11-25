"""The floor of `benchmarks/r2_batch_labels.py` (R2) against `learn_then_test`."""

import importlib.util
import math
import sys
from pathlib import Path

from phasebudget.confidence import Units, learn_then_test

BENCH = Path(__file__).resolve().parent.parent / "benchmarks"


def _module():
    sys.path.insert(0, str(BENCH))
    spec = importlib.util.spec_from_file_location("r2_batch_labels", BENCH / "r2_batch_labels.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_floor_is_where_all_correct_units_start_to_certify():
    m = _module()
    beta, delta = 1e-3, 0.05
    floor = m.floor_units(beta, delta)
    assert 2990 < floor < 3000
    for n, certifies in ((math.floor(floor), False), (math.ceil(floor), True)):
        units = Units.from_pairs([(u, (1.0,), False) for u in range(n)])
        assert (learn_then_test(units, [(1.0,)], beta, delta) is not None) == certifies


def test_checked_never_goes_looser_than_the_pool():
    m = _module()
    rows = [("easy", u, (s,), False) for u in range(5000) for s in (0.2, 0.9)]
    free = m.by_group(rows, 1e-3, 0.05)["easy"]
    capped = m.by_group(rows, 1e-3, 0.05, {"easy": (0.9,)})["easy"]
    assert free < (0.9,) and capped >= (0.9,)
    assert m.by_group(rows, 1e-3, 0.05, {"easy": None})["easy"] is None
