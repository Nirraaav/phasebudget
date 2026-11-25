"""S4's error rule in `benchmarks/stratified_score.py`."""

import importlib.util
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parent.parent / "benchmarks"


def _module():
    sys.path.insert(0, str(BENCH))
    spec = importlib.util.spec_from_file_location("stratified_score", BENCH / "stratified_score.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_fewer_errors_halves_many_and_bounds_few():
    m = _module()
    # (units, wrong units, p, pairs, trusted, wrong, wrong trusted, certified)
    assert m.fewer_errors((0, 0, 1, 1000, 990, 40, 20, True))
    assert not m.fewer_errors((0, 0, 1, 1000, 990, 40, 21, True))
    assert m.fewer_errors((0, 0, 1, 671984, 671984, 1, 1, True))  # 1.5e-6 per trusted pair
    assert not m.fewer_errors((0, 0, 1, 5000, 5000, 2, 2, True))  # 4e-4 per trusted pair
