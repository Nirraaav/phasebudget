import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")  # benchmark dependency, not a library one

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
import l5_unit_independence as l5  # noqa: E402


def test_joins_and_null():
    # chr1 units 0..9 without 7, wrong at 3, 4, 5: two joins over 7 adjacent pairs.
    flags = l5.unit_flags([(("chr1", i), i in (3, 4, 5)) for i in range(10) if i != 7])
    obs, null, adj, rho = l5.joins(flags, 4000, np.random.default_rng(0))
    assert (obs, adj) == (2, 7)
    # Permutation null: 3 wrong of 9 units, E = 7 * (3 * 2) / (9 * 8).
    assert abs(null.mean() - 7 * 6 / 72) < 0.05
    assert rho > 0
    r = l5.test(flags)
    assert r["units"] == 9 and r["wrong"] == 3 and r["p"] < 0.2
    spread = l5.unit_flags([(("chr1", i), i in (0, 4, 8)) for i in range(10)])
    assert l5.joins(spread, 10, np.random.default_rng(0))[0] == 0


def test_widen_flat_and_nested_units():
    from l5_unit_independence import widen

    rows = [("segdup", ("HG002", "chr1", 5), (1.0,), False),
            ("segdup", ("HG002", ("chr1", 5)), (1.0,), True)]
    assert [u for _, u, _, _ in widen(rows, 2)] == [("HG002", "chr1", 2), ("HG002", ("chr1", 2))]
