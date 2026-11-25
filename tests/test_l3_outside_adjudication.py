import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))

from l3_outside_adjudication import corrected, phaser_share, wilson  # noqa: E402


def test_correction_inverts_label_noise():
    p, e = 0.8, 0.05
    f = p * (1 - e) + (1 - p) * e
    assert abs(corrected(f, e) - p) < 1e-12
    assert corrected(0.01, 0.05) == 0.0 and corrected(1.0, 0.05) == 1.0


def test_interval_brackets_point_and_is_conservative():
    p, lo, hi = phaser_share(900, 100, 10, 1000)
    assert lo < p < hi
    assert lo <= corrected(wilson(900, 1000)[0], 0.01) + 1e-12  # e at its lower bound is the low end
    assert phaser_share(0, 0, 0, 0)[1:] == (0.0, 1.0)
