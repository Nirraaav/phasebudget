import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
import l7_score as l7  # noqa: E402


def test_outside_stats_needs_both_variants_inside():
    def inside(chrom, pos):
        return 100 <= pos < 200

    rows = [("chr1", 110, 150, False), ("chr1", 120, 160, True),   # inside: 1 wrong of 2
            ("chr1", 150, 250, True), ("chr1", 300, 400, True),    # outside: 3 wrong of 4
            ("chr1", 300, 410, False), ("chr1", 50, 120, True)]
    s = l7.outside_stats(rows, inside)
    assert s["pairs"] == {"inside": [2, 1], "outside": [4, 3]}
    assert s["share_outside"] == 0.75
    assert s["ratio"] == 1.5
