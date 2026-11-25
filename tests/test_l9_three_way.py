import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
import l9_three_way as l9  # noqa: E402


def row(wh, lp, asm):
    return [0, False, 0, wh, lp, asm]


def test_three_way_needs_all_three_sources_to_agree():
    # cis everywhere: every arm calls cis
    a, b = row([0, 1], [0, 7], [0, 3]), row([0, 1], [0, 7], [0, 3])
    assert [l9.call(a, b, arm) for arm in l9.ARMS] == [1, 1, 1]
    # LongPhase says trans: whatshap+asm still calls, three-way does not
    a, b = row([0, 1], [0, 7], [0, 3]), row([0, 1], [1, 7], [0, 3])
    assert [l9.call(a, b, arm) for arm in l9.ARMS] == [1, None, None]
    # assembly crosses contigs: nothing is called
    a, b = row([0, 1], [0, 7], [0, 3]), row([0, 1], [0, 7], [0, 4])
    assert [l9.call(a, b, arm) for arm in l9.ARMS] == [None, None, None]
    # both phasers agree on trans with the assembly
    a, b = row([0, 1], [1, 7], [0, 3]), row([1, 1], [0, 7], [1, 3])
    assert l9.call(a, b, "three-way") == 0
