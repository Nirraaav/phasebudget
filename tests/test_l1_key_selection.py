import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))

from l1_key_selection import keyed, matched_random, restrict  # noqa: E402


def test_key_restrict_and_random_control():
    keys = {"chr1": keyed([[10, 0, 0], [20, 1, 0], [30, None, None], [40, 1, 1]])}
    assert keys == {"chr1": {10, 40}}
    rows = [("segdup", ("chr1", 0), (1,), True, "chr1", 10, 40),
            ("segdup", ("chr1", 0), (1,), False, "chr1", 10, 20),
            ("easy", ("chr1", 0), (1,), False, "chr1", 30, 40)]
    kept = restrict(rows, keys)
    assert kept == rows[:1]
    rand = matched_random(rows, kept)
    assert all(r[0] == "segdup" for r in rand) and len(rand) <= 2
    # family key: relations compared, contig orientation free
    fam = {"chr1": {"asm": {10: [1, "c1"], 20: [0, "c1"], 30: [0, "c2"]},
                    "tra": {10: 0, 20: 1, 30: 0}}}
    pairs = [("easy", None, (1,), False, "chr1", a, b) for a, b in ((10, 20), (20, 30), (10, 40))]
    assert restrict(pairs, fam) == pairs[:1]
