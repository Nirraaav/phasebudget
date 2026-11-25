import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))

import l6_platform_agreement as l6  # noqa: E402


def test_entries_and_summary():
    own = {"chr1": {"variants": [[10, "A", "G"], [20, "C", "T"], [30, "G", "A"], [40, "T", "C"]],
                    "whatshap": [[0, 1], [0, 1], [1, 1], [0, 1]]}}
    # The other platform: same relation at 10-20, opposite at 20-30, does not call 40.
    other = {"chr1": {"variants": [[10, "A", "G"], [20, "C", "T"], [30, "G", "A"]],
                      "whatshap": [[1, 5], [1, 5], [1, 5]]}}
    rows = [("chr1", 10, (1.0,), True, 20, None),   # wrong, other agrees: shared
            ("chr1", 20, (1.0,), True, 30, None),   # wrong, other opposite: caught
            ("chr1", 30, (1.0,), False, 40, None)]  # other does not call 40: unmatched
    got = l6.entries("P", rows, own, other, "whatshap", lambda c, a, b: "segdup")
    assert [r[4:] for r in got] == [(True, True, True), (False, True, True),
                                    (False, False, False)]
    assert got[0][1] == ("P", "chr1", 0)
    s = l6.summary(got)
    assert (s["wrong"], s["W"], s["S"], s["share"], s["agree"]) == (2, 2, 1, 0.5, 1)


def test_entries_indel_and_split_block():
    own = {"chr1": {"variants": [[10, "A", "AT"], [20, "C", "T"]],
                    "whatshap": [[0, 1], [0, 1]]}}
    other = {"chr1": {"variants": [[10, "A", "AT"], [20, "C", "T"]],
                      "whatshap": [[0, 1], [0, 2]]}}
    got = l6.entries("P", [("chr1", 10, (1.0,), False, 20, None)], own, other, "whatshap",
                     lambda c, a, b: "easy")
    assert got[0][0] == "indel" and got[0][4:] == (False, False, True)
