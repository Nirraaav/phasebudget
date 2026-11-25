import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from l8_key_selection import link_class  # noqa: E402


def test_link_class():
    assert link_class((0, 0), (1, 1)) == "key"
    assert link_class((0, None), (1, 1)) == "missing"
    assert link_class((0, None), (1, None)) == "missing"
    assert link_class((0, 1), (1, 1)) == "disagree"
    assert link_class((0, None), (1, 0)) == "disagree"
