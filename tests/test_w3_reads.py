"""The window of `benchmarks/w3_reads.py` (W3b)."""

import importlib.util
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parent.parent / "benchmarks"


def _module():
    sys.path.insert(0, str(BENCH))
    spec = importlib.util.spec_from_file_location("w3_reads", BENCH / "w3_reads.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_window_is_three_before_the_pair_and_three_after():
    m = _module()
    phased = [0, 2, 3, 5, 7, 8, 9, 11, 12]
    blocks = {v: "A" for v in phased}
    # pair (5, 7) at index 3: slots for indices 0 to 7
    assert m.window(phased, 3, blocks) == [0, 2, 3, 5, 7, 8, 9, 11]


def test_window_leaves_slots_past_the_ends_and_in_other_blocks_empty():
    m = _module()
    phased = [0, 2, 3, 5, 7]
    blocks = {0: "A", 2: "B", 3: "B", 5: "B", 7: "C"}
    # pair (2, 3): slots for indices -2, -1, 0, 1, 2, 3, 4, 5
    assert m.window(phased, 1, blocks) == [-1, -1, -1, 2, 3, 5, -1, -1]
