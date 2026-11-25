"""Liftover, SNV matching and relations of `benchmarks/w4_gate.py` (W4)."""

import importlib.util
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parent.parent / "benchmarks"


def _module():
    sys.path.insert(0, str(BENCH))
    spec = importlib.util.spec_from_file_location("w4_gate", BENCH / "w4_gate.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CHAIN = """chain 100 chrA 1000 + 100 230 chrA 900 + 50 175 1
10 5 0
20 0 3
30

chain 50 chrA 1000 + 500 520 chrA 900 - 100 120 2
20

chain 50 chrA 1000 + 600 610 chrB 900 + 0 10 3
10
"""


def test_lift_follows_blocks_and_gaps(tmp_path):
    m = _module()
    path = tmp_path / "x.chain"
    path.write_text(CHAIN)
    blocks = m.load_chain(path, "chrA")
    starts = [b[0] for b in blocks]
    # blocks: t[100,110)->q50, t[115,135)->q60, t[135,165)->q83
    assert m.lift(blocks, starts, 100) == (50, "+")
    assert m.lift(blocks, starts, 109) == (59, "+")
    assert m.lift(blocks, starts, 112) is None  # target gap
    assert m.lift(blocks, starts, 115) == (60, "+")
    assert m.lift(blocks, starts, 135) == (83, "+")  # after a query gap of 3
    assert m.lift(blocks, starts, 164) == (112, "+")
    assert m.lift(blocks, starts, 165) is None
    # the minus-strand chain: query 100 on the reverse strand is 900 - 1 - 100 forward
    assert m.lift(blocks, starts, 500) == (799, "-")
    assert m.lift(blocks, starts, 519) == (780, "-")
    assert m.lift(blocks, starts, 605) is None  # that chain goes to chrB


def test_match_needs_the_same_snv_alleles():
    m = _module()
    assert m.match("A", "G", "+", ("A", "G")) is False
    assert m.match("A", "G", "+", ("G", "A")) is True
    assert m.match("A", "G", "-", ("T", "C")) is False
    assert m.match("A", "G", "-", ("C", "T")) is True
    assert m.match("A", "G", "+", ("A", "T")) is None
    assert m.match("AT", "A", "+", ("AT", "A")) is None


def _table(alleles, truth, pairs):
    return {"alleles": alleles, "truth": truth,
            "next": {a: b for a, b, _, _ in pairs},
            "support": {a: (g, d) for a, _, g, d in pairs}}


def test_relation_walks_the_block_and_needs_every_pair_confident():
    m = _module()
    # variants 0..4; 2 unphased by the truth; block break between 3 and 4
    t = _table([0, 1, 1, 0, 0], [0, 1, None, 1, 0],
               [(0, 1, 12, 0), (1, 2, 15, 0), (2, 3, 9, 0)])
    assert m.relation(t, 0, 1) == (False, True)
    assert m.relation(t, 0, 3) == (True, False)  # phaser 0/0 same, truth 0/1 differ; 9 reads
    assert m.relation(t, 1, 3) == (True, False)  # phaser 1/0 differ, truth 1/1 same
    assert m.relation(t, 0, 2) is None  # truth does not phase 2
    assert m.relation(t, 3, 4) is None  # different blocks
    t2 = _table([0, 1], [0, 0], [(0, 1, 30, 1)])
    assert m.relation(t2, 0, 1) == (True, False)  # one disagreeing read: not confident
