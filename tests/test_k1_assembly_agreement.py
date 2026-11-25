"""K1's agreement filter and K4's link rewrite on fixed and random examples."""

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))

import k1_assembly_agreement as k1  # noqa: E402


def test_agrees():
    assert k1.agrees(1, 1) and k1.agrees(0, 0)
    assert not k1.agrees(1, 0)
    assert not k1.agrees(1, None)  # the assembly does not relate the pair
    assert not k1.agrees(None, 1)


def test_rewrite_links_fixed():
    # pos1, allele, truth, score, indel; variant 2 unphased, so 1 links to 3.
    raw = [[10, 0, 0, [0.9, 5], False],
           [20, 1, 1, [0.8, 4], False],
           [25, None, 0, None, False],
           [30, 1, 1, [0.7, 3], False],
           [40, 0, 0, None, False]]
    assembly = {10: [1, "a|b"], 20: [0, "a|b"], 30: [0, "a|b"], 40: [1, "c|d"]}
    out = k1.rewrite_links(raw, assembly)
    assert out[0][3] == [0.9, 5]        # 10-20: phaser trans, assembly trans
    assert out[1][3] == [0.8, 4]        # 20-30 (25 unphased): both cis
    assert out[3][3] == list(k1.LOW)    # 30-40: the assembly's contigs differ
    assert out[4][3] is None            # no link to rewrite
    out = k1.rewrite_links(raw, {**assembly, 20: [1, "a|b"]})
    assert out[0][3] == list(k1.LOW)    # assembly now cis where the phaser is trans
    assert out[1][3] == list(k1.LOW)    # and trans where the phaser is cis
    assert raw[0][3] == [0.9, 5]        # input untouched


def brute(raw, assembly):
    out = [list(v) for v in raw]
    phased = [k for k, v in enumerate(raw) if v[1] is not None]
    for a, b in zip(phased, phased[1:]):
        if raw[a][3] is None:
            continue
        x, y = assembly.get(raw[a][0]), assembly.get(raw[b][0])
        same_block = x is not None and y is not None and x[1] == y[1]
        if not (same_block and (x[0] == y[0]) == (raw[a][1] == raw[b][1])):
            out[a][3] = list(k1.LOW)
    return out


def test_rewrite_links_random_against_brute_force():
    rng = random.Random(11)
    for _ in range(500):
        n = rng.randint(0, 8)
        raw = [[10 * (k + 1), rng.choice([0, 1, None]), rng.choice([0, 1, None]),
                rng.choice([None, [rng.random(), rng.randint(0, 9)]]), rng.random() < 0.2]
               for k in range(n)]
        assembly = {}
        for v in raw:
            if rng.random() < 0.8:
                assembly[v[0]] = [rng.randint(0, 1), rng.choice(["a|b", "c|d"])]
        assert k1.rewrite_links(raw, assembly) == brute(raw, assembly)
