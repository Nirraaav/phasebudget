import itertools
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))

import l2_end_calls as l2  # noqa: E402


def test_end_relation_equals_chain_parity():
    rng = random.Random(1)
    for _ in range(200):
        alleles = [rng.randint(0, 1) for _ in range(8)]
        values = [[a, 7] for a in alleles]
        for i, j in itertools.combinations(range(8), 2):
            parity = 1
            for k in range(i, j):  # product of the link relations along the chain
                parity ^= 1 - int(alleles[k] == alleles[k + 1])
            assert l2.relation(values[i], values[j]) == parity
    assert l2.relation([0, 1], [0, 2]) is None and l2.relation([0, None], [0, None]) is None


def test_read_agreement():
    assert l2.read_agrees(1, (2, 0)) and l2.read_agrees(0, (1, 2))
    assert not l2.read_agrees(1, (1, 0)) and not l2.read_agrees(1, (1, 1))
    assert not l2.read_agrees(None, (5, 0)) and not l2.read_agrees(1, None)


def test_certified_floor():
    assert l2.certified(299, 0, 1e-2) and not l2.certified(298, 0, 1e-2)
    assert l2.certified(2995, 0, 1e-3) and not l2.certified(2990, 0, 1e-3)


def test_pairs_brute_force():
    rng = random.Random(2)
    for _ in range(50):
        variants = sorted({(rng.randint(1, 300_000),) for _ in range(40)})
        genes = [(s, s + rng.randint(1, 150_000), f"g{k}")
                 for k, s in enumerate(rng.sample(range(1, 300_000), 4))]
        got = {(i, j): sorted(n) for i, j, n in l2.pairs(variants, genes)}
        want = {}
        for i, j in itertools.combinations(range(len(variants)), 2):
            pi, pj = variants[i][0], variants[j][0]
            names = sorted(g for s, e, g in genes if s <= pi and pj <= e)
            if names and pj - pi <= l2.WINDOW:
                want[(i, j)] = names
        assert got == want
