"""The call, its certificate and the pair enumerator of `benchmarks/q4_clinical_cis_trans.py`."""

import importlib.util
import itertools
import random
import sys
from pathlib import Path

import pytest

BENCH = Path(__file__).resolve().parent.parent / "benchmarks"


def _module():
    sys.path.insert(0, str(BENCH))
    spec = importlib.util.spec_from_file_location("q4_clinical_cis_trans",
                                                  BENCH / "q4_clinical_cis_trans.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


q4 = _module()


def test_b_max_fixed_values():
    betas = q4.BETAS  # 1e-4, 3e-4, 1e-3, 3e-3, 1e-2
    assert betas[q4.b_max(1, 1e-2)] == 1e-2
    assert betas[q4.b_max(2, 1e-2)] == 3e-3
    assert betas[q4.b_max(3, 1e-2)] == 3e-3
    assert betas[q4.b_max(4, 1e-2)] == 1e-3
    assert betas[q4.b_max(10, 1e-2)] == 1e-3
    assert betas[q4.b_max(11, 1e-2)] == 3e-4
    assert betas[q4.b_max(34, 1e-2)] == 1e-4
    assert q4.b_max(101, 1e-2) is None
    assert betas[q4.b_max(1, 1e-3)] == 1e-3
    assert betas[q4.b_max(3, 1e-3)] == 3e-4
    assert betas[q4.b_max(4, 1e-3)] == 1e-4
    assert q4.b_max(11, 1e-3) is None


def test_link_level_and_nesting():
    t = [None, (0.9, 10), (0.8, 5), (0.8, 5), (0.5, 0)]
    q4.check_nested(t)
    assert q4.link_level((0.95, 3), t) == 1
    assert q4.link_level((0.8, 5), t) == 2
    assert q4.link_level((0.6, 9), t) == 4
    assert q4.link_level((0.4, 9), t) == 5
    with pytest.raises(ValueError):
        q4.check_nested([(0.5, 0), None])
    with pytest.raises(ValueError):
        q4.check_nested([(0.5, 0), (0.6, 0)])


def test_chain_parity_and_union_bound_by_enumeration():
    """A call is wrong exactly when an odd number of its links is wrong, and its
    probability of being wrong is at most the sum of its cells' risks, whatever the
    wrong links inside a bad cell are."""
    rng = random.Random(4)
    for _ in range(200):
        n_links = rng.randint(1, 6)
        cell_of = [rng.randint(0, 3) for _ in range(n_links)]
        cells = sorted(set(cell_of))
        risk = {c: rng.uniform(0, 0.3) for c in cells}
        phaser = [rng.randint(0, 1) for _ in range(n_links + 1)]
        worst = 0.0
        for bad in itertools.product((0, 1), repeat=len(cells)):
            p = 1.0
            for c, x in zip(cells, bad):
                p *= risk[c] if x else 1 - risk[c]
            allowed = [k for k in range(n_links) if bad[cells.index(cell_of[k])]]
            wrong_any = 0
            for flips in itertools.product((0, 1), repeat=len(allowed)):
                wrong_links = {k for k, f in zip(allowed, flips) if f}
                truth = [phaser[0]]
                for k in range(n_links):  # the truth's relation differs at a wrong link
                    same = phaser[k] == phaser[k + 1]
                    truth.append(truth[-1] if same != (k in wrong_links) else 1 - truth[-1])
                wrong = (phaser[0] ^ truth[0]) != (phaser[-1] ^ truth[-1])
                assert wrong == (len(wrong_links) % 2 == 1)
                wrong_any = max(wrong_any, wrong)
            worst += p * wrong_any
        assert worst <= sum(risk.values()) + 1e-12


def _random_chromosome(rng):
    n = rng.randint(0, 18)
    pos = sorted(rng.sample(range(1, 400), n))
    allele = [None if rng.random() < 0.2 else rng.randint(0, 1) for _ in range(n)]
    truth = [None if rng.random() < 0.25 else rng.randint(0, 1) for _ in range(n)]
    block = [rng.randint(0, 2) if a is not None else None for a in allele]
    phased = [k for k in range(n) if allele[k] is not None]
    link = [None] * n
    for a, b in zip(phased, phased[1:]):
        if block[a] == block[b]:
            link[a] = (rng.randint(0, 5), (rng.choice("xy"), pos[a] // 60), rng.choice("xy"))
    variants = [(pos[k], allele[k], truth[k], link[k]) for k in range(n)]
    genes = []
    for _ in range(rng.randint(0, 3)):
        s = rng.randint(0, 400)
        genes.append((s, s + rng.randint(0, 200), f"g{len(genes)}"))
    return variants, block, genes


def _brute_force(variants, block, genes, window):
    out = []
    n = len(variants)
    for i in range(n):
        for j in range(i + 1, n):
            (pi, ai, hi, _), (pj, aj, hj, _) = variants[i], variants[j]
            if hi is None or hj is None or pj - pi > window:
                continue
            if not any(s <= pi and pj <= e for s, e, _ in genes):
                continue
            if ai is None or aj is None:
                out.append((i, j, "unphased", None, 0, None, True))
                continue
            chain = [k for k in range(i, j + 1) if variants[k][1] is not None]
            if any(block[a] != block[b] for a, b in zip(chain, chain[1:])):
                out.append((i, j, "break", None, 0, None, True))
                continue
            links = [variants[k][3] for k in chain[:-1]]
            groups = {}
            for lv, _, g in links:
                groups[g] = max(groups.get(g, -1), lv)
            covered = all(variants[k][2] is not None for k in chain[1:-1])
            out.append((i, j, "chain", {c for _, c, _ in links},
                        max((lv for lv, _, _ in links), default=-1), groups, covered))
    return out


def test_pairs_match_brute_force():
    rng = random.Random(17)
    for _ in range(400):
        variants, block, genes = _random_chromosome(rng)
        window = rng.choice((30, 100, 1000))
        reach, holding = q4.gene_index(variants, genes)
        got = [(i, j, s, None if m is None else set(m), lv, None if g is None else dict(g), c)
               for i, j, s, m, lv, g, c in q4.pairs(variants, reach, window)]
        assert got == _brute_force(variants, block, genes, window)
        for k, (pos, *_) in enumerate(variants):
            assert sorted(holding[k]) == sorted((e, g) for s, e, g in genes if s <= pos <= e)


def test_fixed_chain():
    # Variants 0..4; 2 unphased; 3 has no truth; one block; levels 0, 2, 1 on links 0-1, 1-3, 3-4.
    variants = [(10, 0, 0, (0, ("easy", 0), "easy")),
                (20, 1, 1, (2, ("tandem", 0), "tandem")),
                (30, None, 0, None),
                (40, 0, None, (1, ("easy", 0), "easy")),
                (50, 1, 0, None)]
    got = list(q4.pairs(variants, [100] * 5, 100))
    assert [(i, j, s) for i, j, s, *_ in got] == [
        (0, 1, "chain"), (0, 2, "unphased"), (0, 4, "chain"), (1, 2, "unphased"),
        (1, 4, "chain"), (2, 4, "unphased")]
    _, _, _, cells, level, groups, covered = got[2]
    m = len(cells)
    assert (cells, level, covered) == ({("easy", 0), ("tandem", 0)}, 2, False)
    # Call 0-4: phaser trans (0 vs 1), truth (0 vs 0) cis: wrong.
    assert (0 ^ 0) != (1 ^ 0)
    # At r = 1e-2, m = 2 gives beta 3e-3 (index 3): level 2 is trusted, so certified.
    assert level <= q4.b_max(m, 1e-2)
    # At r = 1e-3, m = 2 gives beta 3e-4 (index 1): the tandem link (level 2) is not.
    assert level > q4.b_max(m, 1e-3)
