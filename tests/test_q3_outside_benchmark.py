"""Windows, constrained optima, causes and counts of `benchmarks/q3_outside_benchmark.py` (Q3)."""

import importlib.util
import random
import sys
from itertools import product
from pathlib import Path

from conftest import random_matrix

from phasebudget._types import Fragment
from phasebudget.fragments import FragmentMatrix

BENCH = Path(__file__).resolve().parent.parent / "benchmarks"


def _module():
    sys.path.insert(0, str(BENCH))
    spec = importlib.util.spec_from_file_location("q3_outside_benchmark",
                                                  BENCH / "q3_outside_benchmark.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def oracle(matrix: FragmentMatrix, a: int, b: int, flank: int, relation: int) -> float:
    """Best complementary-model cost over every haplotype of all columns, counting only the
    kept columns, with `h[a] ^ h[b] == relation`. Independent of `window_matrix`."""
    kept = set(range(max(0, a - flank), a + 1)) | set(range(b, min(matrix.n_variants,
                                                                    b + flank + 1)))
    best = float("inf")
    for h in product((0, 1), repeat=matrix.n_variants):
        if h[a] ^ h[b] != relation:
            continue
        total = 0.0
        for f in matrix:
            sides = [0.0, 0.0]
            for j, allele, weight in f.observed():
                if j in kept:
                    sides[0] += weight * (allele != h[j])
                    sides[1] += weight * (allele != 1 - h[j])
            total += min(sides)
        best = min(best, total)
    return best


def test_window_drops_the_columns_between_and_renumbers():
    m = _module()
    matrix = FragmentMatrix([
        Fragment(7, 0, (0, 1, None, 1, 0), (1.0, 2.0, 0.0, 3.0, 4.0)),
        Fragment(8, 2, (1, 0), (5.0, 6.0)),  # calls 2 and 3; 2 is dropped
        Fragment(9, 2, (1,), (1.0,)),  # calls only a dropped column: left out
        Fragment(10, 6, (1, 1), (1.0, 1.0)),  # beyond the flank: left out
    ], 8)
    window, boundary, columns = m.window_matrix(matrix, 1, 3, flank=1)
    assert columns == [0, 1, 3, 4]
    assert boundary == 1
    assert [f.read_id for f in window] == [7, 8]
    assert window[0].start == 0 and window[0].alleles == (0, 1, 1, 0)
    assert window[0].weights == (1.0, 2.0, 3.0, 4.0)
    assert window[1].start == 2 and window[1].alleles == (0,) and window[1].weights == (6.0,)
    # a gap inside the kept columns stays a gap
    gapped = FragmentMatrix([Fragment(1, 0, (0, None, 1, 1), (1.0, 0.0, 1.0, 1.0))], 4)
    w, k, cols = m.window_matrix(gapped, 2, 3, flank=2)
    assert cols == [0, 1, 2, 3] and k == 2
    assert w[0].alleles == (0, None, 1, 1)


def test_candidates_give_the_same_window():
    m = _module()
    rng = random.Random(5)
    for _ in range(200):
        matrix = random_matrix(rng, max_variants=10, max_reads=8)
        if matrix.n_variants < 2:
            continue
        a = rng.randrange(matrix.n_variants - 1)
        b = rng.randrange(a + 1, matrix.n_variants)
        near = [i for i, f in enumerate(matrix) if f.end > a - 2 and f.start <= b + 2]
        full = m.window_matrix(matrix, a, b, 2)
        some = m.window_matrix(matrix, a, b, 2, near)
        assert full[1:] == some[1:]
        assert full[0].fragments == some[0].fragments


def test_relation_costs_match_the_oracle():
    m = _module()
    rng = random.Random(11)
    checked = 0
    for _ in range(300):
        matrix = random_matrix(rng, max_variants=8, max_reads=7)
        if matrix.n_variants < 2:
            continue
        a = rng.randrange(matrix.n_variants - 1)
        b = rng.randrange(a + 1, matrix.n_variants)
        flank = rng.randint(0, 3)
        window, boundary, columns = m.window_matrix(matrix, a, b, flank)
        costs = m.relation_costs(window, boundary, cap=None)
        for relation in (0, 1):
            assert abs(costs[relation] - oracle(matrix, a, b, flank, relation)) < 1e-9
        # a cap at or above the coverage keeps every phase-informative read
        capped = m.relation_costs(window, boundary, cap=max(1, window.max_coverage()),
                                  positions=[100 * c for c in columns])
        assert all(abs(x - y) < 1e-9 for x, y in zip(costs, capped))
        checked += 1
    assert checked > 200


def test_empty_window_has_equal_costs():
    m = _module()
    matrix = FragmentMatrix([Fragment(0, 5, (1,), (1.0,))], 6)
    window, boundary, _ = m.window_matrix(matrix, 0, 1, flank=1)
    assert len(window) == 0
    assert m.relation_costs(window, boundary) == (0.0, 0.0)


def test_votes_optima_and_causes():
    m = _module()
    assert m.majority(0, 0) == "none"
    assert m.majority(5, 2) == "wrong"  # agree = reads backing the phaser's wrong call
    assert m.majority(3, 3) == "tie"
    assert m.majority(1, 4) == "right"
    assert m.optimum((2.0, 5.0), 0) == "right"
    assert m.optimum((2.0, 5.0), 1) == "wrong"
    assert m.optimum((3.0, 3.0), 1) == "tie"
    assert m.cause("none", "tie") == "none"
    assert m.cause("wrong", "right") == "search"
    assert m.cause("right", "right") == "search"
    assert m.cause("wrong", "wrong") == "data"
    assert m.cause("tie", "tie") == "data"
    assert m.cause("right", "wrong") == "other"
    assert m.cause("right", "tie") == "other"


def test_unit_shares_and_gate():
    m = _module()
    rows = [("chr1", 5, False), ("chr1", 99_999, True), ("chr1", 100_000, False),
            ("chr2", 5, False)]
    assert m.unit_shares(rows) == (4, 1, 3, 1)
    inside = (10_000, 10, 100, 5)  # rate 1e-3
    ratio, share, ok = m.gate(inside, (1_000, 4, 50, 3))
    assert ok and abs(ratio - 4.0) < 1e-9 and abs(share - 4 / 14) < 1e-12
    ratio, share, ok = m.gate(inside, (1_000, 2, 50, 2))
    assert not ok and abs(ratio - 2.0) < 1e-12
    ratio, share, ok = m.gate(inside, (100, 2, 5, 2))  # ratio 20, share 2/12 < 0.2
    assert not ok and share < m.GATE_SHARE


def test_truth_check():
    m = _module()
    q100 = {("c", 10, "A", "G"): "0|1", ("c", 20, "C", "T"): "1|0", ("c", 30, "G", "GA"): "0|1",
            ("c", 40, "T", "C"): "0|1"}
    v5 = {("c", 10, "A", "G"): "1|0", ("c", 20, "C", "T"): "0|1",  # same relation (trans)
          ("c", 30, "G", "GA"): "1|1",  # genotype differs
          ("c", 40, "T", "C"): "0/1"}
    counts = m.truth_check([("c", 10, 20), ("c", 20, 30), ("c", 10, 40)], q100, v5)
    assert counts == {"genotype differs": 1, "same genotype, same relation": 1,
                      "same genotype, other relation": 0, "same genotype, v5.0q unphased": 1}
    v5[("c", 20, "C", "T")] = "1|0"
    counts = m.truth_check([("c", 10, 20)], q100, v5)
    assert counts["same genotype, other relation"] == 1
    del v5[("c", 10, "A", "G")]
    assert m.truth_check([("c", 10, 20)], q100, v5)["genotype differs"] == 1
