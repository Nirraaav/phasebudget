"""Window MEC, tail ranks, tau and k-mer spans of `benchmarks/q2_hard_flags.py` (Q2)."""

import math
import random
import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")  # not a core dependency; the whatshap job has it

from phasebudget._types import Fragment
from phasebudget.exact import solve_bruteforce
from phasebudget.fragments import FragmentMatrix

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "benchmarks"))

import q2_hard_flags as m  # noqa: E402


def restricted(matrix: FragmentMatrix, cols: list[int]) -> FragmentMatrix:
    """The oracle's input: every read's calls at `cols`, columns renumbered 0..len(cols)-1."""
    fragments = []
    for f in matrix:
        calls = [(i, f.allele_at(c), f.weight_at(c)) for i, c in enumerate(cols)
                 if f.allele_at(c) is not None]
        if not calls:
            continue
        start, end = calls[0][0], calls[-1][0] + 1
        alleles = [None] * (end - start)
        weights = [0.0] * (end - start)
        for i, a, w in calls:
            alleles[i - start] = a
            weights[i - start] = w
        fragments.append(Fragment(len(fragments), start, tuple(alleles), tuple(weights)))
    return FragmentMatrix(fragments, len(cols))


def random_matrix(rng: random.Random, n: int, reads: int) -> FragmentMatrix:
    fragments = []
    for r in range(reads):
        start = rng.randrange(n)
        end = rng.randrange(start + 1, min(n, start + 6) + 1)
        alleles = [rng.choice((0, 1, None)) if 0 < k < end - start - 1 else rng.choice((0, 1))
                   for k in range(end - start)]
        weights = [0.0 if a is None else float(rng.choice((1, 2, 3, 30))) for a in alleles]
        fragments.append(Fragment(r, start, tuple(alleles), tuple(weights)))
    return FragmentMatrix(fragments, n)


def test_window_mec_fixed():
    # Three reads over two columns: two say "same", one says "opposite" with weight 2.
    alleles = np.array([[0, 0], [1, 1], [0, 1]], np.int8)
    weights = np.array([[1.0, 1.0], [1.0, 1.0], [2.0, 2.0]])
    cost, total = m.window_mec(alleles, weights)
    assert cost == 2.0  # flip one call of the third read, or both calls' worth of the others
    assert total == 8.0
    # A single-call read costs nothing and does not count in the total.
    cost, total = m.window_mec(np.array([[0, -1], [1, 1]], np.int8), np.ones((2, 2)))
    assert (cost, total) == (0.0, 2.0)
    assert m.window_mec(np.array([[0, -1]], np.int8), np.ones((1, 2))) == (0.0, 0.0)


@pytest.mark.parametrize("seed", range(40))
def test_window_mec_matches_exact_on_random_windows(seed):
    rng = random.Random(seed)
    n = rng.randrange(3, 12)
    matrix = random_matrix(rng, n, rng.randrange(1, 14))
    calls = m.flat_calls(matrix)
    a = rng.randrange(n - 1)
    b = rng.randrange(a + 1, n)
    cols = m.window_columns(a, b, n)
    cost, _ = m.window_mec(*m.window_arrays(calls, cols))
    assert math.isclose(cost, solve_bruteforce(restricted(matrix, cols)).cost, abs_tol=1e-9)


def test_flat_calls_balance_counts():
    matrix = FragmentMatrix([Fragment(0, 0, (0, None, 1), (1.0, 0.0, 1.0)),
                             Fragment(1, 1, (1, 1), (1.0, 1.0))], 3)
    offsets, col, read, allele, _ = m.flat_calls(matrix)
    assert list(offsets) == [0, 1, 2, 4]
    assert list(col) == [0, 1, 2, 2]
    assert sorted(zip(read[2:], allele[2:])) == [(0, 1), (1, 1)]


def test_window_columns():
    assert m.window_columns(5, 6, 20) == [3, 4, 5, 6, 7, 8]
    assert m.window_columns(0, 1, 20) == [0, 1, 2, 3]
    assert m.window_columns(18, 19, 20) == [16, 17, 18, 19]
    assert m.window_columns(5, 9, 20) == [3, 4, 5, 9, 10, 11]  # columns inside the pair skipped


def test_tail_rank_and_tau_with_ties():
    ref = np.sort(np.array([0.0, 0.0, 1.0, 2.0, 2.0, 3.0, -math.inf, 5.0, 5.0, 9.0]))
    ranks = m.tail_rank(ref, np.array([9.0, 5.0, 0.0, -math.inf, 10.0]))
    assert list(ranks) == [0.1, 0.3, 0.9, 1.0, 0.0]
    r = np.array([0.1, 0.3, 0.3, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    assert m.choose_tau(r, 0.2) == 0.1  # 0.3 would flag four of ten
    assert m.choose_tau(r, 0.4) == 0.3
    assert m.choose_tau(np.array([0.3, 0.3, 0.3]), 0.09) == -math.inf


def test_flag_uses_the_smallest_rank():
    ref = np.array([[float(i), float(9 - i), 0.0, 0.0] for i in range(10)])
    f = np.array([[9.0, 0.0, 0.0, 0.0], [0.0, 9.0, 0.0, 0.0], [5.0, 5.0, 0.0, 0.0]])
    assert list(m.flag(f, ref, ("depth",), 0.1)) == [True, False, False]
    assert list(m.flag(f, ref, ("depth", "balance"), 0.1)) == [True, True, False]


def test_balance_score():
    pytest.importorskip("scipy")
    s = m.balance_score(np.array([5, 0, 0, 3]), np.array([5, 10, 0, 1]))
    assert s[0] == 0.0
    assert math.isclose(s[1], -math.log10(2 * 0.5 ** 10))
    assert s[2] == 0.0
    assert math.isclose(s[3], -math.log10(2 * 5 / 16))


def test_kmer_spans_and_canonical():
    assert m.kmer_spans(100, "A", k=5) == ((94, 99), (100, 105))
    assert m.kmer_spans(100, "ACG", k=5) == ((94, 99), (102, 107))  # deletion: after the ref
    assert m.canonical("TTTA") == "TAAA"
    assert m.canonical("AAAT") == "AAAT"
