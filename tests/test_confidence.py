"""Read support for phased relations, the Wilson score, and Learn then Test."""

import math
import random

import pytest

from phasebudget import Phasing
from phasebudget._types import Fragment
from phasebudget.confidence import (
    ChainCounts,
    PairSupport,
    Units,
    binomial_cdf,
    chain_counts,
    coverage_grid,
    learn_then_test,
    learn_then_test_by_group,
    pair_support,
    phased_pairs,
    score,
    wilson_lower,
)
from phasebudget.fragments import FragmentMatrix


def _matrix(reads, n):
    """Reads as `{variant: allele}` dicts, turned into fragments with unit weights."""
    fragments = []
    for i, calls in enumerate(reads):
        start, end = min(calls), max(calls) + 1
        alleles = tuple(calls.get(j) for j in range(start, end))
        weights = tuple(1.0 if a is not None else 0.0 for a in alleles)
        fragments.append(Fragment(i, start, alleles, weights))
    return FragmentMatrix(fragments, n)


def _oracle(reads, phasing):
    """Pair support counted pair by pair, read by read."""
    out = []
    for a, b in phased_pairs(phasing):
        agree = disagree = 0
        for calls in reads:
            if a in calls and b in calls:
                if (calls[a] == calls[b]) == (phasing.alleles[a] == phasing.alleles[b]):
                    agree += 1
                else:
                    disagree += 1
        out.append(PairSupport(a, b, agree, disagree))
    return tuple(out)


def test_fixed_example():
    phasing = Phasing(alleles=(0, 1, None, 0, 1), blocks=(0, 0, None, 0, 1))
    reads = [{0: 0, 1: 1}, {0: 1, 1: 0}, {0: 0, 1: 0}, {1: 1, 3: 0}, {3: 0, 4: 1}]
    matrix = _matrix(reads, 5)
    support = pair_support(matrix, phasing)
    # Pairs (0, 1) and (1, 3); variant 2 is unphased, and 3 and 4 sit in different blocks.
    assert [(s.first, s.second) for s in support] == [(0, 1), (1, 3)]
    assert (support[0].agree, support[0].disagree) == (2, 1)
    assert (support[1].agree, support[1].disagree) == (1, 0)


def test_pair_support_matches_oracle_on_random_instances():
    rng = random.Random(7)
    for _ in range(300):
        n = rng.randint(2, 12)
        reads = []
        for _ in range(rng.randint(0, 15)):
            start = rng.randrange(n)
            end = rng.randint(start, n - 1)
            calls = {j: rng.randint(0, 1) for j in range(start, end + 1)
                     if rng.random() < 0.7 or j in (start, end)}
            reads.append(calls)
        alleles = tuple(None if rng.random() < 0.2 else rng.randint(0, 1) for _ in range(n))
        blocks = tuple(None if a is None else rng.randint(0, 2) for a in alleles)
        phasing = Phasing(alleles=alleles, blocks=blocks)
        assert pair_support(_matrix(reads, n), phasing) == _oracle(reads, phasing)


def test_pair_support_checks_lengths():
    with pytest.raises(ValueError):
        pair_support(_matrix([{0: 0, 1: 1}], 2), Phasing(alleles=(0, 1, 0), blocks=(0, 0, 0)))


def test_wilson_known_values_and_order():
    assert wilson_lower(0, 0) == -1.0
    assert wilson_lower(1, 1) == pytest.approx(0.2)
    # p = 1/2, n = 4, z = 2: centre 1, spread 2 * sqrt(1/16 + 1/16), denominator 2.
    assert wilson_lower(2, 4) == pytest.approx((1 - 2 * math.sqrt(1 / 8)) / 2)
    assert wilson_lower(1, 1) < wilson_lower(30, 30) < 1
    assert wilson_lower(24, 38) < wilson_lower(30, 30)
    with pytest.raises(ValueError):
        wilson_lower(3, 2)


def test_score_breaks_ties_by_reads():
    assert score(PairSupport(0, 1, 5, 0)) < score(PairSupport(0, 1, 6, 0))
    assert score(PairSupport(0, 1, 3, 3))[0] < score(PairSupport(0, 1, 3, 0))[0]


def test_binomial_cdf_matches_direct_sum():
    for n in range(1, 40):
        for p in (0.01, 0.2, 0.5):
            for k in range(-1, n + 1):
                direct = sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i)
                             for i in range(max(0, k + 1)))
                assert binomial_cdf(k, n, p) == pytest.approx(min(1.0, direct), abs=1e-12)


def _naive_units_at(rows, threshold):
    trusted, wrong = set(), set()
    for unit, key, bad in rows:
        if key >= threshold:
            trusted.add(unit)
            if bad:
                wrong.add(unit)
    return len(trusted), len(wrong)


def test_units_and_learn_then_test_match_naive_scans():
    rng = random.Random(3)
    for _ in range(100):
        rows = [(rng.randrange(40), (rng.random(), rng.randint(0, 9)), rng.random() < 0.05)
                for _ in range(rng.randint(1, 400))]
        units = Units.from_pairs(rows)
        grid = coverage_grid([key for _, key, _ in rows], start=0.5, step=0.05)
        for threshold in grid:
            assert units.at(threshold) == _naive_units_at(rows, threshold)
        for beta in (0.01, 0.1, 0.3):
            chosen = None
            for threshold in grid:
                n, k = _naive_units_at(rows, threshold)
                if n == 0 or binomial_cdf(k, n, beta) > 0.1:
                    break
                chosen = threshold
            assert learn_then_test(units, grid, beta, 0.1) == chosen


def test_coverage_grid_shares():
    scores = [(float(i), 0) for i in range(100)]
    grid = coverage_grid(scores, start=0.5, step=0.25)
    assert grid == ((50.0, 0), (25.0, 0), (0.0, 0))
    assert coverage_grid([]) == ()


def test_learn_then_test_certifies_clean_units_and_refuses_dirty_ones():
    clean = Units.from_pairs([(u, (1.0, 1), False) for u in range(10_000)])
    assert learn_then_test(clean, [(1.0, 1)], beta=1e-3, delta=0.05) == (1.0, 1)
    dirty = Units.from_pairs([(u, (1.0, 1), u % 10 == 0) for u in range(10_000)])
    assert learn_then_test(dirty, [(1.0, 1)], beta=1e-3, delta=0.05) is None


def test_learn_then_test_by_group_matches_naive_scan_per_group():
    rng = random.Random(5)
    for _ in range(60):
        rows = [(rng.choice("abc"), rng.randrange(30), (rng.random(), rng.randint(0, 9)),
                 rng.random() < rng.choice((0.0, 0.02, 0.2)))
                for _ in range(rng.randint(1, 300))]
        for beta in (0.05, 0.3):
            chosen = learn_then_test_by_group(rows, beta, 0.1, start=0.5, step=0.05)
            assert set(chosen) == {g for g, _, _, _ in rows}
            for group in chosen:
                own = [(u, key, bad) for g, u, key, bad in rows if g == group]
                expected = None
                for threshold in coverage_grid([key for _, key, _ in own], 0.5, 0.05):
                    n, k = _naive_units_at(own, threshold)
                    if n == 0 or binomial_cdf(k, n, beta) > 0.1:
                        break
                    expected = threshold
                assert chosen[group] == expected


def test_learn_then_test_by_group_keeps_a_dirty_group_from_moving_a_clean_one():
    clean = [("easy", u, (1.0, 1), False) for u in range(10_000)]
    dirty = [("segdup", u, (1.0, 1), u % 10 == 0) for u in range(10_000)]
    chosen = learn_then_test_by_group(clean + dirty, beta=1e-3, delta=0.05)
    assert chosen == {"easy": (1.0, 1), "segdup": None}
    pooled = learn_then_test(Units.from_pairs([(u, key, bad) for _, u, key, bad in clean + dirty]),
                             [(1.0, 1)], beta=1e-3, delta=0.05)
    assert pooled is None


def _naive_chain_counts(links, phaser, truth, bins):
    """Every pair of truth-phased variants, judged by the relation's definition."""
    out = [0, 0, 0, 0, 0]
    ends = [k for k, t in enumerate(truth) if t is not None]
    for a, i in enumerate(ends):
        for j in ends[a + 1:]:
            wrong = (phaser[i] == phaser[j]) != (truth[i] == truth[j])
            trusted = all(links[i:j])
            out[0] += 1
            out[1] += wrong
            if trusted:
                out[2] += 1
                out[3] += wrong
                out[4] += bins[j] - bins[i] + 1
    return ChainCounts(*out)


def test_chain_counts_matches_enumeration_of_pairs():
    rng = random.Random(11)
    for _ in range(500):
        n = rng.randint(0, 25)
        phaser = [rng.randint(0, 1) for _ in range(n)]
        truth = [None if rng.random() < 0.2 else rng.randint(0, 1) for _ in range(n)]
        links = [rng.random() < 0.8 for _ in range(max(0, n - 1))]
        bins = sorted(rng.randint(0, 6) for _ in range(n))
        labels = [None if t is None else p ^ t for p, t in zip(phaser, truth)]
        assert chain_counts(links, labels, bins) == _naive_chain_counts(links, phaser, truth, bins)


def test_chain_counts_fixed_example_and_checks():
    # Variants 0..3; the link 1-2 is untrusted; variant 2 has no truth.
    counts = chain_counts([True, False, True], [0, 1, None, 1], [0, 0, 1, 1])
    assert counts == ChainCounts(pairs=3, wrong=2, trusted=1, wrong_trusted=1, trusted_span=1)
    with pytest.raises(ValueError):
        chain_counts([True], [0, 1, 0], [0, 0, 0])
