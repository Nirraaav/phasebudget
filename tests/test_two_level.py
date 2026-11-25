"""Two-level phasing against exhaustive orientation, and on the HG002 failure case."""

import itertools
import random

import pytest

from phasebudget import FragmentMatrix
from phasebudget.dp import CoverageDP
from phasebudget.selection import read_weights, select_optimal
from phasebudget.twolevel import orient_blocks, phased_cost

from conftest import random_matrix


def _one_level(matrix, cap):
    subset = matrix.subset(select_optimal(matrix, cap, read_weights(matrix)).selected)
    return subset, CoverageDP(subset).solve()


def _brute_force(matrix, subset, haplotype):
    blocks = subset.blocks()
    ids = sorted({b for b in blocks if b is not None}, key=blocks.index)
    best = float("inf")
    for flips in itertools.product((0, 1), repeat=max(0, len(ids) - 1)):
        orientation = dict(zip(ids, (0,) + flips))
        candidate = tuple(a ^ orientation[b] if b is not None else a
                          for a, b in zip(haplotype, blocks))
        best = min(best, phased_cost(matrix, candidate, blocks))
    return best


def test_the_chain_program_attains_the_exhaustive_minimum():
    rng = random.Random(7)
    checked = 0
    for _ in range(400):
        matrix = random_matrix(rng, max_variants=9, max_reads=9)
        for cap in (1, 2):
            subset, solution = _one_level(matrix, cap)
            if not len(subset):
                continue
            result = orient_blocks(matrix, subset, solution, max_window=64)
            assert result.dropped_reads == 0
            assert result.cost == pytest.approx(_brute_force(matrix, subset, solution.haplotype))
            checked += result.subset_blocks > 1
    assert checked > 50


def test_orientation_never_raises_the_full_matrix_cost():
    rng = random.Random(8)
    for _ in range(300):
        matrix = random_matrix(rng, max_variants=9, max_reads=9)
        subset, solution = _one_level(matrix, 1)
        if not len(subset):
            continue
        result = orient_blocks(matrix, subset, solution)
        assert result.cost <= phased_cost(matrix, solution.haplotype, subset.blocks()) + 1e-9


def test_the_hg002_junction_case_is_oriented_by_the_spanning_reads():
    # The selected reads leave a break between variants 2 and 3, so the subset
    # optimum's relation across it is arbitrary. Unselected reads spanning the
    # break all say the relation is "different"; two-level phasing must follow
    # them and join the two blocks.
    subset = FragmentMatrix.from_rows(["000---", "111---", "---111", "---000"])
    full = FragmentMatrix.from_rows(
        ["000---", "111---", "---111", "---000", "-0011-", "-1100-", "--01--"])
    solution = CoverageDP(subset).solve()
    result = orient_blocks(full, subset, solution)
    h = result.haplotype
    assert h[2] ^ h[3] == 1
    assert h[0] == h[1] == h[2] and h[3] == h[4] == h[5]
    assert len({b for b in result.phasing.blocks if b is not None}) == 1
    assert result.cost == 0.0


def test_reads_wider_than_the_window_are_dropped_and_counted():
    # Three blocks, with one read calling alleles in the first and the last.
    subset = FragmentMatrix.from_rows(["00-----", "---00--", "-----00"])
    full = FragmentMatrix.from_rows(["00-----", "---00--", "-----00", "0-----1"])
    solution = CoverageDP(subset).solve()
    assert orient_blocks(full, subset, solution, max_window=2).dropped_reads == 1
    assert orient_blocks(full, subset, solution, max_window=3).dropped_reads == 0


def test_variants_no_selected_read_covers_stay_unphased():
    subset = FragmentMatrix.from_rows(["00--", "11--"])
    full = FragmentMatrix.from_rows(["00--", "11--", "--01", "0-1-"])
    result = orient_blocks(full, subset, CoverageDP(subset).solve())
    assert result.phasing.alleles[2] is None and result.phasing.alleles[3] is None
