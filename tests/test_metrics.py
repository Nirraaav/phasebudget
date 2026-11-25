import pytest

from phasebudget import CoverageDP, FragmentMatrix, Phasing
from phasebudget.metrics import (
    an50,
    an50_lengths,
    block_lengths,
    block_n50,
    comparable_pairs,
    flip_errors,
    n50,
    phased_fraction,
    phasing_from_solution,
    qan50,
    qan50_lengths,
    summarize,
    switch_error_rate,
    switch_errors_excluding_flips,
    switch_positions,
)


def one_block(alleles):
    return Phasing(tuple(alleles), tuple(0 for _ in alleles))


def test_a_globally_inverted_block_has_no_switch_error():
    truth = (0, 1, 1, 0)
    assert switch_positions(truth, one_block((1, 0, 0, 1))) == ()
    assert switch_error_rate(truth, one_block((1, 0, 0, 1))) == 0.0


def test_one_switch_is_reported_at_the_variant_after_the_break():
    truth = (0, 0, 0, 0)
    assert switch_positions(truth, one_block((0, 0, 1, 1))) == (2,)


def test_two_adjacent_switches_are_one_flip():
    truth = (0, 0, 0, 0)
    phasing = one_block((0, 0, 1, 0))
    assert switch_positions(truth, phasing) == (2, 3)
    assert flip_errors(truth, phasing) == (2,)
    assert switch_errors_excluding_flips(truth, phasing) == 0


def test_a_long_switch_is_not_a_flip():
    truth = (0, 0, 0, 0, 0)
    phasing = one_block((0, 0, 1, 1, 1))
    assert flip_errors(truth, phasing) == ()
    assert switch_errors_excluding_flips(truth, phasing) == 1


def test_phase_is_compared_only_inside_a_block():
    truth = (0, 0, 0, 0)
    phasing = Phasing((0, 0, 1, 1), (0, 0, 1, 1))
    assert switch_positions(truth, phasing) == ()
    assert comparable_pairs(phasing) == 2


def test_phased_fraction_counts_called_variants():
    assert phased_fraction(Phasing((0, None, 1, 1), (0, None, 1, 1))) == 0.5


def test_block_lengths_in_variants_and_in_base_pairs():
    phasing = Phasing((0, 1, None, 0, 1), (0, 0, None, 1, 1))
    assert sorted(block_lengths(phasing)) == [2.0, 2.0]
    positions = (10, 20, 30, 100, 400)
    assert sorted(block_lengths(phasing, positions)) == [11.0, 301.0]


def test_n50_takes_the_block_that_covers_half_the_total():
    phasing = Phasing(
        (0, 1, 0, None, 1, 0, 1, 0, 1),
        (0, 0, 0, None, 1, 1, 1, 1, 1),
    )
    assert block_lengths(phasing) == (3.0, 5.0)
    assert block_n50(phasing) == 5.0


def test_an50_discounts_a_block_that_skips_variants():
    dense = Phasing((0, 1, 0, 1), (0, 0, 0, 0))
    sparse = Phasing((0, None, None, 1), (0, None, None, 0))
    assert an50(dense) == 4.0
    assert an50(sparse) == 2.0


def test_qan50_cuts_a_block_at_its_switch_error():
    truth = (0, 0, 0, 0)
    clean = one_block((0, 0, 0, 0))
    broken = one_block((0, 0, 1, 1))
    assert qan50(truth, clean) == 4.0
    assert qan50(truth, broken) == 2.0


def test_phasing_from_a_solution_marks_uncovered_variants_unphased():
    matrix = FragmentMatrix.from_rows(["01---", "01---", "---10"])
    phasing = phasing_from_solution(matrix, CoverageDP(matrix).solve())
    assert phasing.alleles[2] is None
    assert phasing.blocks == (0, 0, None, 1, 1)


def test_summarize_reports_every_metric():
    matrix = FragmentMatrix.from_rows(["0110", "0110", "1001"])
    solution = CoverageDP(matrix).solve()
    phasing = phasing_from_solution(matrix, solution)
    report = summarize((0, 1, 1, 0), phasing)
    assert report["switch_errors"] == 0.0
    assert report["phased_fraction"] == 1.0
    assert set(report) == {
        "phased_fraction",
        "comparable_pairs",
        "switch_errors",
        "switch_error_rate",
        "flip_errors",
        "switch_errors_excluding_flips",
        "blocks",
        "block_n50",
        "an50",
        "qan50",
    }


def test_truth_of_the_wrong_length_is_an_error():
    with pytest.raises(ValueError, match="same variants"):
        switch_positions((0, 1), one_block((0, 1, 0)))


def test_a_variant_the_truth_leaves_unphased_is_skipped_not_scored():
    # The prediction is wrong at variant 2 only, and the truth does not phase
    # it, so no switch may be charged there or at the pair after it.
    truth = (0, 1, None, 1)
    phasing = one_block((0, 1, 1, 1))
    assert switch_positions(truth, phasing) == ()
    assert comparable_pairs(phasing, truth) == 2
    assert comparable_pairs(phasing) == 3


def test_a_real_switch_across_a_skipped_variant_still_counts():
    truth = (0, 1, None, 0)
    phasing = one_block((0, 1, 0, 1))
    assert switch_positions(truth, phasing) == (3,)
    assert switch_error_rate(truth, phasing) == 0.5


def test_a_complete_truth_gives_the_same_answer_as_before():
    truth = (0, 1, 1, 0, 1)
    phasing = one_block((0, 1, 0, 0, 1))
    assert comparable_pairs(phasing, truth) == comparable_pairs(phasing)
    assert switch_positions(truth, phasing) == (2, 3)


def test_pooled_length_lists_reproduce_the_single_chromosome_n50s():
    truth = [0, 0, 1, 1, 0, None, 1]
    phasing = Phasing((0, 0, 1, 0, 0, 1, 1), (0, 0, 0, 0, 1, 1, 1))
    positions = [10, 20, 35, 50, 70, 90, 100]
    assert n50(an50_lengths(phasing, positions)) == an50(phasing, positions)
    assert n50(qan50_lengths(truth, phasing, positions)) == qan50(truth, phasing, positions)
