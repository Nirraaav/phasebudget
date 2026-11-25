import random

import pytest

from phasebudget import FragmentMatrix
from phasebudget._types import Margin
from phasebudget._validation import CoverageLimitExceeded
from phasebudget.adaptive import (
    MARGIN_SOURCES,
    adaptive_selection,
    fragile_boundaries,
    static_selection,
)
from phasebudget.dp import CoverageDP
from phasebudget.selection import is_feasible, read_weights, select_greedy, select_optimal
from phasebudget.simulate import simulate

from conftest import random_matrix


def margins(values, spanned=None):
    if spanned is None:
        spanned = [True] * len(values)
    return tuple(
        Margin(boundary=index, relation=0, margin=value, spanned=flag)
        for index, (value, flag) in enumerate(zip(values, spanned))
    )


def test_fragile_boundaries_take_the_weakest_decisions():
    assert fragile_boundaries(margins([5.0, 1.0, 4.0, 0.5]), quantile=0.5) == (1, 3)


def test_fragile_boundaries_ignore_unspanned_gaps():
    chosen = fragile_boundaries(
        margins([0.0, 3.0, 2.0], spanned=[False, True, True]), quantile=0.5
    )
    assert 0 not in chosen


def test_fragile_boundaries_always_report_at_least_one():
    assert len(fragile_boundaries(margins([2.0] * 10), quantile=0.01)) == 1


def test_fragile_boundaries_of_a_matrix_without_spans_is_empty():
    assert fragile_boundaries(()) == ()


def test_the_loop_never_returns_a_worse_solution_than_its_first_round():
    generator = random.Random(31)
    for _ in range(40):
        matrix = random_matrix(generator, max_variants=8, max_reads=8)
        coverage = 2
        first = static_selection(matrix, coverage)
        loop = adaptive_selection(matrix, coverage, iterations=4)
        assert loop.full_cost <= first.full_cost + 1e-9


def test_every_round_respects_the_coverage_cap():
    instance = simulate(16, coverage=9.0, read_span=6, error_rate=0.06, seed=33)
    result = adaptive_selection(instance.matrix, 3, iterations=5)
    for step in result.history:
        assert is_feasible(instance.matrix, step.selection.selected, 3)


def test_the_loop_stops_when_the_selection_repeats():
    instance = simulate(14, coverage=8.0, read_span=5, error_rate=0.05, seed=34)
    result = adaptive_selection(instance.matrix, 3, iterations=10)
    assert result.iterations <= 10
    selections = [step.selection.selected for step in result.history]
    assert len(selections) == len(set(selections)) or selections[-1] in selections[:-1]


def test_the_reported_round_is_the_best_round():
    instance = simulate(18, coverage=10.0, read_span=6, error_rate=0.08, seed=35)
    result = adaptive_selection(instance.matrix, 3, iterations=5)
    assert result.full_cost == min(step.full_cost for step in result.history)
    assert result.selection in [step.selection for step in result.history]


def test_boosted_reads_span_a_fragile_boundary():
    instance = simulate(16, coverage=9.0, read_span=6, error_rate=0.05, seed=36)
    result = adaptive_selection(instance.matrix, 3, iterations=3)
    step = result.history[0]
    index = {fragment.read_id: fragment for fragment in instance.matrix}
    for read_id in step.boosted:
        assert any(index[read_id].spans(boundary) for boundary in step.fragile)


def test_a_single_iteration_is_static_selection():
    instance = simulate(14, coverage=8.0, read_span=5, error_rate=0.05, seed=37)
    loop = adaptive_selection(instance.matrix, 3, iterations=1)
    static = static_selection(instance.matrix, 3)
    assert loop.selection.selected == static.selection.selected
    assert loop.iterations == 1


def test_a_different_selector_can_be_plugged_in():
    instance = simulate(14, coverage=8.0, read_span=5, error_rate=0.05, seed=38)
    result = adaptive_selection(
        instance.matrix, 3, iterations=3, selector=select_greedy
    )
    assert all(step.selection.method == "greedy" for step in result.history)


def test_a_boost_of_one_cannot_move_the_loop():
    matrix = FragmentMatrix.from_rows(["0110", "1001"])
    with pytest.raises(ValueError, match="boost"):
        adaptive_selection(matrix, 2, boost=1.0)


def test_zero_iterations_is_an_error():
    matrix = FragmentMatrix.from_rows(["0110"])
    with pytest.raises(ValueError, match="at least 1"):
        adaptive_selection(matrix, 1, iterations=0)


def divergent_instance():
    """Three reads where selection at cap 1 drops the only read spanning gap 1.

    Reads 0 and 1 are heavy and disjoint, so exact selection takes both and
    leaves read 2 out. Read 2 is the only one covering variants 1 and 2
    together, so boundary 1 is a real decision in the full matrix and a block
    boundary in the subset. This is the OQ-10 confound at minimum size.
    """
    return FragmentMatrix.from_rows(
        ["01--", "--01", "-10-"],
        [[3.0] * 4, [3.0] * 4, [0.5] * 4],
    )


def test_the_subset_can_hide_a_boundary_the_full_matrix_spans():
    matrix = divergent_instance()
    selection = select_optimal(matrix, 1, read_weights(matrix))
    assert selection.selected == (0, 1)
    subset = matrix.subset(selection.selected)
    assert matrix.spanned_boundaries()[1] is True
    assert subset.spanned_boundaries()[1] is False
    full = CoverageDP(matrix).margins()[1]
    partial = CoverageDP(subset).margins()[1]
    assert full.spanned and full.margin == pytest.approx(0.5)
    assert not partial.spanned and partial.margin == pytest.approx(0.0)


def test_the_two_margin_sources_boost_different_reads():
    matrix = divergent_instance()
    subset_loop = adaptive_selection(matrix, 1, iterations=2, margin_source="subset")
    full_loop = adaptive_selection(matrix, 1, iterations=2, margin_source="full")
    assert subset_loop.history[0].fragile == (0,)
    assert subset_loop.history[0].boosted == (0,)
    # Read 2 is the one that could firm up boundary 1, and only the full
    # source ever asks for it.
    assert full_loop.history[0].fragile == (1,)
    assert full_loop.history[0].boosted == (2,)


def test_full_margins_are_the_margins_of_the_whole_matrix():
    instance = simulate(12, coverage=6.0, read_span=5, error_rate=0.06, seed=51)
    expected = CoverageDP(instance.matrix).margins()
    result = adaptive_selection(
        instance.matrix, 3, iterations=3, margin_source="full"
    )
    assert result.history[0].margins == expected


def test_full_margins_do_not_change_between_rounds():
    instance = simulate(16, coverage=9.0, read_span=6, error_rate=0.08, seed=52)
    result = adaptive_selection(
        instance.matrix, 3, iterations=4, margin_source="full"
    )
    first = result.history[0]
    for step in result.history[1:]:
        assert step.margins == first.margins
        assert step.fragile == first.fragile


def test_the_first_round_ignores_the_margin_source():
    instance = simulate(14, coverage=8.0, read_span=5, error_rate=0.07, seed=53)
    subset_loop = adaptive_selection(instance.matrix, 3, margin_source="subset")
    full_loop = adaptive_selection(instance.matrix, 3, margin_source="full")
    assert (
        subset_loop.history[0].selection.selected
        == full_loop.history[0].selection.selected
    )


def test_the_full_margin_loop_respects_the_coverage_cap():
    instance = simulate(16, coverage=9.0, read_span=6, error_rate=0.06, seed=54)
    result = adaptive_selection(
        instance.matrix, 3, iterations=5, margin_source="full"
    )
    for step in result.history:
        assert is_feasible(instance.matrix, step.selection.selected, 3)


def test_the_full_margin_loop_never_returns_a_worse_first_round():
    generator = random.Random(41)
    for _ in range(40):
        matrix = random_matrix(generator, max_variants=8, max_reads=8)
        first = static_selection(matrix, 2)
        loop = adaptive_selection(matrix, 2, iterations=4, margin_source="full")
        assert loop.full_cost <= first.full_cost + 1e-9


def test_the_margin_source_is_recorded_on_the_result_and_every_step():
    instance = simulate(12, coverage=7.0, read_span=5, error_rate=0.05, seed=55)
    for source in MARGIN_SOURCES:
        result = adaptive_selection(
            instance.matrix, 3, iterations=3, margin_source=source
        )
        assert result.margin_source == source
        assert all(step.margin_source == source for step in result.history)


def test_an_unknown_margin_source_is_an_error():
    matrix = FragmentMatrix.from_rows(["0110", "1001"])
    with pytest.raises(ValueError, match="margin_source"):
        adaptive_selection(matrix, 2, margin_source="oracle")


def test_a_full_margin_source_refuses_a_matrix_over_the_state_limit():
    instance = simulate(10, coverage=12.0, read_span=6, error_rate=0.05, seed=56)
    with pytest.raises(CoverageLimitExceeded):
        adaptive_selection(
            instance.matrix, 2, margin_source="full", max_active_reads=4
        )


def test_the_signal_is_validated():
    matrix = simulate(12, coverage=6, read_span=4, error_rate=0.1, seed=2).matrix
    with pytest.raises(ValueError):
        adaptive_selection(matrix, 3, signal="entropy")


def test_the_posterior_signal_never_reports_worse_than_static():
    for seed in range(15):
        matrix = simulate(20, coverage=8, read_span=4, error_rate=0.1, seed=300 + seed).matrix
        static = static_selection(matrix, 3)
        steered = adaptive_selection(matrix, 3, iterations=3, signal="posterior")
        assert steered.full_cost <= static.full_cost + 1e-9
        assert all(len(step.margins) == matrix.n_variants - 1 for step in steered.history)


def test_an_absolute_threshold_takes_every_boundary_below_it_and_no_other():
    margins = tuple(
        Margin(boundary=b, relation=0, margin=value, spanned=spanned)
        for b, (value, spanned) in enumerate([(0.5, True), (3.0, True), (0.1, False), (2.9, True)])
    )
    assert fragile_boundaries(margins, threshold=2.94) == (0, 3)
    assert fragile_boundaries(margins, threshold=0.01) == ()


def test_the_probability_rule_needs_the_posterior_signal():
    matrix = simulate(12, coverage=6, read_span=4, error_rate=0.1, seed=2).matrix
    with pytest.raises(ValueError):
        adaptive_selection(matrix, 3, error_probability=0.05)
    with pytest.raises(ValueError):
        adaptive_selection(matrix, 3, signal="posterior", error_probability=0.7)
    result = adaptive_selection(matrix, 3, signal="posterior", error_probability=0.05)
    assert result.full_cost <= static_selection(matrix, 3).full_cost + 1e-9


def test_stitching_one_round_returns_it():
    from phasebudget.adaptive import stitch_haplotypes

    matrix = simulate(30, coverage=6, read_span=4, error_rate=0.1, seed=4).matrix
    h = CoverageDP(matrix).solve().haplotype
    stitched, choice = stitch_haplotypes(matrix, [h], window=7)
    assert stitched == h
    assert set(choice) == {0}


def test_stitching_picks_the_better_piece_in_each_window():
    from phasebudget.adaptive import stitch_haplotypes

    # Error-free reads confined to each half, plus reads across the middle.
    # Each candidate is right on one half and scrambled on the other, so the
    # stitch must take window 0 from the first and window 1 from the second.
    truth = (0, 1, 1, 0, 1, 0, 0, 1)
    text = "".join(str(a) for a in truth)
    comp = "".join(str(1 - a) for a in truth)
    rows = [text[:4] + "----", comp[:4] + "----", "----" + text[4:], "----" + comp[4:],
            "--" + text[2:6] + "--", "--" + comp[2:6] + "--"] * 2
    matrix = FragmentMatrix.from_rows(rows)
    scramble = lambda part: tuple(1 - a if i % 2 else a for i, a in enumerate(part))
    left_right = truth[:4] + scramble(truth[4:])
    right_right = scramble(truth[:4]) + truth[4:]
    stitched, choice = stitch_haplotypes(matrix, [left_right, right_right], window=4)
    assert choice == (0, 1)
    assert matrix.cost(stitched) == 0.0


def test_stitched_acceptance_is_never_worse_than_global():
    for seed in range(12):
        matrix = simulate(60, coverage=10, read_span=5, error_rate=0.15, seed=700 + seed).matrix
        kwargs = dict(iterations=4, signal="posterior", error_probability=0.05)
        global_rule = adaptive_selection(matrix, 3, **kwargs)
        stitched = adaptive_selection(matrix, 3, acceptance="stitched", window=10, **kwargs)
        assert stitched.full_cost <= global_rule.full_cost + 1e-9
        assert stitched.full_cost == pytest.approx(matrix.cost(stitched.solution.haplotype))
        assert set(global_rule.support) == set(global_rule.selected)
        assert stitched.acceptance == "stitched"


def test_the_acceptance_rule_is_validated():
    matrix = simulate(12, coverage=6, read_span=4, error_rate=0.1, seed=2).matrix
    with pytest.raises(ValueError):
        adaptive_selection(matrix, 3, acceptance="majority")


def test_the_shuffled_control_keeps_the_count_and_moves_the_places():
    from phasebudget.posterior import posterior_margins

    matrix = simulate(60, coverage=10, read_span=5, error_rate=0.15, seed=21).matrix
    kwargs = dict(iterations=2, error_probability=0.05)
    real = adaptive_selection(matrix, 3, signal="posterior", **kwargs)
    control = adaptive_selection(matrix, 3, signal="shuffled", shuffle_seed=5, **kwargs)
    first_real, first_control = real.history[0], control.history[0]
    assert len(first_real.fragile) == len(first_control.fragile)
    assert sorted(m.margin for m in first_real.margins) == sorted(
        m.margin for m in first_control.margins
    )
    assert posterior_margins(matrix.subset(first_real.selection.selected))


def test_oriented_stitching_never_costs_more_than_plain_stitching():
    from phasebudget.simulate import simulate_genomic

    for seed in range(6):
        instance = simulate_genomic(300, coverage=20.0, heterozygosity=0.001,
                                    read_length=10_000, error_rate=0.1, quality_spread=2.0,
                                    seed=seed)
        common = dict(signal="posterior", error_probability=0.05, iterations=3, window=10)
        plain = adaptive_selection(instance.matrix, 3, acceptance="stitched", **common)
        oriented = adaptive_selection(instance.matrix, 3, acceptance="stitched_oriented",
                                      **common)
        assert oriented.full_cost <= plain.full_cost + 1e-9
        assert oriented.selection.selected == plain.selection.selected


def test_a_covered_variant_no_support_read_calls_does_not_break_a_segment():
    # Review of 1 October: after uncalled variants became unphased, `_segments`
    # cut at them. Variant 4 lies inside both support reads with a gap there.
    from phasebudget.adaptive import _segments

    matrix = FragmentMatrix.from_rows(["0000-0", "1111-1", "---000", "---000", "---000",
                                       "---000", "---000"])
    assert _segments(matrix, (0, 1), 3) == [0, 0, 0, 1, 1, 1]
