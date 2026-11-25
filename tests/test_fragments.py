import pytest

from phasebudget import Fragment, FragmentMatrix


def test_from_rows_trims_leading_and_trailing_gaps():
    matrix = FragmentMatrix.from_rows(["--01-", "1-1--"])
    first, second = matrix.fragments
    assert (first.start, first.end) == (2, 4)
    assert (second.start, second.end) == (0, 3)
    assert second.alleles == (1, None, 1)
    assert matrix.n_variants == 5


def test_to_rows_round_trips_the_allele_pattern():
    rows = ("0-1--", "-11--", "---01")
    assert FragmentMatrix.from_rows(list(rows)).to_rows() == rows


def test_coverage_counts_interior_gaps():
    matrix = FragmentMatrix.from_rows(["0--1"])
    assert matrix.coverage() == (1, 1, 1, 1)
    assert matrix.max_coverage() == 1


def test_blocks_split_where_no_read_spans_the_boundary():
    matrix = FragmentMatrix.from_rows(["01--", "--10"])
    assert matrix.blocks() == (0, 0, 1, 1)
    assert matrix.spanned_boundaries() == (True, False, True)


def test_blocks_are_none_at_uncovered_variants():
    matrix = FragmentMatrix.from_rows(["01---", "---01"])
    assert matrix.blocks() == (0, 0, None, 1, 1)


def test_subset_keeps_the_variant_coordinates():
    matrix = FragmentMatrix.from_rows(["0101", "1010", "0011"])
    subset = matrix.subset([0, 2])
    assert subset.n_variants == 4
    assert subset.read_ids == (0, 2)


def test_subset_rejects_unknown_reads():
    matrix = FragmentMatrix.from_rows(["01"])
    with pytest.raises(KeyError):
        matrix.subset([5])


def test_cost_takes_the_cheaper_side_per_read():
    matrix = FragmentMatrix.from_rows(["01", "10"])
    assert matrix.cost((0, 1)) == 0.0
    assert matrix.best_assignment((0, 1)) == (0, 1)


def test_cost_with_a_forced_assignment_is_not_the_minimum():
    matrix = FragmentMatrix.from_rows(["01", "10"])
    assert matrix.cost((0, 1), (0, 0)) == 2.0


def test_weighted_cost_uses_the_given_weights():
    # Against the haplotype pair 00 / 11 the read 01 disagrees with the first
    # side at variant 1, weight 5, and with the second at variant 0, weight 3.
    matrix = FragmentMatrix.from_rows(["01"], [[3.0, 5.0]])
    assert matrix.cost((0, 0)) == 3.0
    assert matrix.cost((0, 0), (0,)) == 5.0


def test_untrimmed_fragment_is_rejected():
    with pytest.raises(ValueError, match="trimmed"):
        FragmentMatrix([Fragment(0, 0, (None, 1), (0.0, 1.0))])


def test_gap_must_carry_zero_weight():
    with pytest.raises(ValueError, match="weight 0.0"):
        FragmentMatrix([Fragment(0, 0, (1, None, 1), (1.0, 2.0, 1.0))])


def test_duplicate_read_ids_are_rejected():
    fragment = Fragment(0, 0, (1,), (1.0,))
    with pytest.raises(ValueError, match="duplicate"):
        FragmentMatrix([fragment, fragment])


def test_row_without_a_called_allele_is_rejected():
    with pytest.raises(ValueError, match="no called allele"):
        FragmentMatrix.from_rows(["---"])


def test_n_variants_smaller_than_the_data_is_rejected():
    fragment = Fragment(0, 0, (1, 0), (1.0, 1.0))
    with pytest.raises(ValueError, match="smaller than"):
        FragmentMatrix([fragment], n_variants=1)


def test_empty_matrix_is_allowed():
    matrix = FragmentMatrix([])
    assert (len(matrix), matrix.n_variants, matrix.max_coverage()) == (0, 0, 0)


def test_a_variant_no_read_calls_is_unphased_but_does_not_break_the_block():
    # Variant 1 lies inside both reads' spans, and both have a gap there.
    matrix = FragmentMatrix.from_rows(["0-0", "1-1"])
    assert matrix.coverage() == (2, 2, 2)
    assert matrix.blocks() == (0, None, 0)
