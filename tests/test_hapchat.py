"""HapCHAT's merge-then-select reduction, as fixed in phasebudget.hapchat."""

import math
import random

from phasebudget.fragments import FragmentMatrix
from phasebudget.hapchat import agreement, log_ratio, merge_reads, select_hapchat
from phasebudget.selection import is_feasible
from phasebudget.simulate import simulate_genomic


def test_agreement_counts_only_shared_called_sites():
    matrix = FragmentMatrix.from_rows(["0110-", "-11-1", "---00"])
    assert agreement(matrix[0], matrix[1]) == (2, 0)
    assert agreement(matrix[0], matrix[2]) == (1, 0)
    assert agreement(matrix[1], matrix[2]) == (0, 1)


def test_the_ratio_follows_the_paper():
    p = 0.15
    for m, x in ((10, 0), (8, 1), (3, 3), (0, 5)):
        expected = math.log(p ** x * (1 - p / 3 ** m)) - math.log(p ** m * (1 - p / 3 ** x))
        assert math.isclose(log_ratio(m, x, p), expected, rel_tol=1e-12)
    # About eight net agreements are needed at p = 0.15 to reach 10^6.
    assert log_ratio(8, 0) >= math.log(1e6) > log_ratio(7, 0)


def test_every_pair_in_a_cluster_passes_and_every_read_is_placed():
    instance = simulate_genomic(
        200, coverage=20, read_length=15000, heterozygosity=0.002, error_rate=0.05, seed=4
    )
    matrix = instance.matrix
    merged = merge_reads(matrix)
    placed = sorted(r for cluster in merged.members for r in cluster)
    assert placed == sorted(matrix.read_ids)
    by_id = {f.read_id: f for f in matrix}
    limit = math.log(1e6)
    for cluster in merged.members:
        for i, a in enumerate(cluster):
            for b in cluster[i + 1:]:
                assert log_ratio(*agreement(by_id[a], by_id[b])) >= limit
    assert len(merged.matrix) < len(matrix)


def test_merged_reads_rarely_mix_haplotypes():
    instance = simulate_genomic(
        200, coverage=20, read_length=15000, heterozygosity=0.002, error_rate=0.05, seed=6
    )
    origin = dict(zip(instance.matrix.read_ids, instance.origin))
    merged = merge_reads(instance.matrix)
    mixed = sum(len({origin[r] for r in cluster}) > 1 for cluster in merged.members)
    assert mixed <= 0.02 * len(merged.members)


def test_the_selection_is_feasible_on_the_merged_matrix():
    rng = random.Random(1)
    for case in range(5):
        instance = simulate_genomic(
            150, coverage=25, read_length=12000, heterozygosity=0.0015,
            error_rate=rng.choice([0.01, 0.1]), seed=40 + case,
        )
        merged, selection = select_hapchat(instance.matrix, 8, seed=case)
        assert is_feasible(merged.matrix, list(selection.selected), 8)
        assert selection.method == "hapchat"


def test_log_ratio_survives_reads_that_share_thousands_of_sites():
    # ONT chromosome 6, 1 October: two ultra-long reads shared more than 646
    # sites, and 3**agree overflowed a float.
    from phasebudget.hapchat import log_ratio

    small = log_ratio(40, 1)
    assert log_ratio(2000, 1) > small
    assert log_ratio(2000, 0) > log_ratio(2000, 1)
