import random
import statistics

import pytest

from phasebudget import CoverageDP
from phasebudget.simulate import (
    COIN_FLIP_PHRED,
    fragments_from_assignment,
    jittered_error,
    geometric_gap,
    log_odds_weight,
    phred_to_error,
    simulate,
    simulate_connected,
    simulate_genomic,
    simulate_genomic_connected,
    simulate_with_genotype_errors,
)


def test_the_same_seed_gives_the_same_instance():
    first = simulate(12, coverage=8.0, read_span=5, error_rate=0.03, seed=42)
    second = simulate(12, coverage=8.0, read_span=5, error_rate=0.03, seed=42)
    assert first.matrix.to_rows() == second.matrix.to_rows()
    assert first.haplotype == second.haplotype


def test_a_different_seed_gives_a_different_instance():
    first = simulate(20, coverage=8.0, read_span=5, seed=1)
    second = simulate(20, coverage=8.0, read_span=5, seed=2)
    assert first.matrix.to_rows() != second.matrix.to_rows()


def test_coverage_rises_with_the_requested_coverage():
    low = simulate(30, coverage=4.0, read_span=6, seed=3).matrix
    high = simulate(30, coverage=16.0, read_span=6, seed=3).matrix
    assert high.max_coverage() > low.max_coverage()
    assert len(high) > len(low)


def test_reads_do_not_run_past_the_last_variant():
    instance = simulate(9, coverage=12.0, read_span=7, seed=5)
    assert all(fragment.end <= 9 for fragment in instance.matrix)


def test_gaps_appear_only_inside_a_span():
    instance = simulate(25, coverage=10.0, read_span=8, gap_rate=0.4, seed=6)
    gapped = 0
    for fragment in instance.matrix:
        assert fragment.alleles[0] is not None
        assert fragment.alleles[-1] is not None
        gapped += any(allele is None for allele in fragment.alleles)
    assert gapped > 0


def test_error_free_reads_agree_with_the_truth():
    instance = simulate(15, coverage=6.0, read_span=5, error_rate=1e-12, seed=7)
    for index, fragment in enumerate(instance.matrix):
        side = instance.origin[index]
        for position, allele, _ in fragment.observed():
            assert allele == instance.haplotype[position] ^ side


def test_raising_the_error_rate_raises_the_optimal_cost():
    clean = simulate(12, coverage=8.0, read_span=6, error_rate=0.001, seed=8)
    noisy = simulate(12, coverage=8.0, read_span=6, error_rate=0.2, seed=8)
    assert CoverageDP(noisy.matrix).solve().cost > CoverageDP(clean.matrix).solve().cost


def test_uniform_weights_give_unweighted_mec():
    instance = simulate(10, coverage=6.0, read_span=5, seed=9, weight_model="uniform")
    for fragment in instance.matrix:
        assert set(fragment.weights) <= {0.0, 1.0}


def test_simulate_connected_returns_one_block():
    instance = simulate_connected(20, coverage=10.0, read_span=6, seed=10)
    assert set(instance.matrix.blocks()) == {0}


def test_simulate_connected_gives_up_rather_than_lying():
    with pytest.raises(RuntimeError, match="single-block"):
        simulate_connected(60, coverage=0.5, read_span=2, seed=11, attempts=5)


def test_log_odds_weight_grows_as_the_call_gets_certain():
    assert log_odds_weight(0.5) == pytest.approx(0.0)
    assert log_odds_weight(0.001) > log_odds_weight(0.1) > 0.0
    with pytest.raises(ValueError):
        log_odds_weight(0.0)


def test_phred_conversion_matches_the_definition():
    assert phred_to_error(20.0) == pytest.approx(0.01)
    assert phred_to_error(30.0) == pytest.approx(0.001)


def test_hand_built_matrix_is_error_free():
    haplotype = (0, 1, 1, 0, 1)
    matrix = fragments_from_assignment(
        haplotype, [(0, 3), (2, 5), (0, 5)], [0, 1, 1]
    )
    assert CoverageDP(matrix).solve().cost == 0.0


def test_geometric_gaps_are_at_least_one():
    generator = random.Random(11)
    assert all(geometric_gap(generator, 0.01) >= 1 for _ in range(500))
    assert geometric_gap(generator, 1.0) == 1


def test_geometric_gaps_have_the_requested_mean():
    generator = random.Random(12)
    for rate in (0.001, 0.01, 0.25):
        gaps = [geometric_gap(generator, rate) for _ in range(40000)]
        assert statistics.fmean(gaps) == pytest.approx(1.0 / rate, rel=0.05)


def test_a_gap_rate_outside_the_unit_interval_is_an_error():
    generator = random.Random(13)
    for rate in (0.0, -0.1, 1.5):
        with pytest.raises(ValueError, match="rate"):
            geometric_gap(generator, rate)


def test_genomic_positions_are_strictly_increasing():
    instance = simulate_genomic(80, heterozygosity=1e-3, seed=14)
    assert len(instance.positions) == 80
    assert all(b > a for a, b in zip(instance.positions, instance.positions[1:]))
    assert instance.positions[0] == 1


def test_expected_variants_per_read_is_density_times_read_length():
    for heterozygosity, read_length in ((1e-3, 10000), (5e-4, 10000), (1e-3, 4000)):
        lengths = []
        for seed in range(25):
            instance = simulate_genomic(
                400,
                coverage=8.0,
                read_length=read_length,
                heterozygosity=heterozygosity,
                seed=seed,
            )
            lengths.extend(len(fragment) for fragment in instance.matrix)
        # Two biases pull in opposite directions. Reads overhanging either end
        # of the variant range cover fewer sites, and reads covering none are
        # dropped, which size-biases the survivors upward. The second wins at
        # low expected counts, so the band is two-sided.
        expected = heterozygosity * read_length
        assert statistics.fmean(lengths) == pytest.approx(expected, rel=0.15)


def test_genomic_coverage_matches_the_requested_depth():
    depths = []
    for seed in range(25):
        instance = simulate_genomic(
            300, coverage=12.0, read_length=10000, heterozygosity=1e-3, seed=seed
        )
        depths.extend(instance.matrix.coverage())
    assert statistics.fmean(depths) == pytest.approx(12.0, rel=0.05)


def test_genomic_coverage_does_not_fall_off_at_the_ends():
    """Read starts reach back a full read length, so no variant is under-covered."""
    ends: list[int] = []
    middle: list[int] = []
    for seed in range(40):
        instance = simulate_genomic(
            120, coverage=15.0, read_length=10000, heterozygosity=1e-3, seed=seed
        )
        coverage = instance.matrix.coverage()
        ends.extend(coverage[:5] + coverage[-5:])
        middle.extend(coverage[55:65])
    assert statistics.fmean(ends) == pytest.approx(statistics.fmean(middle), rel=0.08)


def test_halving_the_density_halves_the_variants_per_read():
    dense = simulate_genomic(400, heterozygosity=2e-3, read_length=10000, seed=15)
    sparse = simulate_genomic(400, heterozygosity=1e-3, read_length=10000, seed=15)
    dense_mean = statistics.fmean(len(f) for f in dense.matrix)
    sparse_mean = statistics.fmean(len(f) for f in sparse.matrix)
    assert dense_mean > 1.6 * sparse_mean
    # The matrix width is held fixed, so the genome stretches instead.
    assert dense.n_variants == sparse.n_variants == 400
    assert sparse.positions[-1] > 1.6 * dense.positions[-1]


def test_a_genomic_read_covering_one_variant_is_kept():
    """Density low enough that most reads touch a single site, and they stay."""
    instance = simulate_genomic(
        60, coverage=10.0, read_length=1000, heterozygosity=1e-4, seed=16
    )
    assert len(instance.matrix) > 0
    assert any(len(fragment) == 1 for fragment in instance.matrix)


def test_genomic_instances_are_reproducible_from_the_seed():
    first = simulate_genomic(50, heterozygosity=1e-3, error_rate=0.05, seed=17)
    second = simulate_genomic(50, heterozygosity=1e-3, error_rate=0.05, seed=17)
    assert first.matrix.to_rows() == second.matrix.to_rows()
    assert first.haplotype == second.haplotype
    assert first.positions == second.positions


def test_genomic_bad_parameters_are_errors():
    with pytest.raises(ValueError, match="n_variants"):
        simulate_genomic(0)
    with pytest.raises(ValueError, match="read_length"):
        simulate_genomic(10, read_length=0)
    with pytest.raises(ValueError, match="heterozygosity"):
        simulate_genomic(10, heterozygosity=0.0)
    with pytest.raises(ValueError, match="weight model"):
        simulate_genomic(10, weight_model="quadratic")


def test_the_connected_genomic_generator_returns_one_block():
    instance = simulate_genomic_connected(
        40, coverage=10.0, read_length=10000, heterozygosity=1e-3, seed=18
    )
    assert all(block == 0 for block in instance.matrix.blocks())


def test_the_connected_genomic_generator_gives_up_loudly():
    with pytest.raises(RuntimeError, match="no single-block instance"):
        simulate_genomic_connected(
            60,
            coverage=2.0,
            read_length=200,
            heterozygosity=1e-4,
            seed=19,
            attempts=5,
        )


def test_genotype_errors_make_the_two_copies_agree_at_miscalled_sites():
    instance = simulate_with_genotype_errors(
        20, homozygous_rate=0.3, error_rate=0.0, weight_model="uniform", seed=21
    )
    first, second = instance.haplotype, instance.second_haplotype()
    miscalled = set(instance.homozygous)
    assert miscalled
    for index in range(instance.n_variants):
        if index in miscalled:
            assert first[index] == second[index]
        else:
            assert first[index] != second[index]


def test_heterozygous_sites_exclude_the_miscalled_ones():
    instance = simulate_with_genotype_errors(20, homozygous_rate=0.3, seed=22)
    sites = set(instance.heterozygous_sites())
    assert sites.isdisjoint(instance.homozygous)
    assert len(sites) + len(instance.homozygous) == instance.n_variants


def test_a_zero_homozygous_rate_leaves_every_site_heterozygous():
    instance = simulate_with_genotype_errors(20, homozygous_rate=0.0, seed=23)
    assert instance.homozygous == ()
    assert instance.second_haplotype() == instance.complement()
    assert len(instance.heterozygous_sites()) == 20


def test_error_free_reads_match_the_copy_they_came_from():
    instance = simulate_with_genotype_errors(
        16, homozygous_rate=0.25, error_rate=0.0, weight_model="uniform", seed=24
    )
    copies = (instance.haplotype, instance.second_haplotype())
    for fragment, side in zip(instance.matrix, instance.origin):
        for position, allele, _ in fragment.observed():
            assert allele == copies[side][position]


def test_the_miscalled_rate_is_roughly_the_requested_one():
    counts = []
    for seed in range(60):
        instance = simulate_with_genotype_errors(50, homozygous_rate=0.2, seed=seed)
        counts.append(len(instance.homozygous) / 50)
    assert statistics.fmean(counts) == pytest.approx(0.2, rel=0.15)


def test_an_error_free_instance_needs_uniform_weights():
    with pytest.raises(ValueError, match="log-odds"):
        simulate_with_genotype_errors(10, error_rate=0.0)


def test_genotype_error_bad_parameters_are_errors():
    with pytest.raises(ValueError, match="homozygous_rate"):
        simulate_with_genotype_errors(10, homozygous_rate=1.5)
    with pytest.raises(ValueError, match="n_variants"):
        simulate_with_genotype_errors(0)
    with pytest.raises(ValueError, match="weight model"):
        simulate_with_genotype_errors(10, weight_model="cubic")


def test_quality_jitter_never_goes_below_a_coin_flip():
    """A base worse than a coin flip would carry a negative log-odds weight.

    `FragmentMatrix` rejects a negative weight, so an unfloored draw crashes
    the simulator rather than producing a bad instance. The floor is the phred
    score of error probability 0.5.
    """
    generator = random.Random(3)
    drawn = [jittered_error(generator, 0.05, 12.0) for _ in range(20000)]
    assert max(drawn) <= 0.5 + 1e-12
    assert min(drawn) > 0.0
    assert min(log_odds_weight(rate) for rate in drawn) >= 0.0
    assert phred_to_error(COIN_FLIP_PHRED) == pytest.approx(0.5)


def test_quality_jitter_is_off_at_spread_zero():
    generator = random.Random(4)
    assert jittered_error(generator, 0.05, 0.0) == 0.05


def test_a_wide_quality_spread_still_builds_an_instance():
    for spread in (2.0, 4.0, 8.0, 16.0):
        instance = simulate(40, coverage=10.0, read_span=5, error_rate=0.05,
                            seed=11, quality_spread=spread)
        weights = [w for f in instance.matrix for _, _, w in f.observed()]
        assert min(weights) >= 0.0


def test_quality_spread_separates_read_weights():
    """Without it every base weighs the same and the additive objective ties.

    At `quality_spread=0` a read's weight is its called-allele count times one
    constant, so selection is choosing between a handful of distinct values.
    """
    flat = simulate(60, coverage=10.0, read_span=6, error_rate=0.05, seed=12)
    spread = simulate(60, coverage=10.0, read_span=6, error_rate=0.05, seed=12,
                      quality_spread=2.0)
    distinct = lambda matrix: len(
        {round(fragment.total_weight, 9) for fragment in matrix}
    )
    assert distinct(flat.matrix) < 10
    assert distinct(spread.matrix) == len(spread.matrix)


def _fingerprint(instance):
    import hashlib

    # Weights are rounded: with a quality spread they go through gauss, pow and
    # log, whose last bits differ between the libm of Linux, macOS and Windows.
    # The integer structure, which is what a changed random stream would move,
    # is hashed exactly.
    blob = repr((
        instance.matrix.to_rows(),
        [tuple(round(w, 6) for w in f.weights) for f in instance.matrix],
        instance.haplotype,
        instance.origin,
        instance.positions,
    )).encode()
    return hashlib.sha256(blob).hexdigest()[:16]


def test_default_genomic_instances_are_unchanged_by_the_correlated_error_options():
    # Taken before the options existed; they must not move any earlier result.
    from phasebudget.simulate import simulate_genomic

    pinned = {
        (1, ()): "9cd59dcb9c8806c5",
        (7, (("quality_spread", 2.0), ("gap_rate", 0.1))): "b95002c84f5ce5a1",
        (8000, (("error_rate", 0.1), ("read_length", 15000))): "cdfa2279236d37ad",
    }
    for (seed, options), expected in pinned.items():
        instance = simulate_genomic(150, coverage=12, seed=seed, **dict(options))
        assert _fingerprint(instance) == expected
        assert instance.systematic == ()


def _observed_error(instance, sites=None):
    wrong = total = 0
    for fragment, side in zip(instance.matrix, instance.origin):
        for offset, allele in enumerate(fragment.alleles):
            j = fragment.start + offset
            if allele is None or (sites is not None and j not in sites):
                continue
            total += 1
            wrong += allele != instance.haplotype[j] ^ side
    return wrong / total


def test_systematic_sites_pull_both_haplotypes_to_one_allele():
    from phasebudget.simulate import simulate_genomic

    instance = simulate_genomic(
        400, coverage=20, read_length=10000, error_rate=0.01, seed=3,
        systematic_rate=0.1, systematic_strength=0.6,
    )
    sites = set(instance.systematic)
    assert 20 <= len(sites) <= 60
    # Half the reads already carry the biased allele, so the observed error at
    # a systematic site is about strength / 2 plus the base rate.
    assert 0.2 < _observed_error(instance, sites) < 0.4
    others = set(range(instance.n_variants)) - sites
    assert _observed_error(instance, others) < 0.03


def test_bursts_raise_the_true_error_and_leave_the_weights_alone():
    from phasebudget.simulate import simulate_genomic

    plain = simulate_genomic(300, coverage=15, error_rate=0.02, seed=5)
    bursty = simulate_genomic(
        300, coverage=15, error_rate=0.02, seed=5, burst_rate=0.5, burst_length=5,
        burst_error=0.4,
    )
    assert _observed_error(bursty) > 2 * _observed_error(plain)
    weights = {w for f in bursty.matrix for w in f.weights if w > 0}
    assert len(weights) == 1


def test_correlated_error_options_are_validated():
    import pytest as _pytest

    from phasebudget.simulate import simulate_genomic

    for bad in ({"systematic_rate": 1.5}, {"burst_rate": -0.1}, {"burst_length": 0}):
        with _pytest.raises(ValueError):
            simulate_genomic(20, seed=1, **bad)


def test_single_block_instances_cover_every_variant():
    # Review of 1 October: the single-block helpers must keep their original
    # meaning, every variant covered and in block 0, after uncalled variants
    # became unphased.
    from phasebudget.simulate import simulate_connected, simulate_genomic_connected

    for seed in range(20):
        for instance in (simulate_connected(12, coverage=4, read_span=4, seed=seed, gap_rate=0.3),
                         simulate_genomic_connected(40, seed=seed, gap_rate=0.3)):
            assert all(instance.matrix.coverage())
            assert all(block in (0, None) for block in instance.matrix.blocks())
