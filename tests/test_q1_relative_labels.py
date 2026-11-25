"""Sibling IBD labels of `benchmarks/q1_relative_labels.py` (Q1)."""

import importlib.util
import itertools
import math
import random
import sys
from pathlib import Path

import pytest

BENCH = Path(__file__).resolve().parent.parent / "benchmarks"


def _module():
    sys.path.insert(0, str(BENCH))
    spec = importlib.util.spec_from_file_location("q1_relative_labels",
                                                  BENCH / "q1_relative_labels.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_category_of_genotype_pairs():
    m = _module()
    assert m.category(0, 0) is None
    assert m.category(-1, 1) is None
    assert m.category(0, 2) == m.category(2, 0) == 0
    assert m.category(1, 0) == m.category(2, 1) == m.category(0, 1) == 1
    assert m.category(1, 1) == 2
    assert m.category(2, 2) == 3


def test_transitions_are_stochastic_and_reach_the_sibling_prior():
    m = _module()
    for d in (0, 1_000, 1e6, 1e8):
        for row in m.transition(d):
            assert math.isclose(sum(row), 1.0)
    assert m.transition(0) == ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    for row in m.transition(1e11):
        assert all(math.isclose(x, p, abs_tol=1e-9) for x, p in zip(row, m.INIT))


def test_forward_backward_matches_path_enumeration():
    """Posterior marginals against summing over all 3^n state paths."""
    m = _module()
    rng = random.Random(3)
    emit = m.EMIT0
    for _ in range(5):
        n = 5
        pos = sorted(rng.sample(range(10_000_000), n))
        obs = [rng.randrange(4) for _ in range(n)]
        marg = [[0.0] * 3 for _ in range(n)]
        for path in itertools.product(range(3), repeat=n):
            w = m.INIT[path[0]] * emit[path[0]][obs[0]]
            for t in range(1, n):
                w *= m.transition(pos[t] - pos[t - 1])[path[t - 1]][path[t]] * emit[path[t]][obs[t]]
            for t, s in enumerate(path):
                marg[t][s] += w
        post = m.forward_backward(pos, obs, emit)
        for t in range(n):
            z = sum(marg[t])
            for s in range(3):
                assert post[s][t] == pytest.approx(marg[t][s] / z, rel=1e-9, abs=1e-12)


def test_site_labels_need_a_confident_run_and_the_margin():
    m = _module()
    pos = [0, 100, 200, 300, 400, 500, 600]
    ga = [1, 1, 1, 1, 1, 1, 1]
    gb = [0, 2, 0, 1, 2, 0, 2]
    p1 = [0.999, 0.999, 0.999, 0.999, 0.5, 0.999, 0.999]
    labels = m.site_labels(pos, ga, gb, p1, tau=0.99, margin=0)
    assert labels == {0: (0, 0), 100: (1, 0), 200: (0, 0), 500: (0, 1), 600: (1, 1)}
    trimmed = m.site_labels(pos, ga, gb, p1, tau=0.99, margin=100)
    assert trimmed == {100: (1, 0), 200: (0, 0)}  # the second run is 100 bp long
    assert m.sibling_relation(labels, 0, 200) is True
    assert m.sibling_relation(labels, 0, 100) is False
    assert m.sibling_relation(labels, 200, 500) is None  # different runs
    assert m.sibling_relation(labels, 0, 300) is None


def test_both_needs_both_labels():
    m = _module()
    hmm = {0: (0, 0), 10: (1, 0), 30: (0, 1)}
    window = {10: (1, 4), 30: (0, 4), 40: (1, 4)}
    assert m.both(hmm, window) == {10: (1, (0, 4)), 30: (0, (1, 4))}
    assert m.sibling_relation(m.both(hmm, window), 10, 30) is None  # different HMM runs


def test_siblings_are_combined_only_when_they_agree():
    m = _module()
    assert m.combined([None, True, True]) is True
    assert m.combined([False, None]) is False
    assert m.combined([True, False]) is None
    assert m.combined([None, None]) is None


def test_window_labels_fixed():
    m = _module()
    # A het everywhere; B: hom at 0, 10, 30, 1000, 1010; opposite homozygote (IBS0) at 500
    pos = [0, 10, 20, 30, 500, 1000, 1010, 1020]
    ga = [1, 1, 1, 1, 0, 1, 1, 1]
    gb = [0, 2, 1, 0, 2, 2, 0, 1]
    assert m.category(0, 2) == 0
    labels = m.window_labels(pos, ga, gb, w=100, min_ibs1=0.25)
    assert labels == {0: (0, 0), 10: (1, 0), 30: (0, 0), 1000: (1, 1), 1010: (0, 1)}
    assert m.sibling_relation(labels, 0, 10) is False
    assert m.sibling_relation(labels, 0, 30) is True
    assert m.sibling_relation(labels, 30, 1000) is None  # an IBS0 site lies between
    # a window reaching the IBS0 site (at 500) labels nothing; 1010's window starts at 510
    assert m.window_labels(pos, ga, gb, w=500, min_ibs1=0.25) == {1010: (0, 1)}
    # IBD2-like: every site both heterozygous or both hom alt, one stray IBS1 site
    pos2 = list(range(0, 1000, 10))
    ga2 = [1] * 100
    gb2 = [1] * 99 + [2]
    assert m.window_labels(pos2, ga2, gb2, w=200, min_ibs1=0.25) == {}


@pytest.mark.parametrize("labeller", ["both", "window", "hmm"])
def test_simulated_siblings_give_correct_labels(labeller):
    """Randomized: exact genotypes on simulated sibships label about an eighth of the child's
    pairs and never wrongly (the HMM and both rules), or wrongly in under 3e-3 of labels (the
    window rule alone).
    The simulation draws mostly rare alleles, so a stretch shared by neither parent can lack
    opposite homozygotes over 250 kb; in the Pedigree data 7.7% of informative sites in such
    stretches are opposite homozygotes (amendment 1)."""
    m = _module()
    total = labelled = wrong = 0
    for seed in range(4):
        n, k, w = m.simulate_check(120_000, 300_000_000, seed, labeller=labeller)
        total, labelled, wrong = total + n, labelled + k, wrong + w
    assert wrong == 0 if labeller != "window" else wrong < 3e-3 * labelled
    assert 0.03 < labelled / total < 0.3


def test_genotype_errors_reach_the_labels():
    """A 1% genotype error rate gives label errors of the order of the rate: the HMM cannot
    see a single wrong homozygous call inside an IBD1 run."""
    m = _module()
    n, k, w = m.simulate_check(120_000, 300_000_000, 1, rate=1e-2, labeller="hmm")
    assert k > 0 and 2e-3 < w / k < 0.1
