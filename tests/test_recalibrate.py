"""Base-quality recalibration: counting against the truth, pooling bins, weights."""

import importlib.util
import math
import random
from collections import Counter
from pathlib import Path

import pytest

from phasebudget._types import Fragment
from phasebudget.fragments import FragmentMatrix

ROOT = Path(__file__).parents[1]


def _load():
    path = ROOT / "benchmarks" / "recalibrate.py"
    spec = importlib.util.spec_from_file_location("recalibrate", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _simulated(rng, truth, rates, n_reads=4000, length=8):
    """Reads of known haplotype whose allele errors follow `rates` per quality."""
    fragments = []
    qualities = list(rates)
    for read_id in range(n_reads):
        start = rng.randrange(len(truth) - length)
        side = rng.randint(0, 1)
        alleles, weights = [], []
        for j in range(start, start + length):
            q = rng.choice(qualities)
            wrong = rng.random() < rates[q]
            alleles.append(truth[j] ^ side ^ int(wrong))
            weights.append(float(q))
        fragments.append(Fragment(read_id, start, tuple(alleles), tuple(weights)))
    return FragmentMatrix(fragments, len(truth))


def test_counts_recover_the_simulated_error_rates():
    module = _load()
    rng = random.Random(3)
    truth = [rng.randint(0, 1) for _ in range(500)]
    rates = {20: 0.05, 40: 0.01, 93: 0.001}
    alleles, wrong = module.count_errors(_simulated(rng, truth, rates), truth)
    for q, rate in rates.items():
        assert wrong[q] / alleles[q] == pytest.approx(rate, abs=0.35 * rate + 0.001)


def test_reads_with_few_truth_sites_or_a_tie_are_skipped():
    module = _load()
    truth = [0, 1, None, None, 0, 1]
    matrix = FragmentMatrix([
        Fragment(0, 0, (0, 1), (30.0, 30.0)),          # two truth sites: skipped
        Fragment(1, 0, (0, 0, 0, 0), (30.0,) * 4),       # sites 0, 1 only: skipped
        Fragment(2, 0, (0, 1, 0, 0, 0, 1), (30.0,) * 6),  # four sites, all agree
        Fragment(3, 0, (0, 0, 1, 1, 1, 1), (30.0,) * 6),  # two agree, two disagree: tie
    ], 6)
    alleles, wrong = module.count_errors(matrix, truth)
    assert alleles == Counter({30: 4}) and wrong == Counter()


def test_bins_pool_from_the_top_until_fifty_wrong():
    module = _load()
    alleles = Counter({93: 10000, 80: 5000, 60: 2000, 30: 500})
    wrong = Counter({93: 20, 80: 40, 60: 60, 30: 10})
    table = module.fit_table(alleles, wrong)
    assert [(b["high"], b["low"]) for b in table] == [(93, 80), (60, 30)]
    assert table[0]["wrong"] == 60 and table[1]["wrong"] == 70
    first = table[0]
    rate = (60 + 1) / (15000 + 2)
    assert first["error"] == pytest.approx(rate)
    assert first["weight"] == pytest.approx(math.log((1 - rate) / rate))


def test_apply_table_maps_every_allele_and_keeps_gaps_at_zero():
    module = _load()
    table = [{"high": 93, "low": 60, "weight": 7.0}, {"high": 59, "low": 3, "weight": 3.0}]
    matrix = FragmentMatrix([Fragment(0, 0, (0, None, 1), (93.0, 0.0, 20.0))], 3)
    assert module.apply_table(matrix, table)[0].weights == (7.0, 0.0, 3.0)
