"""The port of WhatsHap's read selection, against WhatsHap itself.

`select_whatshap` claims to select exactly the reads that WhatsHap 2.8's
`readselection` selects. `tests/fixtures/whatshap_readselection.json` holds
WhatsHap's own output on 160 simulated instance and cap pairs, recorded with
`benchmarks/whatshap_bridge.py`, so the claim is checked without WhatsHap
installed. When WhatsHap is importable, a live comparison runs as well.
"""

import importlib.util
import json
import random
from pathlib import Path

import pytest

from phasebudget import FragmentMatrix
from phasebudget.selection import (
    _StdIntSet,
    _WhatsHapHeap,
    is_feasible,
    is_phase_informative,
    select_fischer_marschall,
    select_whatshap,
    whatshap_score_drop,
)
from phasebudget.simulate import simulate

from conftest import random_matrix

ROOT = Path(__file__).parents[1]


def _bridge():
    path = ROOT / "benchmarks" / "whatshap_bridge.py"
    spec = importlib.util.spec_from_file_location("whatshap_bridge", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _recorded():
    return json.loads((ROOT / "tests" / "fixtures" / "whatshap_readselection.json").read_text())


def _instance(case, bridge):
    params = {key: case[key] for key in
              ("coverage", "read_span", "error_rate", "quality_spread", "gap_rate", "seed")}
    instance = simulate(case["n_variants"], **params)
    return instance, bridge.integer_weights(instance.matrix)


def test_the_port_reproduces_recorded_whatshap_selections():
    bridge = _bridge()
    data = _recorded()
    assert data["whatshap_version"] == "2.8"
    mismatched = []
    for case in data["cases"]:
        instance, matrix = _instance(case, bridge)
        selected = select_whatshap(matrix, case["cap"], positions=instance.positions).selected
        if list(selected) != case["selected"]:
            mismatched.append((case["seed"], case["cap"]))
    assert mismatched == []


def test_the_published_pseudocode_does_not_reproduce_them():
    # This is the finding that made the port necessary: the paper's procedure
    # and the deployed code choose different reads on most instances.
    bridge = _bridge()
    cases = _recorded()["cases"]
    same = 0
    for case in cases:
        _, matrix = _instance(case, bridge)
        same += list(select_fischer_marschall(matrix, case["cap"]).selected) == case["selected"]
    assert same < len(cases) // 2


def test_live_agreement_with_whatshap():
    pytest.importorskip("whatshap")
    bridge = _bridge()
    for seed in range(20):
        instance = simulate(40, coverage=8, read_span=4, error_rate=0.1,
                            quality_spread=2.0, gap_rate=0.1, seed=500 + seed)
        matrix = bridge.integer_weights(instance.matrix)
        informative = bridge.phase_informative(matrix)
        for cap in (2, 4):
            expected = bridge.whatshap_selection(informative, cap, instance.positions)
            assert select_whatshap(matrix, cap, positions=instance.positions).selected == expected


def test_live_dp_cost_agreement_with_whatshap():
    pytest.importorskip("whatshap")
    from phasebudget.dp import CoverageDP

    bridge = _bridge()
    for seed in range(10):
        instance = simulate(25, coverage=5, read_span=4, error_rate=0.1,
                            quality_spread=2.0, gap_rate=0.1, seed=900 + seed)
        matrix = bridge.phase_informative(bridge.integer_weights(instance.matrix))
        assert CoverageDP(matrix).solve().cost == pytest.approx(bridge.whatshap_cost(matrix))


def test_the_score_drop_counts_variants_not_newly_covered():
    # A queued read calling variants 1..5, after an admitted read newly
    # covered 0 and 1: WhatsHap subtracts 4, the four variants outside {0, 1}.
    assert whatshap_score_drop([1, 2, 3, 4, 5], {0, 1}) == 4
    assert whatshap_score_drop([0, 1], {0, 1}) == 0


def test_ties_follow_the_heap_layout_not_the_read_order():
    # Reads 0 and 1 tie after read 2 is admitted. WhatsHap's heap then yields
    # read 1; a read-id tie-break would yield read 0.
    matrix = FragmentMatrix.from_rows(
        ["-10", "-01", "10-"],
        weights=[[0, 32, 16], [0, 16, 27], [19, 17, 0]],
    )
    assert select_whatshap(matrix, 2).selected == (1, 2)
    assert select_fischer_marschall(matrix, 2).selected == (0, 2)


def test_the_heap_pops_in_non_increasing_order_under_score_changes():
    rng = random.Random(4)
    heap = _WhatsHapHeap()
    scores = {}
    for item in range(60):
        scores[item] = (rng.randint(0, 5), rng.randint(0, 5), rng.randint(0, 3))
        heap.push(scores[item], item)
    for item in rng.sample(range(60), 30):
        scores[item] = (scores[item][0] - rng.randint(0, 3),) + scores[item][1:]
        heap.change_score(item, scores[item])
    popped = []
    while heap:
        score, item = heap.pop()
        assert score == scores[item]
        popped.append(score)
    assert popped == sorted(popped, reverse=True)


def test_unordered_set_order_matches_libstdcxx():
    # Recorded from g++ 13.3 / libstdc++: bucket count 13 after the first
    # insertion, a new element at the front of its bucket or of the list.
    values = _StdIntSet()
    for value in (30, 70, 10, 50, 120, 40, 80):
        values.insert(value)
    assert values.order == [80, 40, 120, 50, 10, 70, 30]
    values.clear()
    for value in (60, 20):
        values.insert(value)
    assert values.order == [20, 60]


def test_the_port_respects_the_cap_and_selects_informative_reads():
    rng = random.Random(11)
    for _ in range(300):
        matrix = random_matrix(rng, max_variants=9, max_reads=9)
        for cap in (1, 2, 3):
            selected = select_whatshap(matrix, cap).selected
            assert is_feasible(matrix, selected, cap)
            by_id = {fragment.read_id: fragment for fragment in matrix}
            assert all(is_phase_informative(by_id[read_id]) for read_id in selected)


def test_positions_must_match_the_variants():
    matrix = FragmentMatrix.from_rows(["01", "10"])
    with pytest.raises(ValueError):
        select_whatshap(matrix, 1, positions=[5])
