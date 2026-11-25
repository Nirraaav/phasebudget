import random

import pytest

from phasebudget import FragmentMatrix
from phasebudget.exact import solve_bruteforce
from phasebudget.ilp import build_model, solve

from conftest import random_matrix



def test_the_model_has_one_variable_per_read_variant_and_observation():
    matrix = FragmentMatrix.from_rows(["011", "0-1"])
    model = build_model(matrix)
    observations = sum(fragment.n_observed for fragment in matrix)
    assert model.n_variables == len(matrix) + matrix.n_variants + observations
    assert model.n_constraints == 4 * observations


def test_the_symmetry_breaking_constraint_is_present_by_default():
    model = build_model(FragmentMatrix.from_rows(["01"]))
    assert model.fixed == {"h_0": 0}
    assert "fix_h_0: h_0 = 0" in model.to_lp_string()
    assert build_model(FragmentMatrix.from_rows(["01"]), break_symmetry=False).fixed == {}


def test_the_lp_file_is_well_formed():
    text = build_model(FragmentMatrix.from_rows(["011", "100"])).to_lp_string()
    assert text.startswith("\\ weighted minimum error correction")
    assert "Minimize" in text
    assert "Subject To" in text
    assert text.rstrip().endswith("End")
    assert "Binary" in text


def test_the_objective_carries_the_weights():
    matrix = FragmentMatrix.from_rows(["01"], [[2.0, 3.0]])
    model = build_model(matrix)
    assert model.objective["y_0_0"] == 2.0
    assert model.objective["y_0_1"] == -3.0
    assert model.constant == 3.0


def test_the_objective_at_the_optimum_matches_the_brute_force_cost():
    # Evaluate the written-out program by hand at the brute-force optimum,
    # which checks the linearization without needing a solver.
    generator = random.Random(41)
    for _ in range(60):
        matrix = random_matrix(generator, max_variants=6, max_reads=5)
        optimum = solve_bruteforce(matrix)
        model = build_model(matrix)
        assignment = {}
        for index in range(len(matrix)):
            assignment[f"x_{index}"] = optimum.assignment[index]
        for position in range(matrix.n_variants):
            assignment[f"h_{position}"] = optimum.haplotype[position]
        for index, fragment in enumerate(matrix):
            for position, _, _ in fragment.observed():
                assignment[f"y_{index}_{position}"] = (
                    assignment[f"x_{index}"] ^ assignment[f"h_{position}"]
                )
        value = model.constant + sum(
            coefficient * assignment[name]
            for name, coefficient in model.objective.items()
        )
        assert value == pytest.approx(optimum.cost)
        for name, terms, sense, rhs in model.constraints:
            total = sum(
                coefficient * assignment[variable]
                for variable, coefficient in terms.items()
            )
            assert total <= rhs + 1e-9, name


def test_solve_matches_brute_force_when_pulp_is_available():
    pytest.importorskip("pulp")
    generator = random.Random(43)
    for _ in range(10):
        matrix = random_matrix(generator, max_variants=6, max_reads=5)
        assert solve(matrix).cost == pytest.approx(solve_bruteforce(matrix).cost)


def test_solve_says_what_to_install_when_pulp_is_missing(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name == "pulp":
            raise ImportError("no pulp")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    with pytest.raises(ImportError, match="phasebudget\\[ilp\\]"):
        solve(FragmentMatrix.from_rows(["01"]))
