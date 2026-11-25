"""The 0-1 integer program for weighted minimum error correction.

Variables are `x_r` for the haplotype side of read `r`, `h_j` for the allele of
the first haplotype at variant `j`, and one linearization variable `y_rj` per
observed allele standing for `x_r XOR h_j`. Read `r` carries the value `y_rj`
at variant `j`, so an observation of allele `a` with weight `w` contributes
`w * y_rj` when `a = 0` and `w * (1 - y_rj)` when `a = 1`.

The XOR is linearized by the four inequalities

    y <= x + h,  y >= x - h,  y >= h - x,  y <= 2 - x - h,

which are tight at every integral point. Fixing `h_0 = 0` breaks the complement
symmetry, halving the search tree without removing an optimum.

The model is written in CPLEX LP format with no third-party dependency, so it
can be handed to any solver. `solve` additionally runs it through PuLP when
PuLP is installed.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import _pulp
from ._types import WMECSolution
from .fragments import FragmentMatrix


@dataclass
class ILPModel:
    """A written-out 0-1 program, independent of any solver."""

    matrix: FragmentMatrix
    objective: dict[str, float] = field(default_factory=dict)
    constant: float = 0.0
    constraints: list[tuple[str, dict[str, float], str, float]] = field(
        default_factory=list
    )
    binaries: list[str] = field(default_factory=list)
    fixed: dict[str, int] = field(default_factory=dict)

    @property
    def n_variables(self) -> int:
        return len(self.binaries)

    @property
    def n_constraints(self) -> int:
        return len(self.constraints)

    def to_lp_string(self) -> str:
        lines = [
            "\\ weighted minimum error correction",
            f"\\ {len(self.matrix)} reads, {self.matrix.n_variants} variants",
            "Minimize",
            " obj: " + _linear(self.objective, self.constant),
            "Subject To",
        ]
        for name, terms, sense, rhs in self.constraints:
            lines.append(f" {name}: {_linear(terms)} {sense} {_number(rhs)}")
        for name, value in sorted(self.fixed.items()):
            lines.append(f" fix_{name}: {name} = {value}")
        lines.append("Binary")
        lines.extend(f" {name}" for name in self.binaries)
        lines.append("End")
        return "\n".join(lines) + "\n"


def _number(value: float) -> str:
    text = f"{value:.10g}"
    return text


def _linear(terms: dict[str, float], constant: float = 0.0) -> str:
    parts: list[str] = []
    for name, coefficient in terms.items():
        if coefficient == 0.0:
            continue
        sign = "-" if coefficient < 0 else "+"
        magnitude = _number(abs(coefficient))
        parts.append(f"{sign} {magnitude} {name}")
    if constant:
        sign = "-" if constant < 0 else "+"
        parts.append(f"{sign} {_number(abs(constant))}")
    if not parts:
        return "0"
    text = " ".join(parts)
    return text[2:] if text.startswith("+ ") else text


def build_model(matrix: FragmentMatrix, *, break_symmetry: bool = True) -> ILPModel:
    """Construct the program. No solving, no solver dependency."""
    model = ILPModel(matrix=matrix)
    for index in range(len(matrix)):
        model.binaries.append(f"x_{index}")
    for position in range(matrix.n_variants):
        model.binaries.append(f"h_{position}")
    counter = 0
    for index, fragment in enumerate(matrix):
        for position, allele, weight in fragment.observed():
            name = f"y_{index}_{position}"
            model.binaries.append(name)
            x, h = f"x_{index}", f"h_{position}"
            model.constraints.append(
                (f"c{counter}a", {name: 1.0, x: -1.0, h: -1.0}, "<=", 0.0)
            )
            model.constraints.append(
                (f"c{counter}b", {name: -1.0, x: 1.0, h: -1.0}, "<=", 0.0)
            )
            model.constraints.append(
                (f"c{counter}c", {name: -1.0, x: -1.0, h: 1.0}, "<=", 0.0)
            )
            model.constraints.append(
                (f"c{counter}d", {name: 1.0, x: 1.0, h: 1.0}, "<=", 2.0)
            )
            counter += 1
            if allele == 0:
                model.objective[name] = model.objective.get(name, 0.0) + weight
            else:
                model.objective[name] = model.objective.get(name, 0.0) - weight
                model.constant += weight
    if break_symmetry and matrix.n_variants:
        model.fixed["h_0"] = 0
    return model


def solve(
    matrix: FragmentMatrix,
    *,
    break_symmetry: bool = True,
    solver=None,
) -> WMECSolution:
    """Solve the program with PuLP. Raises `ImportError` when PuLP is absent.

    The returned solution is rebuilt from the `h` variables alone and its cost
    is recomputed from the matrix, so a solver that reports an inconsistent
    objective is caught by `phasebudget.exact.solution_is_optimal` rather than
    trusted.
    """
    try:
        import pulp
    except ImportError as error:  # pragma: no cover - depends on the environment
        raise ImportError(
            "solving the ILP needs PuLP; install with pip install 'phasebudget[ilp]', "
            "or write the model with build_model(...).to_lp_string()"
        ) from error

    model = build_model(matrix, break_symmetry=break_symmetry)
    problem = pulp.LpProblem("wmec", pulp.LpMinimize)
    variables = {
        name: _pulp.variable(problem, name, cat="Binary") for name in model.binaries
    }
    problem += (
        pulp.lpSum(
            coefficient * variables[name]
            for name, coefficient in model.objective.items()
        )
        + model.constant
    )
    for name, terms, sense, rhs in model.constraints:
        expression = pulp.lpSum(
            coefficient * variables[variable] for variable, coefficient in terms.items()
        )
        if sense == "<=":
            problem += expression <= rhs, name
        elif sense == ">=":
            problem += expression >= rhs, name
        else:
            problem += expression == rhs, name
    for name, value in model.fixed.items():
        problem += variables[name] == value, f"fix_{name}"
    status = _pulp.solve(problem, solver)
    if status != "Optimal":
        raise RuntimeError(f"solver returned status {status}")
    haplotype = tuple(
        int(round(variables[f"h_{position}"].value() or 0.0))
        for position in range(matrix.n_variants)
    )
    return WMECSolution(
        cost=matrix.cost(haplotype),
        haplotype=haplotype,
        assignment=matrix.best_assignment(haplotype),
    )
