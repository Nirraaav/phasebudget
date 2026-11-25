"""PuLP 2, 3 and 4 behind one small interface.

PuLP 4.0 changed the three calls this package uses: variables are created by
`LpProblem.add_variable` instead of the `LpVariable` constructor,
`PULP_CBC_CMD` and `LpStatus` are gone, and `solve` returns an `LpSolveStats`
whose `status` is an enum. The CI job that installs PuLP picked up 4.0 and
failed on all three. Everything that builds a PuLP program goes through here.
"""

from __future__ import annotations

from typing import Any, Optional


def variable(
    problem: Any,
    name: str,
    *,
    low: Optional[float] = None,
    up: Optional[float] = None,
    cat: str = "Continuous",
) -> Any:
    """A variable of `problem`, created the way the installed PuLP expects."""
    import pulp

    add = getattr(problem, "add_variable", None)
    if add is not None:
        return add(name, lowBound=low, upBound=up, cat=cat)
    return pulp.LpVariable(name, lowBound=low, upBound=up, cat=cat)


#: Tried in order. PuLP up to 3.x bundles CBC as `PULP_CBC_CMD`; PuLP 4 does
#: not, and needs a `cbc`, `highs` or `glpsol` executable on the path.
SOLVER_NAMES = ("PULP_CBC_CMD", "COIN_CMD", "HiGHS_CMD", "GLPK_CMD")


def default_solver() -> Any:
    """The first available solver of `SOLVER_NAMES`, with its log switched off."""
    import pulp

    for name in SOLVER_NAMES:
        factory = getattr(pulp, name, None)
        if factory is None:
            continue
        solver = factory(msg=False)
        if solver.available():
            return solver
    raise RuntimeError(
        "PuLP found no MILP solver; PuLP 4 no longer bundles CBC, so install "
        "one, for example the coinor-cbc package, and put `cbc` on the path"
    )


def solve(problem: Any, solver: Any = None) -> str:
    """Solve and return the status name, `"Optimal"` on success."""
    import pulp

    result = problem.solve(solver or default_solver())
    status = getattr(result, "status", None)
    if status is not None and hasattr(status, "name"):
        return str(status.name)
    return str(pulp.LpStatus[problem.status])
