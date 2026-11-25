"""PhaseBudget: exact read selection and certified phase relations for haplotype phasing."""

from ._types import (
    Fragment,
    Margin,
    Phasing,
    SelectionResult,
    WMECSolution,
)
from ._validation import CoverageLimitExceeded
from .adaptive import (
    MARGIN_SOURCES,
    AdaptiveResult,
    AdaptiveStep,
    adaptive_selection,
    static_selection,
)
from .dp import CoverageDP
from .exact import bruteforce_margins, constrained_optimum, solve_bruteforce
from .fragments import FragmentMatrix
from .metrics import phasing_from_solution, summarize
from .selection import (
    component_count,
    distinct_variants_covered,
    min_coverage,
    read_weights,
    select_fischer_marschall,
    select_greedy,
    select_max_min_coverage,
    select_optimal,
    select_random,
    select_whatshap,
    selection_objective,
)
from .simulate import (
    SimulatedInstance,
    simulate,
    simulate_connected,
    simulate_genomic,
    simulate_genomic_connected,
)

__version__ = "0.1.0"

__all__ = [
    "MARGIN_SOURCES",
    "AdaptiveResult",
    "AdaptiveStep",
    "CoverageDP",
    "CoverageLimitExceeded",
    "Fragment",
    "FragmentMatrix",
    "Margin",
    "Phasing",
    "SelectionResult",
    "SimulatedInstance",
    "WMECSolution",
    "adaptive_selection",
    "bruteforce_margins",
    "component_count",
    "constrained_optimum",
    "distinct_variants_covered",
    "min_coverage",
    "phasing_from_solution",
    "read_weights",
    "select_fischer_marschall",
    "select_greedy",
    "select_max_min_coverage",
    "select_optimal",
    "select_random",
    "select_whatshap",
    "selection_objective",
    "simulate",
    "simulate_connected",
    "simulate_genomic",
    "simulate_genomic_connected",
    "solve_bruteforce",
    "static_selection",
    "summarize",
    "__version__",
]
