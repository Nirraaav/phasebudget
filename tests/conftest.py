"""Shared instance generators for the randomized differential tests."""

from __future__ import annotations

import random

import pytest

from phasebudget.fragments import FragmentMatrix


def random_matrix(
    rng: random.Random,
    *,
    max_variants: int = 8,
    max_reads: int = 7,
    weighted: bool = True,
    gap_rate: float = 0.25,
) -> FragmentMatrix:
    """A small matrix with ragged spans, interior gaps, and uneven weights."""
    n_variants = rng.randint(1, max_variants)
    rows = []
    for _ in range(rng.randint(1, max_reads)):
        start = rng.randint(0, n_variants - 1)
        stop = rng.randint(start + 1, n_variants)
        row = ["-"] * n_variants
        for position in range(start, stop):
            interior = stop - start > 2 and position not in (start, stop - 1)
            if interior and rng.random() < gap_rate:
                continue
            row[position] = rng.choice("01")
        if all(symbol == "-" for symbol in row):
            row[start] = rng.choice("01")
        rows.append("".join(row))
    weights = None
    if weighted:
        weights = [
            [round(rng.uniform(0.25, 4.0), 3) for _ in range(n_variants)]
            for _ in rows
        ]
    return FragmentMatrix.from_rows(rows, weights)


@pytest.fixture
def rng() -> random.Random:
    return random.Random(607)
