"""The mpileup parsing and the genotype rule of `benchmarks/grandparent_genotypes.py`."""

import importlib.util
import sys
from pathlib import Path

import pytest

BENCH = Path(__file__).resolve().parent.parent / "benchmarks"


@pytest.fixture(scope="module")
def gg():
    sys.path.insert(0, str(BENCH))
    spec = importlib.util.spec_from_file_location("grandparent_genotypes",
                                                  BENCH / "grandparent_genotypes.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_counts_skip_read_starts_ends_and_indels(gg):
    # Two reference reads, a read start with mapping quality '~', an alt read in both
    # strands, an insertion of 3 bases containing the alt letter, and a deletion marker.
    assert gg.counts("..^~Ta+3TTA,$t*-2AC", "T") == (3, 2)
    assert gg.counts("", "A") == (0, 0)


def test_genotype_rule(gg):
    assert gg.genotype(7, 0) is None
    assert gg.genotype(9, 1) == (0, 0)
    assert gg.genotype(1, 9) == (1, 1)
    assert gg.genotype(5, 5) == (0, 1)
    assert gg.genotype(8, 2) is None  # 0.2: between homozygous and heterozygous
