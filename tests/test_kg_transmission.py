"""The trio rule of `genome_common.transmitted` (W1 amendment 1)."""

import importlib.util
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parent.parent / "benchmarks"


def _module():
    sys.path.insert(0, str(BENCH))
    spec = importlib.util.spec_from_file_location("genome_common", BENCH / "genome_common.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_transmitted_follows_the_homozygous_parent():
    t = _module().transmitted
    assert t(1, 0, 1) == 0 and t(1, 2, 1) == 1  # father homozygous
    assert t(1, 1, 0) == 1 and t(1, 1, 2) == 0  # mother homozygous
    assert t(1, 0, 2) == 0 and t(1, 2, 0) == 1  # both, opposite


def test_transmitted_refuses_what_it_cannot_decide():
    t = _module().transmitted
    assert t(1, 1, 1) is None  # both heterozygous
    assert t(1, 0, 0) is None and t(1, 2, 2) is None  # not Mendelian
    assert t(0, 0, 1) is None and t(2, 2, 1) is None  # child not heterozygous
    assert t(1, None, 0) is None and t(None, 0, 2) is None


def test_alt_count_reads_genotypes_only():
    a = _module()._alt_count
    assert a("0|1") == 1 and a("1/1:30") == 2 and a("0|0") == 0
    assert a("./.") is None and a("1|2") is None and a("1") is None


def test_depth_seed_reads_seeded_and_plain_labels():
    d = _module().depth_seed
    assert d("20") == (20.0, 0) and d("10s3") == (10.0, 3) and d("30s12") == (30.0, 12)
