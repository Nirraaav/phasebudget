"""Eligibility and reasons of `benchmarks/r3_label_eligible.py` (R3)."""

import importlib.util
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parent.parent / "benchmarks"


def _module():
    sys.path.insert(0, str(BENCH))
    spec = importlib.util.spec_from_file_location("r3_label_eligible",
                                                  BENCH / "r3_label_eligible.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_pair_is_eligible_only_when_both_sites_are():
    m = _module()
    assert m.pair_class("eligible", "eligible") == "eligible"
    assert m.pair_class("eligible", "other") == "other"
    assert m.pair_class("other", "uninformative") == "uninformative"


def test_site_is_uninformative_only_when_both_replicates_say_so():
    m = _module()
    key = ("chr1", 10, "A", "G")
    het = {key: 1}
    hom = {key: 2}
    assert m.site_reason(*key, {"A": (het, het), "B": (het, het)}) == "uninformative"
    assert m.site_reason(*key, {"A": (het, het), "B": (hom, het)}) == "other"
    assert m.site_reason(*key, {"A": ({}, het), "B": (het, het)}) == "other"  # missing reads 0/0
