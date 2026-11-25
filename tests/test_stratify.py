"""Stratum membership in `benchmarks/stratify.py`."""

import gzip
import importlib.util
import sys
from pathlib import Path


def _module():
    path = Path(__file__).parents[1] / "benchmarks" / "stratify.py"
    spec = importlib.util.spec_from_file_location("stratify", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["stratify"] = module
    spec.loader.exec_module(module)
    return module


def test_bed_is_merged_and_half_open(tmp_path):
    stratify = _module()
    bed = tmp_path / "x.bed.gz"
    with gzip.open(bed, "wt") as handle:
        handle.write("chr1\t10\t20\nchr1\t15\t30\nchr1\t40\t41\nchr2\t0\t100\n")
    beds = stratify.load_bed(bed, ["chr1"])
    assert beds["chr1"] == ([10, 40], [30, 41])
    assert stratify.inside(beds["chr1"], [9, 10, 29, 30, 40, 41]) == \
        [False, True, True, False, True, False]
