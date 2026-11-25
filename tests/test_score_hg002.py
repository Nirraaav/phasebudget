"""The HG002 scorer's parsers and its sorted AUC, on hand-written inputs."""

import importlib.util
import random
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]


def _load(name):
    path = ROOT / "benchmarks" / name
    spec = importlib.util.spec_from_file_location(name[:-3], path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


VARIANTS = [(100, "A", "G"), (200, "C", "T"), (300, "G", "A"), (400, "T", "C")]
HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"


def test_parental_truth_keeps_only_parent_oriented_phased_calls(tmp_path):
    score = _load("score_hg002.py")
    vcf = tmp_path / "truth.vcf"
    vcf.write_text(HEADER
                   + "chr20\t100\t.\tA\tG\t.\tPASS\t.\tGT:PS\t1|0:PATMAT\n"
                   + "chr20\t200\t.\tC\tT\t.\tPASS\t.\tGT:PS\t0|1:1234\n"     # local phase set
                   + "chr20\t300\t.\tG\tC\t.\tPASS\t.\tGT:PS\t0|1:PATMAT\n"   # other ALT
                   + "chr20\t400\t.\tT\tC\t.\tPASS\t.\tGT:PS\t0/1:PATMAT\n")  # unphased
    assert score.parental_truth(vcf, VARIANTS) == [1, None, None, None]


def test_whatshap_output_becomes_alleles_and_blocks(tmp_path):
    score = _load("score_hg002.py")
    vcf = tmp_path / "phased.vcf"
    vcf.write_text(HEADER
                   + "chr20\t100\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:100\n"
                   + "chr20\t200\t.\tC\tT\t.\tPASS\t.\tGT:PS\t1|0:100\n"
                   + "chr20\t300\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\n"
                   + "chr20\t400\t.\tT\tC\t.\tPASS\t.\tGT:PS\t1|0:400\n")
    phasing = score.whatshap_phasing(vcf, VARIANTS)
    assert phasing.alleles == (0, 1, None, 1)
    assert phasing.blocks == (100, 100, None, 400)


def test_hapcut2_blocks_become_alleles_and_blocks(tmp_path):
    score = _load("score_hg002.py")
    blocks = tmp_path / "out.blocks"
    blocks.write_text(
        "BLOCK: offset: 1 len: 2 phased: 2 SPAN: 101 fragments 3\n"
        "1\t1\t0\tchr20\t100\tA\tG\t0/1\t0\t.\t100.00\t10\n"
        "2\t0\t1\tchr20\t200\tC\tT\t0/1\t0\t.\t100.00\t10\n"
        "******** \n"
        "BLOCK: offset: 3 len: 2 phased: 1 SPAN: 1 fragments 1\n"
        "3\t-\t-\tchr20\t300\tG\tA\t0/1\t0\t.\t0.00\t0\n"
        "4\t0\t1\tchr20\t400\tT\tC\t0/1\t0\t.\t100.00\t10\n"
        "******** \n")
    phasing = score.hapcut2_phasing(blocks, VARIANTS)
    assert phasing.alleles == (1, 0, None, 0)
    assert phasing.blocks == (0, 0, None, 1)


def test_sorted_auc_matches_the_pairwise_definition():
    score = _load("score_hg002.py")
    rank_probability = _load("benchmark_objective_truth.py").rank_probability
    rng = random.Random(2)
    for _ in range(50):
        bad = [rng.choice([0.0, 1.0, 2.5, rng.random() * 5]) for _ in range(rng.randint(1, 15))]
        good = [rng.choice([0.0, 1.0, 2.5, rng.random() * 5]) for _ in range(rng.randint(1, 30))]
        assert score.auc(bad, good) == pytest.approx(rank_probability(bad, good))


def test_time_logs_are_summed_and_maximized(tmp_path):
    score = _load("score_hg002.py")
    log = tmp_path / "run.log"
    log.write_text("\tElapsed (wall clock) time (h:mm:ss or m:ss): 1:02.50\n"
                   "\tMaximum resident set size (kbytes): 2048\n"
                   "\tElapsed (wall clock) time (h:mm:ss or m:ss): 1:00:00\n"
                   "\tMaximum resident set size (kbytes): 1024\n")
    seconds, peak = score.time_log(log)
    assert seconds == pytest.approx(62.5 + 3600)
    assert peak == pytest.approx(2.0)
