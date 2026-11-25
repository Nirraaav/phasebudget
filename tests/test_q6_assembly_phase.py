"""Q6's contig rule, relations and gate arithmetic on fixed and random examples."""

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))

import q6_assembly_phase as q6  # noqa: E402


def brute_contig(intervals, pos0):
    cover = [n for s, e, n in intervals if s <= pos0 < e]
    return cover[0] if len(cover) == 1 else None


def test_segments_fixed():
    segs = q6.segments([(0, 10, "a"), (5, 20, "b"), (30, 40, "c")])
    assert q6.contig_at(segs, 0) == "a"
    assert q6.contig_at(segs, 4) == "a"
    assert q6.contig_at(segs, 5) is None  # a and b overlap
    assert q6.contig_at(segs, 9) is None
    assert q6.contig_at(segs, 10) == "b"
    assert q6.contig_at(segs, 19) == "b"
    assert q6.contig_at(segs, 20) is None  # gap
    assert q6.contig_at(segs, 35) == "c"
    assert q6.contig_at(segs, 40) is None
    assert q6.contig_at(q6.segments([]), 3) is None


def test_segments_random_against_brute_force():
    rng = random.Random(7)
    for _ in range(300):
        intervals = []
        for k in range(rng.randint(0, 6)):
            s = rng.randint(0, 50)
            intervals.append((s, s + rng.randint(1, 20), f"c{k}"))
        segs = q6.segments(intervals)
        for pos in range(-2, 75):
            assert q6.contig_at(segs, pos) == brute_contig(intervals, pos)


def test_contig_of_needs_the_whole_span():
    segs = q6.segments([(0, 10, "a"), (10, 20, "b")])
    assert q6.contig_of(segs, 5, "A") == "a"  # 0-based 4
    assert q6.contig_of(segs, 10, "A") == "a"  # 0-based 9, last base of a
    assert q6.contig_of(segs, 10, "AC") is None  # spans 9 and 10: two contigs
    assert q6.contig_of(segs, 11, "ACG") == "b"


def test_assembly_values():
    variants = [(5, "A", "C"), (6, "A", "C"), (7, "A", "C"), (15, "A", "C"), (8, "AT", "A")]
    calls = {(5, "A", "C"): "0|1", (6, "A", "C"): "1|0", (7, "A", "C"): "1|1",
             (15, "A", "C"): "0|1", (8, "AT", "A"): "1|0"}
    hap1 = q6.segments([(0, 10, "h1a"), (10, 30, "h1b")])
    hap2 = q6.segments([(0, 30, "h2a")])
    out = q6.assembly_values(variants, calls, lambda p, r: p != 6, hap1, hap2)
    assert out[0] == [0, "h1a|h2a"]
    assert out[1] is None  # outside the BED
    assert out[2] is None  # homozygous
    assert out[3] == [0, "h1b|h2a"]
    assert out[4] == [1, "h1a|h2a"]
    # different hap1 contigs relate nothing
    assert q6.relation(out[0], out[3]) is None
    assert q6.relation(out[0], out[4]) == 0


def test_relation():
    assert q6.relation(None, 1) is None
    assert q6.relation(1, 1) == 1
    assert q6.relation(0, 1) == 0
    assert q6.relation([1, "x"], [1, "x"]) == 1
    assert q6.relation([1, "x"], [0, "x"]) == 0
    assert q6.relation([1, "x"], [1, "y"]) is None
    assert q6.relation([1, None], [1, None]) is None


def pair(group, ws, lp, asm, v5, trio=None):
    return ("chr2", 1, 2, group, {"whatshap": ws, "longphase": lp, "assembly": asm}, v5, trio)


def test_gates():
    pairs = ([pair("segdup", 0, 1, 1, 1)] * 17 + [pair("segdup", 0, 1, 0, 1)] * 3
             + [pair("lowmap", 0, 1, 0, 1)] * 5 + [pair("tandem", 1, 1, None, 1)] * 30)
    ra = q6.gate_a(pairs, 5)
    assert ra == {"segdup": (20, 17), "lowmap": (5, 0), "tandem": (0, 0)}
    assert q6.a_holds(ra) == (True, ["segdup"])  # 17/20 = 0.85; lowmap not judged
    rb = q6.gate_b(pairs, 5)
    assert rb["segdup"] == (20, 20, 0, 3)
    assert rb["tandem"] == (0, 0, 0, 0)  # the assembly does not phase them
    assert q6.b_holds(rb) == (True, ["segdup"])  # 3 * 3 <= 20
    edge = pairs + [pair("segdup", 0, 1, 0, 1)] * 5
    assert q6.b_holds(q6.gate_b(edge, 5))[0] is True  # 3 * 8 <= 25
    worse = pairs + [pair("segdup", 0, 1, 0, 1)] * 6
    assert q6.a_holds(q6.gate_a(worse, 5))[0] is False  # 17/26 < 0.8
    assert q6.b_holds(q6.gate_b(worse, 5))[0] is False  # 3 * 9 > 26
    assert q6.a_holds({"segdup": (3, 3)}) == (False, [])  # nothing judged
