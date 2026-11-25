"""Score the registered ONT prediction: R1, R2 and R3 of `docs/design/ont-prediction.md`.

Reads committed files and nothing else: the baselines B0 to B2 (amendment 3),
the ONT scoring's O7 rows, and the per-chromosome prediction sums, pooled
here at full precision. The printed prediction has two decimals, and a
rounded value can cross R2's threshold of 0.5.

Usage: `python3 benchmarks/score_prediction.py BASELINES ONT_SCORE PREDICTION_TSV...`, or
`python3 benchmarks/score_prediction.py --hg005 - HG005_SCORE PREDICTION_TSV` for H5-3,
where R3 is not registered and its lines are to be ignored.
"""

from __future__ import annotations

import csv
import re
import sys

DEPTHS = ("full", "30", "20", "10")
CAPS = (3, 5, 10, 15)
CELL = re.compile(r"(\w+)x ([+-][\d.]+)(?: \[([+-][\d.]+), ([+-][\d.]+)\])?")


def rows(text: str, header: str) -> dict:
    """{(cap, depth): (value, low, high)} from the `cap N:` lines after `header`."""
    out, inside = {}, False
    for line in text.splitlines():
        if line.startswith(header):
            inside = True
            continue
        if inside:
            m = re.match(r"\s+cap\s+(\d+): (.*)", line)
            if not m:
                break
            for depth, value, low, high in CELL.findall(m.group(2)):
                out[(int(m.group(1)), depth)] = (float(value), low and float(low),
                                                 high and float(high))
    return out


def sign(x: float) -> int:
    return (x > 0) - (x < 0)


def within_two(p: float, o: float) -> bool:
    return p * o > 0 and 0.5 <= p / o <= 2


def score(prediction: dict, baselines: dict, observed: dict) -> list[str]:
    out = ["cell        observed [interval]       model    B0     B1     B2"]
    for c in CAPS:
        for d in DEPTHS:
            o, lo, hi = observed[(c, d)]
            out.append(f"cap {c:>2} {d:>4}x {o:+7.2f} [{lo:+.2f}, {hi:+.2f}]  "
                       f"{prediction[(c, d)]:+8.4f} " + " ".join(
                           f"{baselines[b][(c, d)]:+6.2f}" for b in ("B0", "B1", "B2")))
    misses = [f"c{c} {d}x" for c in (3, 5) for d in DEPTHS
              if not (observed[(c, d)][1] <= prediction[(c, d)] <= observed[(c, d)][2]
                      or within_two(prediction[(c, d)], observed[(c, d)][0]))]
    r1 = not misses
    out.append(f"R1 level, caps 3 and 5: misses {', '.join(misses) or 'none'} -> "
               f"{'holds' if r1 else 'fails'} (rule: P inside the interval or within a factor of 2)")
    disagree = [f"c{c} {d}x" for c in CAPS for d in DEPTHS
                if (prediction[(c, d)] >= 0.5) != (observed[(c, d)][1] > 0)]
    r2 = not disagree
    out.append(f"R2 regime, {len(observed)} cells: disagreements {', '.join(disagree) or 'none'} -> "
               f"{'holds' if r2 else 'fails'} (rule: P >= 0.5 iff the interval's low end > 0)")
    decided = [k for k, (_, lo, hi) in observed.items() if lo > 0 or hi < 0]
    predictors = {"model": prediction, **baselines}
    mae = {n: sum(abs(p[k] - observed[k][0]) for k in observed) / len(observed)
           for n, p in predictors.items()}
    agree = {n: sum(sign(p[k]) == sign(observed[k][0]) for k in decided)
             for n, p in predictors.items()}
    for n in predictors:
        out.append(f"   {n:>5}: mean absolute error {mae[n]:.3f}; sign agreement "
                   f"{agree[n]}/{len(decided)}")
    others = [n for n in predictors if n != "model"]
    r3 = all(mae["model"] < mae[n] and agree["model"] >= agree[n] for n in others)
    out.append(f"R3 against B0, B1, B2 -> {'holds' if r3 else 'fails'} (rule: lowest mean "
               "absolute error, and sign agreement at least each baseline's)")
    verdict = "supported" if r1 and r2 else "partly supported" if r2 else "not supported"
    out.append(f"model {verdict} (R1 and R2); R3 {'holds' if r3 else 'fails'}")
    return out


def main() -> None:
    global DEPTHS
    args = sys.argv[1:]
    header = "O7 "
    if args[0] == "--hg005":
        # H5-3 of `hg005-run-plan.md`: rules R1 and R2 only, depths full and 10, no baselines.
        DEPTHS, header, args = ("full", "10"), "regime (descriptive)", args[1:]
    baselines_path, score_path, *prediction_paths = args
    sums = {}
    for path in prediction_paths:
        with open(path) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                acc = sums.setdefault((int(row["cap"]), row["depth"]), [0, 0.0])
                acc[0] += int(row["boundaries"])
                acc[1] += float(row["summed_difference"])
    prediction = {k: 1e4 * total / n for k, (n, total) in sums.items()}
    baselines = {b: {} for b in ("B0", "B1", "B2")}
    for line in open(baselines_path) if baselines_path != "-" else ():
        if line[0].isdigit():
            cap, depth, *values = line.split("\t")
            for b, v in zip(("B0", "B1", "B2"), values):
                baselines[b][(int(cap), depth)] = float(v)
    observed = rows(open(score_path).read(), header)
    if baselines_path == "-":
        baselines = {b: {k: 0.0 for k in observed} for b in ("B0", "B1", "B2")}
    print(f"registered prediction score (docs/design/ont-prediction.md): {baselines_path}, "
          f"{score_path}, {', '.join(prediction_paths)}; prediction pooled at full precision")
    lines = score(prediction, baselines, observed)
    if baselines_path == "-":
        lines = [line for line in lines if not line.startswith(("   ", "R3"))]
        lines[-1] = lines[-1].split(";")[0] + "; R3 not registered"
        lines = [line.replace("  B0     B1     B2", "") for line in lines]
        lines = [re.sub(r"(  [+-]\d+\.\d{4})(  [+-]\d+\.\d{2}){3}$", r"\1", line) for line in lines]
    print("\n".join(lines))


def demo() -> None:
    cells = [(c, d) for c in CAPS for d in DEPTHS]
    observed = {k: (1.0, 0.5, 1.5) if k[0] == 3 else (0.0, -0.1, 0.1) for k in cells}
    exact = {k: observed[k][0] for k in cells}
    zero = {k: 0.0 for k in cells}
    lines = score(exact, {"B0": zero, "B1": zero, "B2": zero}, observed)
    assert lines[-1] == "model supported (R1 and R2); R3 holds", lines[-1]
    lines = score(zero, {"B0": exact, "B1": exact, "B2": zero}, observed)
    assert lines[-1] == "model not supported (R1 and R2); R3 fails", lines[-1]


if __name__ == "__main__":
    if sys.argv[1:] == ["--demo"]:
        demo()
        print("demo ok")
    else:
        main()
