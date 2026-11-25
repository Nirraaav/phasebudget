"""Figures and a summary table from `benchmark_phase_diagram.py`'s CSV.

1. `phase_diagram.png`: small multiples, base-quality spread by rows and
   variant density by columns; each panel is depth against error rate, shaded
   by the fraction of instances whose truth is a wMEC optimum, every cell
   labelled, with a dark line where that fraction crosses one half.
2. `oq5_informativeness.png`: every cell's informativeness, for the MEC margin
   and the posterior margin, against its mean objective-truth gap, with the
   median per gap bin drawn as a line; 0.5 is no information.

It also prints, per panel, the largest error rate at which the truth is an
optimum in at least half the instances, which is the boundary in numbers.
matplotlib is needed here and nowhere in the library.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

# Reference palette of the dataviz method: sequential blue, light to dark, and
# the first two categorical slots, validated for colour-vision deficiency.
RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
BLUE, ORANGE = "#2a78d6", "#eb6834"
SURFACE, INK, SECONDARY, MUTED, GRID, AXIS = (
    "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
)


def read_rows(path: Path):
    comments, lines = [], []
    with open(path) as handle:
        for line in handle:
            (comments if line.startswith("#") else lines).append(line)
    rows = list(csv.DictReader(lines))
    for row in rows:
        for key, value in row.items():
            row[key] = float(value)
    return [c.lstrip("# ").rstrip() for c in comments], rows


def style(ax):
    ax.set_facecolor(SURFACE)
    for spine in ax.spines.values():
        spine.set_color(AXIS)
        spine.set_linewidth(0.8)
    ax.tick_params(colors=SECONDARY, labelsize=8, length=0)


def staircase(ax, grid):
    """The one-half boundary drawn along cell edges, so it never crosses a label.

    Cells are on a categorical grid, so the boundary is where a row changes
    from at least one half to below it; vertical segments separate the two
    sides of each row and horizontal segments join rows whose boundary moves.
    """
    edges = []
    for y, row in enumerate(grid):
        passing = [x for x, value in enumerate(row) if not math.isnan(value) and value >= 0.5]
        edges.append(max(passing) + 0.5 if passing else -0.5)
    for y, x in enumerate(edges):
        ax.plot([x, x], [y - 0.5, y + 0.5], color=INK, linewidth=2, solid_capstyle="butt")
        if y + 1 < len(edges) and edges[y + 1] != x:
            ax.plot([x, edges[y + 1]], [y + 0.5, y + 0.5], color=INK, linewidth=2,
                    solid_capstyle="butt")


def phase_diagram(rows, out: Path, provenance):
    densities = sorted({r["variants_per_read"] for r in rows})
    spreads = sorted({r["quality_spread"] for r in rows})
    errors = sorted({r["error_rate"] for r in rows})
    depths = sorted({r["coverage"] for r in rows})
    cmap = LinearSegmentedColormap.from_list("seq_blue", RAMP)
    fig, axes = plt.subplots(
        len(spreads), len(densities),
        figsize=(3.4 * len(densities) + 1.2, 3.0 * len(spreads) + 1.0),
        squeeze=False, facecolor=SURFACE,
    )
    image = None
    for i, spread in enumerate(spreads):
        for k, density in enumerate(densities):
            ax = axes[i][k]
            style(ax)
            grid = [[math.nan] * len(errors) for _ in depths]
            for r in rows:
                if r["quality_spread"] == spread and r["variants_per_read"] == density:
                    grid[depths.index(r["coverage"])][errors.index(r["error_rate"])] = r["truth_optimal"]
            image = ax.imshow(grid, origin="lower", cmap=cmap, vmin=0, vmax=1, aspect="auto")
            staircase(ax, grid)
            for y, row in enumerate(grid):
                for x, value in enumerate(row):
                    if not math.isnan(value):
                        ax.text(x, y, f"{value:.2f}", ha="center", va="center", fontsize=7,
                                color="#ffffff" if value > 0.55 else INK)
            ax.set_xticks(range(len(errors)), [f"{e:g}" for e in errors])
            ax.set_yticks(range(len(depths)), [f"{d:g}" for d in depths])
            ax.set_title(f"{density:g} variants per read, quality spread {spread:g}",
                         fontsize=9, color=INK)
            if i == len(spreads) - 1:
                ax.set_xlabel("per-base error rate", fontsize=8, color=SECONDARY)
            if k == 0:
                ax.set_ylabel("sequencing depth", fontsize=8, color=SECONDARY)
    bar = fig.colorbar(image, ax=axes, shrink=0.8, pad=0.02)
    bar.set_label("fraction of instances whose truth is a wMEC optimum", fontsize=8, color=SECONDARY)
    bar.ax.tick_params(colors=SECONDARY, labelsize=7, length=0)
    bar.outline.set_edgecolor(AXIS)
    fig.suptitle("Where the wMEC optimum is the truth (dark line: edge of one half)",
                 fontsize=11, color=INK)
    fig.text(0.01, 0.005, "; ".join(provenance[:2]), fontsize=5.5, color=MUTED)
    fig.savefig(out, dpi=180, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)


def binned_median(points, edges):
    line = []
    for low, high in zip(edges, edges[1:]):
        inside = [y for x, y in points if low <= x < high]
        if len(inside) >= 3:
            line.append(((low * high) ** 0.5 if low > 0 else high / 2, statistics.median(inside)))
    return line


def informativeness(rows, out: Path, provenance, minimum: int):
    usable = [r for r in rows if r["instances_with_switch"] >= minimum]
    fig, ax = plt.subplots(figsize=(7.2, 4.4), facecolor=SURFACE)
    style(ax)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    gaps = [max(r["gap"], 1e-3) for r in usable]
    edges = [0.0] + [10 ** (e / 4) for e in range(-8, 10)]
    for key, color, name in (("margin_info", BLUE, "MEC margin"),
                             ("posterior_info", ORANGE, "posterior margin")):
        points = [(g, r[key]) for g, r in zip(gaps, usable) if not math.isnan(r[key])]
        ax.scatter([p[0] for p in points], [p[1] for p in points], s=16, color=color,
                   alpha=0.45, edgecolors=SURFACE, linewidths=0.8, label=name)
        line = binned_median(points, edges)
        if line:
            ax.plot([p[0] for p in line], [p[1] for p in line], color=color, linewidth=2)
            ax.annotate(name, line[-1], xytext=(6, 0), textcoords="offset points",
                        fontsize=8, color=SECONDARY, va="center")
    ax.axhline(0.5, color=MUTED, linewidth=1, linestyle=(0, (4, 3)))
    ax.set_xscale("log")
    ticks = [t for t in (0.1, 0.3, 1, 3, 10, 30) if min(gaps) <= t <= max(gaps) * 1.5]
    ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.text(max(gaps), 0.505, "no information", fontsize=7, color=MUTED, va="bottom", ha="right")
    ax.set_xlabel("mean objective-truth gap of the cell (wMEC cost of truth minus optimum)",
                  fontsize=8, color=SECONDARY)
    ax.set_ylabel("P(signal at a switch error < signal elsewhere)", fontsize=8, color=SECONDARY)
    ax.set_title(f"Signal informativeness against the objective-truth gap "
                 f"({len(usable)} cells with at least {minimum} instances carrying a switch)",
                 fontsize=10, color=INK)
    legend = ax.legend(frameon=False, fontsize=8, loc="upper right")
    for text in legend.get_texts():
        text.set_color(SECONDARY)
    fig.text(0.01, 0.005, "; ".join(provenance[:2]), fontsize=5.5, color=MUTED)
    fig.savefig(out, dpi=180, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    return usable


def boundary_table(rows):
    print("largest error rate with truth optimal in at least half the instances")
    print(f"{'variants/read':>13} {'spread':>7} " + " ".join(
        f"{'depth ' + format(d, 'g'):>9}" for d in sorted({r['coverage'] for r in rows})))
    for density in sorted({r["variants_per_read"] for r in rows}):
        for spread in sorted({r["quality_spread"] for r in rows}):
            cells = []
            for depth in sorted({r["coverage"] for r in rows}):
                passing = [r["error_rate"] for r in rows
                           if r["variants_per_read"] == density and r["quality_spread"] == spread
                           and r["coverage"] == depth and r["truth_optimal"] >= 0.5]
                cells.append(f"{max(passing):>9g}" if passing else f"{'none':>9}")
            print(f"{density:>13g} {spread:>7g} " + " ".join(cells))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=Path("benchmarks/results/phase_diagram.csv"))
    parser.add_argument("--out-dir", type=Path, default=Path("benchmarks/results"))
    parser.add_argument("--minimum", type=int, default=5)
    args = parser.parse_args()
    provenance, rows = read_rows(args.csv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    phase_diagram(rows, args.out_dir / "phase_diagram.png", provenance)
    usable = informativeness(rows, args.out_dir / "oq5_informativeness.png", provenance, args.minimum)
    for line in provenance:
        print(line)
    boundary_table(rows)
    hard = [r["margin_info"] for r in usable]
    soft = [r["posterior_info"] for r in usable]
    better = sum(s > h for h, s in zip(hard, soft))
    print(f"cells with a usable informativeness: {len(usable)}; posterior higher in {better}")
    correlation = statistics.correlation([math.log(max(r["gap"], 1e-3)) for r in usable], hard)
    print(f"correlation of log gap with MEC-margin informativeness: {correlation:+.3f}")
    correlation = statistics.correlation([math.log(max(r["gap"], 1e-3)) for r in usable], soft)
    print(f"correlation of log gap with posterior informativeness: {correlation:+.3f}")


if __name__ == "__main__":
    main()
