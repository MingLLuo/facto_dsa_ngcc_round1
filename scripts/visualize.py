"""Local figures: pipeline, stage timings, and the complexity picture.

Everything is written to ``figures/`` as PNG; matplotlib is the only extra
dependency and the script degrades to a text summary if it is missing.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FIGURES = REPO / "figures"

# keep matplotlib from writing into an unwritable home directory
os.environ.setdefault("MPLCONFIGDIR", str(REPO / ".mplconfig"))

LEVELS = {128: (10, 13), 256: (17, 32), 512: (32, 62)}
F1_SETS = [(14, 15, 160), (17, 18, 192), (22, 23, 256)]
PRODUCERS = {"standard-128-seed0.json": "scripts/run_standard_128.py",
             "structure-scaling.json": "scripts/measure_structure_scaling.py"}


def recorded(name):
    """Load one recorded result file, with a clear message when it is missing."""
    path = REPO / "results" / name
    if not path.exists():
        raise SystemExit(f"{path} is missing; run {PRODUCERS.get(name, 'the matching script')} first")
    return json.loads(path.read_text())


def bezout_reduced(n, m):
    return n ** (m - n)


def main():
    FIGURES.mkdir(exist_ok=True)
    standard = recorded("standard-128-seed0.json")
    timings = {"recover": standard["recover_seconds"], "flag": standard["flag_seconds"],
               "solve": standard["solve_seconds"]}
    scaling = recorded("structure-scaling.json")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
    except Exception as exc:
        print(f"matplotlib unavailable ({exc}); writing a text summary instead")
        summary = {
            "stage_timings_standard_128": timings,
            "structure_scaling": scaling,
            "levels": {str(k): {"n": v[0], "m": v[1], "f": v[1] - v[0],
                                "bezout_reduced": bezout_reduced(*v)}
                       for k, v in LEVELS.items()},
        }
        (REPO / "results" / "visualization-summary.json").write_text(
            json.dumps(summary, indent=2) + "\n")
        return

    # 1. pipeline diagram
    fig, ax = plt.subplots(figsize=(11, 2.6))
    ax.set_axis_off()
    steps = ["public key\n(e coefficients)", "structural\nrecovery",
             "triangular flag\n(rank-1 peeling)", "reduced bilinear system\n"
             "(f = m - n unknowns)", "reduced solve\n(msolve parametrisation)",
             "image condition\n+ rescaling", "triangular\ninversion", "verify\nP(z) = h"]
    x = 0.0
    width = 0.100
    gap = 0.028
    for index, label in enumerate(steps):
        ax.add_patch(FancyBboxPatch((x, 0.36), width, 0.32, boxstyle="round,pad=0.006",
                                    linewidth=1.0, edgecolor="#333", facecolor="#eef3fb",
                                    zorder=2))
        ax.text(x + width / 2, 0.52, label, ha="center", va="center", fontsize=7,
                zorder=3)
        if index + 1 < len(steps):
            ax.add_patch(FancyArrowPatch((x + width + 0.002, 0.52),
                                         (x + width + gap - 0.002, 0.52),
                                         arrowstyle="-|>", mutation_scale=10,
                                         linewidth=1.1, color="#333", zorder=4))
        x += width + gap
    ax.set_xlim(-0.01, x + 0.01)
    ax.set_ylim(0.22, 0.82)
    ax.set_title("Facto-DSA attack chain (public key only)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGURES / "pipeline.png", dpi=200)
    plt.close(fig)

    # 2. stage timings, standard 128
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    names = list(timings)
    values = [timings[k] for k in names]
    bars = ax.bar(names, values, color=["#8fbcd4", "#f0a868", "#9fd6a0"])
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.4, f"{value:.2f}s",
                ha="center", fontsize=8)
    ax.set_ylabel("seconds")
    ax.set_title("Standard 128-bit instance (n=10, m=13)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGURES / "stage-timings-128.png", dpi=200)
    plt.close(fig)

    # 3. structure recovery scaling: time and memory
    ns = [row["n"] for row in scaling]
    secs = [row["seconds"] for row in scaling]
    mems = [row["peak_mib"] for row in scaling]
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    ax.loglog(ns, secs, "o-", label="seconds")
    ax.loglog(ns, mems, "s--", label="peak MiB")
    ax.set_xlabel("n")
    ax.set_title("Structural recovery cost (measured)", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURES / "structure-scaling.png", dpi=200)
    plt.close(fig)

    # 4. Bezout of the reduced system vs f = m - n
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    fs = list(range(1, 16))
    for n, color in [(10, "#1f77b4"), (17, "#d62728"), (32, "#2ca02c")]:
        values = [math.log10(max(bezout_reduced(n, n + f), 1)) for f in fs]
        ax.plot(fs, values, "-", color=color, label=f"n={n}")
    ax.axhline(math.log10(1e9), color="gray", linestyle=":", linewidth=1)
    ax.text(9.6, math.log10(1e9) + 0.3, "10^9 (practical)", fontsize=7, color="gray")
    offsets = [(8, 6), (8, -8), (8, -22)]
    for (n, m, level), offset in zip(F1_SETS, offsets):
        value = math.log10(bezout_reduced(n, m))
        ax.plot([m - n], [value], "k*", markersize=9, zorder=5)
        ax.annotate(f"f=1 (n={n}), {level}-bit", (m - n, value),
                    textcoords="offset points", xytext=offset, fontsize=7,
                    arrowprops=dict(arrowstyle="-", linewidth=0.6, color="#444"))
    ax.set_xlabel("f = m - n  (free parameters after the y-only reduction)")
    ax.set_ylabel("log10 Bezout number of the reduced system")
    ax.set_title("Why f = m - n decides the attack cost", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURES / "bezout-vs-f.png", dpi=200)
    plt.close(fig)

    # 5. why it breaks, in two pictures
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.4, 4.0),
                                      gridspec_kw={"width_ratios": [1, 1.35]})
    ambient, derivative = 210, 165          # n = 10: n(2n+1) and 3n(n+1)/2
    defect = ambient - derivative           # = n(n-1)/2 = 45
    left.bar(["all quadratic\nforms", "the derivative\nspace"], [ambient, derivative],
             color=["#c9d8e8", "#8fbcd4"])
    left.set_ylabel("dimension")
    left.set_title("The public derivative space is missing\nexactly n(n-1)/2 = 45 dimensions",
                   fontsize=9.5)
    left.annotate("", xy=(1, derivative), xytext=(1, ambient),
                  arrowprops=dict(arrowstyle="<->", color="#c0392b", linewidth=1.4))
    left.text(1.12, (ambient + derivative) / 2, "45\nthe X-squared block\nthat names the\nhidden splitting",
              color="#c0392b", fontsize=7.6, va="center")
    left.set_ylim(0, 235)
    left.grid(axis="y", alpha=0.3)

    right.set_axis_off()
    right.set_title("and then the target equation shrinks", fontsize=9.5)
    rows = [
        ("public target: 13 cubic equations in 20 unknowns", "#c9d8e8"),
        ("always 7 free dimensions  (20 - 13)", "#dbe7f3"),
        ("slice only the y-block, eliminate the a-block", "#dbe7f3"),
        ("3 equations of degree <= 10 in 3 unknowns", "#9fd6a0"),
    ]
    for index, (text, colour) in enumerate(rows):
        y = 0.78 - index * 0.21
        right.add_patch(FancyBboxPatch((0.03, y - 0.075), 0.94, 0.15,
                                       boxstyle="round,pad=0.008",
                                       linewidth=1.0, edgecolor="#333", facecolor=colour))
        right.text(0.5, y, text, ha="center", va="center", fontsize=8.4)
    right.set_xlim(0, 1)
    right.set_ylim(0, 1)
    fig.suptitle("Why the construction leaks (standard 128-bit set: n = 10, m = 13)",
                 fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIGURES / "why-it-breaks.png", dpi=200)
    plt.close(fig)

    print(f"wrote figures to {FIGURES}")
    for path in sorted(FIGURES.glob("*.png")):
        print("  ", path.name)


if __name__ == "__main__":
    main()
