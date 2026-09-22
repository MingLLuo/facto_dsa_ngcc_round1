"""End-to-end forgery on small parameters (fast regression of the whole chain)."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _common import peak_mib  # noqa: E402
from factodsa.attack import run_case  # noqa: E402

SIZES = [(2, 3), (3, 5), (4, 7), (5, 9), (6, 11)]
Q = 65519


def main():
    output = REPO / "results"
    output.mkdir(exist_ok=True)
    records = []
    for n, m in SIZES:
        forged = 0
        times = []
        for seed in range(3):
            started = time.perf_counter()
            report = run_case(n, m, Q, seed=seed, target_seed=seed + 1, solve_seed=seed + 2,
                              timeout=300, flag_timeout=120, slices=6)
            times.append(time.perf_counter() - started)
            forged += report["status"] == "forged"
        row = {"n": n, "m": m, "f": m - n, "forged": f"{forged}/3",
               "median_seconds": round(sorted(times)[1], 2)}
        records.append(row)
        print(json.dumps(row), flush=True)
    summary = {"q": Q, "sizes": records, "peak_mib": peak_mib(),
               "route": "structural (y-only slice) + triangular inversion",
               "verification": "every public cubic coefficient of P(z) compared with h"}
    (output / "small-validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"forged {sum(int(r['forged'].split('/')[0]) for r in records)}/{3 * len(SIZES)}")


if __name__ == "__main__":
    main()
