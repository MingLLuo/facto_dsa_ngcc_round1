"""Standard 128-bit instance (q=65519, n=10, m=13) end to end, with timings."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _common import peak_mib  # noqa: E402
from factodsa.attack import run_case  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--slices", type=int, default=6)
    args = parser.parse_args()
    started = time.perf_counter()
    report = run_case(10, 13, 65519, seed=args.seed, target_seed=1, solve_seed=1,
                      timeout=args.timeout, flag_timeout=args.timeout,
                      slices=args.slices, rescalings=48)
    report["wall_seconds"] = round(time.perf_counter() - started, 2)
    report["peak_mib"] = peak_mib()
    output = REPO / "results"
    output.mkdir(exist_ok=True)
    (output / f"standard-128-seed{args.seed}.json").write_text(
        json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "z"}, indent=2))
    print("forged" if report["status"] == "forged" else report["status"])


if __name__ == "__main__":
    main()
