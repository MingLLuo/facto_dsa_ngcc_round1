"""Measure structural recovery on the four parameter sets the docs quote.

    python3 scripts/measure_structure_scaling.py

Sizes: the official 128, 256 and 512-bit sets, plus the 160-bit ``f = 1`` set.
Each size runs in its own child process, so the peak-memory reading belongs to
that size alone.  Writes ``results/structure-scaling.json``.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SETS = [(10, 13), (14, 15), (17, 32), (32, 62)]        # 128, 160 (f=1), 256, 512
CHILD = """
import json, sys, time
sys.path.insert(0, {src!r})
sys.path.insert(0, {scripts!r})
from _common import peak_mib
from factodsa.model import generate
from factodsa.recovery import extract_structure
n, m = {n}, {m}
P, _ = generate(n, m, 65519, 0)
started = time.perf_counter()
_, info = extract_structure(P, 65519)
print(json.dumps({{"n": n, "m": m, "seconds": round(time.perf_counter() - started, 2),
                   "peak_mib": peak_mib(), "status": info["status"]}}))
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=REPO / "results" / "structure-scaling.json")
    args = parser.parse_args()

    records = []
    for n, m in SETS:
        child = CHILD.format(src=str(REPO / "src"), scripts=str(REPO / "scripts"), n=n, m=m)
        completed = subprocess.run([sys.executable, "-c", child], check=True,
                                   capture_output=True, text=True)
        record = json.loads(completed.stdout.strip().splitlines()[-1])
        records.append(record)
        print(f"n={n:<3} m={m:<3} {record['seconds']:>6.2f} s  "
              f"{record['peak_mib']:>8.1f} MiB  {record['status']}", flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(records, indent=2) + "\n")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
