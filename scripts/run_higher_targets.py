"""Higher target levels with f = m - n = 1.

Under the authors' own criteria there are parameter sets with m = n + 1 for
targets above 128 bits.  There the y-only reduction leaves a single free
parameter, so the solve stage degenerates to the roots of one univariate
polynomial; the binding stage becomes the triangular-flag recovery, which
depends only on n.  This script measures each stage per set.

    python3 scripts/run_higher_targets.py --sets 160
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

from _common import peak_mib  # noqa: E402
from factodsa.attack import recover_quadratic_space  # noqa: E402
from factodsa.flag_recovery import recover_flag  # noqa: E402
from factodsa.inversion import bilinear_residuals  # noqa: E402
from factodsa.model import generate  # noqa: E402
from factodsa.polynomials import monomials as mono  # noqa: E402
from factodsa.recovery import extract_structure  # noqa: E402
from factodsa.structural_solve import (interpolate_residuals, reduce_bilinear,  # noqa: E402
                                       solve_reduced)

Q = 65519
SETS = {160: (14, 15), 192: (17, 18), 256: (22, 23)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sets", type=int, nargs="*", default=[160],
                        help="target levels to run, e.g. --sets 160 192; default "
                             "160, the only one driven end to end here")
    parser.add_argument("--flag-timeout", type=int, default=3600)
    parser.add_argument("--solve-timeout", type=int, default=1800)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    output = REPO / "results"
    output.mkdir(exist_ok=True)
    records = []
    for target in args.sets:
        if target not in SETS:
            print(f"unknown set {target}", file=sys.stderr)
            continue
        n, m = SETS[target]
        record = {"target_bits": target, "n": n, "m": m, "f": m - n}
        print(json.dumps(record), flush=True)

        P, _oracle = generate(n, m, Q, args.seed)
        started = time.perf_counter()
        arrays, info = extract_structure(P, Q)
        record["structure_seconds"] = round(time.perf_counter() - started, 2)
        record["structure_peak_mib"] = peak_mib()
        print(json.dumps(record), flush=True)

        quadratic_space = recover_quadratic_space(arrays["Qspace"] % Q, Q, n)
        started = time.perf_counter()
        L, Q_tri, ok = recover_flag(quadratic_space, Q, timeout=args.flag_timeout)
        record["flag_seconds"] = round(time.perf_counter() - started, 2)
        record["flag_ok"] = bool(ok)
        print(json.dumps(record), flush=True)
        if not ok:
            records.append(record)
            continue

        frame = arrays["frame"] % Q
        rng = np.random.default_rng(args.seed + 1)
        x0 = rng.integers(Q, size=n)
        y0 = rng.integers(Q, size=n)
        z0 = (frame @ np.concatenate([x0, y0])) % Q
        triples = mono(2 * n, 3)
        h = np.array([sum(int(P[i, k]) * int(z0[triples[k][0]]) % Q
                          * int(z0[triples[k][1]]) % Q * int(z0[triples[k][2]])
                          for k in range(len(triples))) % Q for i in range(m)])
        started = time.perf_counter()
        reduction = reduce_bilinear(arrays["G"] % Q, h, Q, n, m, seed=args.seed + 1)
        polynomials = interpolate_residuals(reduction, h, Q, n, m, seed=args.seed + 1)
        pairs, method = solve_reduced(reduction, polynomials, Q, n, m,
                                      timeout=args.solve_timeout, seed=args.seed + 1)
        record["solve_method"] = method
        record["free_parameters"] = reduction["free"]
        record["residual_degree_bound"] = polynomials["degree"]
        record["gr_seconds"] = round(time.perf_counter() - started, 2)
        G = arrays["G"] % Q
        record["pairs"] = len(pairs)
        record["pairs_valid"] = sum(1 for a, y in pairs
                                    if not np.any(bilinear_residuals(
                                        G, h, tuple(a) + tuple(y), Q, n)))
        record["peak_mib"] = peak_mib()
        records.append(record)
        print(json.dumps(record), flush=True)

    (output / "higher-targets.json").write_text(json.dumps(records, indent=2) + "\n")
    print(f"wrote {output / 'higher-targets.json'}")


if __name__ == "__main__":
    main()
