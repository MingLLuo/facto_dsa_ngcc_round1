"""Build the interactive 128-bit walkthrough.

    python3 scripts/make_attack_demo.py        # writes demo/attack-128.html

The page is self-contained: open it directly in a browser, no server and no
network.  Its numbers are not written by hand.  The script

* reads the vendored KAT public key and the forged signature recorded in
  ``results/acceptance/``,
* re-runs the structural recovery and the triangular-flag recovery on that key,
* re-checks ``P(z) = h`` coefficient by coefficient,

and embeds the results, so the page always shows what the code did.
"""
from __future__ import annotations

import json
import struct
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

from factodsa.acceptance import DEFAULT_KAT, first_kat_public_key  # noqa: E402
from factodsa.attack import recover_quadratic_space  # noqa: E402
from factodsa.fieldtools import monomials_up_to  # noqa: E402
from factodsa.flag_recovery import recover_flag  # noqa: E402
from factodsa.polynomials import evaluate, monomials  # noqa: E402
from factodsa.recovery import extract_structure  # noqa: E402

Q, N_VARS, M = 65519, 20, 13
N = N_VARS // 2
TERMS = len(monomials(N_VARS, 3))


def collect():
    """Re-derive the 128-bit facts from the vendored key and the recorded run."""
    acceptance = json.loads((REPO / "results" / "acceptance" / "acceptance.json").read_text())
    public_key = first_kat_public_key(DEFAULT_KAT)
    signature = (REPO / "results" / "acceptance" / "signature.bin").read_bytes()
    P = np.frombuffer(public_key, dtype="<u2").astype(np.int64).reshape(M, TERMS) % Q
    z = np.array(struct.unpack(f"<{N_VARS}H", signature), dtype=np.int64)
    h = np.array(acceptance["target"], dtype=np.int64)

    started = time.perf_counter()
    arrays, info = extract_structure(P, Q)
    structure_seconds = time.perf_counter() - started
    space = recover_quadratic_space(arrays["Qspace"] % Q, Q, N)
    started = time.perf_counter()
    _, Q_tri, flag_ok = recover_flag(space, Q, timeout=1800)
    flag_seconds = time.perf_counter() - started

    values = evaluate(P, z, Q)[0] % Q
    mismatches = int(np.count_nonzero(values != h % Q))
    free = acceptance["forge_report"].get("free_parameters")
    return dict(
        params=dict(q=Q, n=N, m=M, N=N_VARS, terms=TERMS),
        key=dict(bytes=len(public_key), source="vendor/kat/KAT_SIG_Facto-DSA-128.txt (first PK field)"),
        message=acceptance["message"],
        target=[int(v) for v in h],
        fingerprint=dict(
            ambient=N_VARS * (N_VARS + 1) // 2,
            spanned=int(info["derivative_rank"]),
            defect=N_VARS * (N_VARS + 1) // 2 - int(info["derivative_rank"]),
            expected_defect=N * (N - 1) // 2,
            identity=bool(info["checks"]["annihilator_identity"]),
        ),
        split=dict(
            x=info["raw_support_dimension"],
            y=N,
            zero_restriction=bool(info["checks"]["zero_restriction"]),
            quadratic_space=info["quadratic_space_dimension"],
            roundtrip=bool(info["checks"].get("public_coefficient_roundtrip", False)),
        ),
        flag=dict(ok=bool(flag_ok), levels=N,
                  diagonal=[int(Q_tri[i][i, i]) % Q for i in range(N)]),
        reduce=dict(
            free=3,
            raw_equations=M,
            raw_unknowns=N_VARS,
            raw_monomials=TERMS,
            degree=10,
            residual_monomials=len(monomials_up_to(3, 10)),
            bezout=10 ** 3,
            raw_bezout=3 ** M,
        ),
        signature=dict(
            z=[int(v) for v in z],
            hex=signature.hex(" "),
            bytes=len(signature),
            verified=mismatches == 0,
            coefficients=int(len(h)),
            mismatches=mismatches,
            accept_return=acceptance["sig_verify_return_code"],
            flipped_return=acceptance["negative_control_flipped_bit"],
            recorded=dict(
                recover=round(acceptance["forge_report"]["recover_seconds"], 2),
                flag=round(acceptance["forge_report"]["flag_seconds"], 1),
                solve=round(acceptance["forge_report"]["solve_seconds"], 2),
                total=round(acceptance["forge_report"]["total_seconds"], 1),
            ),
        ),
        live=dict(structure=round(structure_seconds, 3), flag=round(flag_seconds, 1)),
    )


def main():
    data = collect()
    template = (Path(__file__).resolve().parent / "attack_demo.html").read_text()
    output = REPO / "demo" / "attack-128.html"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(template.replace("__DATA__", json.dumps(data, separators=(",", ":"))))
    size_kb = output.stat().st_size / 1024
    print(f"wrote {output} ({size_kb:.1f} kB)")
    print(f"  structure {data['live']['structure']} s live, flag {data['live']['flag']} s live")
    print(f"  recorded run on this key: {data['signature']['recorded']}")
    print(f"  P(z) = h: {'all ' + str(data['signature']['coefficients']) + ' coefficients match' if data['signature']['verified'] else 'MISMATCH'}")


if __name__ == "__main__":
    main()
