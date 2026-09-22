"""Independent re-verification of the standard 128-bit result, in SageMath.

Run:  sage -python scripts/verify_with_sage.py      (from the repository root)

This is a second toolchain, not a second implementation of the attack.  It
imports nothing from ``factodsa`` and uses no secret key; it loads the authors'
KAT public key and the forged signature recorded in ``results/acceptance/`` and
re-checks the three claims of the finding with Sage objects:

  1. the forged ``z`` satisfies ``P(z) = h`` (multivariate polynomial
     substitution, exactly, in GF(q));
  2. the public first-derivative space has rank ``3n(n+1)/2 = 165`` out of 210,
     so it is ``n(n-1)/2 = 45`` dimensions short (matrices over GF(q));
  3. the images of the trace annihilator of that space span exactly the
     ``n = 10`` hidden X-directions, and every public cubic vanishes on them.

The heavy solving stage stays with msolve: Sage's ``ideal.variety()`` needs a
lex Groebner basis and ran for minutes on the reduced system that msolve
parametrises in under a second.
"""
import json
import struct
import sys
from itertools import combinations_with_replacement
from pathlib import Path

try:
    from sage.all import GF, PolynomialRing, matrix, vector
except ImportError:                                  # pragma: no cover
    sys.exit("SageMath not available: run this with `sage -python scripts/"
             "verify_with_sage.py` (or `make sage-check`)")

ROOT = Path(__file__).resolve().parents[1]
ACCEPTANCE = ROOT / "results" / "acceptance"

Q, N, M = 65519, 20, 13
n = N // 2
CUBIC, QUADRATIC = 1540, 210                      # C(22,3) and C(21,2)
D_STAR = 3 * n * (n + 1) // 2                     # 165

F = GF(Q)
triples = list(combinations_with_replacement(range(N), 3))
pairs = list(combinations_with_replacement(range(N), 2))
pair_index = {pair: i for i, pair in enumerate(pairs)}
failed = []


def check(label, ok, detail=""):
    print(f"[{'ok  ' if ok else 'FAIL'}] {label}{'  ' + detail if detail else ''}")
    if not ok:
        failed.append(label)


if not (ACCEPTANCE / "signature.bin").exists():
    sys.exit("run `python3 scripts/run_acceptance.py --source kat` first")

record = json.loads((ACCEPTANCE / "acceptance.json").read_text())
public_key = (ACCEPTANCE / "pk.bin").read_bytes()
signature = (ACCEPTANCE / "signature.bin").read_bytes()
check("public key size", len(public_key) == 2 * M * CUBIC, f"{len(public_key)} bytes")
check("signature size", len(signature) == 2 * N, f"{len(signature)} bytes")

P = matrix(F, M, CUBIC, struct.unpack(f"<{M * CUBIC}H", public_key))
z = vector(F, struct.unpack(f"<{N}H", signature))
h = vector(F, record["target"])

# 1. the forged signature meets the recorded protocol target ------------------
R = PolynomialRing(F, N, "z")
Z = R.gens()
images = []
for i in range(M):
    images.append(sum((P[i, k] * Z[t[0]] * Z[t[1]] * Z[t[2]]
                       for k, t in enumerate(triples)), R(0)))
check("forged signature satisfies P(z) = h",
      vector(F, [poly(*z) for poly in images]) == h)
check("recorded sig_verify return code", record["sig_verify_return_code"] == "0")
check("flipped-bit control", record["negative_control_flipped_bit"] == "-1")

# 2. the derivative space is 45 dimensions short of the ambient 210 -----------
derivative = matrix(F, M * N, QUADRATIC)
for i in range(M):
    for k, t in enumerate(triples):
        if not P[i, k]:
            continue
        for v in set(t):
            rest = list(t)
            rest.remove(v)
            derivative[i * N + v, pair_index[tuple(rest)]] += P[i, k] * t.count(v)
rank = derivative.rank()
check("derivative rank", rank == D_STAR, f"{rank} of {QUADRATIC}")
check("missing dimensions", QUADRATIC - rank == n * (n - 1) // 2,
      f"{QUADRATIC - rank} = n(n-1)/2")

# 3. the annihilator images span exactly the hidden X-half -------------------
kernel = derivative.right_kernel_matrix()
annihilator = [matrix(F, N, N) for _ in range(kernel.nrows())]
for r in range(kernel.nrows()):
    for a, b in pairs:
        annihilator[r][a, b] = annihilator[r][b, a] = kernel[r, pair_index[(a, b)]]
support = matrix(F, [list(A.column(b)) for A in annihilator for b in range(N)]).T
check("annihilator image dimension", support.rank() == n, f"{support.rank()} = n")

W = support.column_space().basis_matrix().T                   # N x n
SU = PolynomialRing(F, n, "u")
U = SU.gens()
image_of_z = []
for i in range(N):
    image_of_z.append(sum((W[i, j] * U[j] for j in range(n)), SU(0)))
restrict_to_W = R.hom(image_of_z, SU)                         # z = W u
check("every public cubic vanishes on W",
      all(restrict_to_W(poly) == SU(0) for poly in images), "P|_W = 0")

print()
if failed:
    print("SAGE CHECK FAILED: " + ", ".join(failed))
else:
    print("SAGE CHECK PASSED: P(z) = h, d* = 165 of 210, "
          "and W is the n = 10 hidden X-half")
sys.exit(1 if failed else 0)
