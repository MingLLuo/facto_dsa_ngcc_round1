"""Public-coefficient-only extraction with exact algebraic certificates.

Nothing in this module accepts or imports a model, an oracle or a secret key.
The input is the public coefficient matrix; the output is a separated
polynomial representation of the same map, together with the certificates that
prove the representation is exact - every public coefficient is rebuilt at the
end.

The stages are the ones described in ``docs/attack.md``:

1. the first-derivative space and its trace annihilator, which expose the hidden
   splitting of the variables;
2. one linear system per coordinate, which pins the complement of that split;
3. the quadratic factor space, the multiplication tensor and the pure-``y``
   cubic, read off the bidegree split;
4. a right symmetriser of the tensor and the quadratic modifier of the map.
"""
from __future__ import annotations

import time
from math import isqrt

from .finite_field import (column_basis, inverse, kernel, matmul, np, rank, rref,
                           solve)
from .modifier import recover_modifier
from .polynomials import (annihilator_matrices, derivatives, join_bidegrees,
                          monomials, quadratic_matrices, split_bidegrees,
                          substitute_cubics, substitute_quadratics,
                          variables_times_quadratics)
from .symmetrizer import recover_symmetrizer


class StructureInconclusive(RuntimeError):
    """The public coefficients do not satisfy the structural conditions."""

    def __init__(self, message, details=None):
        super().__init__(message)
        self.details = details or {}


def _separate(C21, C12, C03, n, q, max_trials):
    """Recover the complement of ``ker L2`` from the ``x y^2`` coefficients.

    ``Hess[i, j]`` is the ``x``-Hessian of the coefficient of ``y_j``.  For a
    probe direction ``y`` the matrix ``sum_j Hess[i, j] y_j`` has full column
    rank, and its equation against ``C12`` is solved exactly; ``n`` independent
    probes glued together give the map ``H`` whose substitution removes the
    whole ``x y^2`` block.  The final identity check is the certificate.
    """
    Hess = 2 * quadratic_matrices(C21.transpose(0, 2, 1), n, q) % q
    rng = np.random.default_rng(20260921)
    directions = []
    solutions = []
    trials = 0
    for it in range(max_trials):
        y = (np.eye(n, dtype=np.int64)[it] if it < n
             else rng.integers(q, size=n, dtype=np.int64))
        trials += 1
        if rank(np.asarray(directions + [y], dtype=np.int64), q) <= len(directions):
            continue
        M = np.einsum("ijab,j->iab", Hess, y, optimize=True) % q
        M = M.reshape(-1, n)
        if rank(M, q) != n:
            continue
        quadratic_values = np.array([int(y[a]) * int(y[b]) % q for a, b in monomials(n, 2)],
                                    dtype=np.int64)
        rhs = -(C12.reshape(-1, len(quadratic_values)) @ quadratic_values) % q
        solution = solve(M, rhs, q)
        directions.append(y)
        solutions.append(solution)
        if len(directions) == n:
            break
    if len(directions) < n:
        raise StructureInconclusive("Insufficient independent full-rank complement probes",
                                    dict(probes=trials, accepted=len(directions)))

    # Glue the probe solutions into the complement map H.
    H = matmul(np.array(solutions).T, inverse(np.array(directions).T, q), q)
    # Exact x y^2 coefficients of D_x C21[H y] must cancel C12 entirely.
    polar = np.einsum("ijab,bk->iajk", Hess, H, optimize=True) % q
    correction = np.stack([polar[:, :, a, b] if a == b
                           else (polar[:, :, a, b] + polar[:, :, b, a]) % q
                           for a, b in monomials(n, 2)], axis=-1)
    if np.any((C12 + correction) % q):
        raise StructureInconclusive(
            "Complement probe solution fails the full coefficient identity")
    c21_at_H = substitute_quadratics(C21.transpose(0, 2, 1), H, q)
    c12_at_H = np.einsum("iap,ak->ikp", C12, H, optimize=True) % q
    new_C03 = (C03 + variables_times_quadratics(c21_at_H, q)
               + variables_times_quadratics(c12_at_H, q)) % q
    return H, new_C03, dict(probes=trials, independent_probes=n,
                            each_probe_full_column_rank=True,
                            xy2_coefficient_identity=True)


def extract_structure(P, q, *, max_probe_trials=512, roundtrip=True):
    """Recover the separated form of the public map from its coefficients only.

    Returns ``(arrays, info)``.  ``arrays`` holds the recovered subspaces and
    tensors; ``info`` records the certificate results and the per-stage timings.
    Raises :class:`StructureInconclusive` when a structural condition fails, and
    :class:`ValueError` when the input itself is malformed.
    """
    started = time.perf_counter()
    stages = {}
    checks = {}

    if isinstance(q, (bool, np.bool_)) or not isinstance(q, (int, np.integer)):
        raise TypeError("q must be an integer prime")
    q = int(q)
    if not 3 < q <= 65519 or any(q % d == 0 for d in range(2, isqrt(q) + 1)):
        raise ValueError("Expected a prime 3 < q <= 65519")
    raw = np.asarray(P)
    if raw.dtype.kind not in "iu":
        raise TypeError("Public coefficients must be integers; floating inputs are not rounded")
    P = np.ascontiguousarray(raw % q, dtype=np.int64)
    if P.ndim != 2:
        raise ValueError("Expected output rows and canonical cubic columns")
    N = next((N for N in range(4, 129, 2) if N * (N + 1) * (N + 2) // 6 == P.shape[1]), None)
    if N is None:
        raise ValueError("Cubic column count does not match an even input dimension")
    if q in (2, 3):
        raise ValueError("This implementation requires characteristic > 3")
    n = N // 2
    m = len(P)
    dstar = 3 * n * (n + 1) // 2
    if not n < m <= 2 * n - 1:
        raise ValueError("Expected n < m <= 2n-1")

    # 1. Derivative space and its trace annihilator.
    clock = time.perf_counter()
    D = derivatives(P, N, q)
    K = kernel(D, q)
    d_rank = D.shape[1] - len(K)
    stages["derivative_and_kernel"] = time.perf_counter() - clock
    if d_rank > dstar:
        raise StructureInconclusive("Derivative rank exceeds the model bound",
                                    dict(derivative_rank=d_rank, structural_bound=dstar))
    clock = time.perf_counter()
    A = annihilator_matrices(K, N, q)
    support = A.transpose(1, 0, 2).reshape(N, -1)
    W = column_basis(support, q)
    details = dict(derivative_shape=list(D.shape), derivative_rank=d_rank,
                   structural_bound=dstar, annihilator_dimension=len(K),
                   raw_support_dimension=W.shape[1])
    if d_rank != dstar or W.shape[1] != n:
        raise StructureInconclusive(
            "Derivative support does not identify an n-dimensional space", details)
    checks["annihilator_identity"] = not np.any(matmul(D, K.T, q))
    stages["support"] = time.perf_counter() - clock
    del D, K, A, support

    # 2. Complete W to a basis by coordinate vectors outside its pivot rows, and
    #    substitute it into the public cubic; the pure x^3 block must vanish.
    _, pivots = rref(W.T, q)
    missing = [i for i in range(N) if i not in pivots]
    C0 = np.eye(N, dtype=np.int64)[:, missing]
    B0 = np.concatenate([W, C0], axis=1)
    clock = time.perf_counter()
    rewritten = substitute_cubics(P, B0, q)
    C30, C21, C12, C03 = split_bidegrees(rewritten, n)
    checks["zero_restriction"] = not np.any(C30)
    if not checks["zero_restriction"]:
        raise StructureInconclusive("Candidate support fails P restricted to W = 0",
                                    details)
    stages["public_cubic_substitution"] = time.perf_counter() - clock

    # 3. One linear system per coordinate removes the x y^2 block.
    clock = time.perf_counter()
    H, C03, probe_checks = _separate(C21, C12, C03, n, q, max_probe_trials)
    C = (C0 + matmul(W, H, q)) % q
    frame = np.concatenate([W, C], axis=1)
    checks.update(probe_checks)
    stages["unique_complement"] = time.perf_counter() - clock

    # 4. The quadratic factor space and the multiplication tensor.
    clock = time.perf_counter()
    coeffs = C21.transpose(0, 2, 1).reshape(m * n, -1)
    reduced, pivots = rref(coeffs, q)
    Qspace = reduced[:len(pivots)]
    if len(Qspace) != n:
        raise StructureInconclusive("Quadratic coefficient space does not have dimension n")
    G = solve(Qspace.T, coeffs.T, q).T.reshape(m, n, n).transpose(0, 2, 1)
    reconstructed = matmul(G.transpose(0, 2, 1).reshape(m * n, n), Qspace, q)
    checks["quadratic_tensor_identity"] = np.array_equal(reconstructed, coeffs)
    stages["quadratic_space_and_tensor"] = time.perf_counter() - clock

    # 5. A right symmetriser of the tensor, then the quadratic modifier.
    clock = time.perf_counter()
    sym = recover_symmetrizer(G, q)
    stages["symmetrizer"] = time.perf_counter() - clock
    checks["symmetrizer_unique"] = sym["nullity"] == 1
    checks["symmetrizer_invertible"] = sym["checks"].get("generator_invertible", False)
    checks["exact_symmetry"] = sym["checks"].get("generator_exact_symmetry", False)
    if not all(checks[k] for k in ["symmetrizer_unique", "symmetrizer_invertible",
                                   "exact_symmetry"]):
        raise StructureInconclusive("Tensor symmetrizer certificate failed")
    clock = time.perf_counter()
    modifier = recover_modifier(G, C03, q, max_trials=max_probe_trials + 4 * n * n)
    stages["quadratic_modifier"] = time.perf_counter() - clock
    if modifier["status"] != "recovered":
        raise StructureInconclusive("Modifier reconstruction was not certified", modifier)
    checks["modifier_identity"] = modifier["checks"]["coefficient_identity"]
    checks["modifier_unique"] = modifier["checks"]["unique_polynomial_solution"]

    # 6. Rebuild the public coefficients from the recovered form, and check the
    #    linear involution that the construction leaves behind.
    clock = time.perf_counter()
    separated = join_bidegrees(C21, C03, n)
    if roundtrip:
        recovered_public = substitute_cubics(separated, inverse(frame, q), q)
        checks["public_coefficient_roundtrip"] = np.array_equal(P, recovered_public)
    stages["public_roundtrip"] = time.perf_counter() - clock
    # This involution is a structure of P, not the original secret matrix S.
    F_inv = inverse(frame, q)
    J = matmul(frame * np.array([-1] * n + [1] * n, dtype=np.int64), F_inv, q)
    checks["involution_square"] = np.array_equal(matmul(J, J, q), np.eye(N, dtype=np.int64))
    if not all(v for v in checks.values() if isinstance(v, (bool, np.bool_))):
        raise ArithmeticError("An exact public certificate failed")

    arrays = dict(W=W, C=C, frame=frame, Qspace=Qspace, G=G, C03=C03,
                  symmetrizer=sym["generator"], symmetric_space=sym["symmetric_matrices"],
                  modifier=modifier["r"], involution=J)
    info = dict(status="certified", q=q, n=n, m=m, **details,
                quadratic_space_dimension=n, symmetrizer_rank=sym["rank"],
                symmetrizer_nullity=sym["nullity"], checks=checks, stages_seconds=stages,
                modifier_probes=modifier.get("probe_counts"),
                modifier_uniqueness_witness=modifier.get("uniqueness_witness"),
                total_seconds=time.perf_counter() - started)
    return arrays, info
