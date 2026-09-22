"""Canonical homogeneous coefficients, exact derivatives and substitutions.

Every polynomial here is a coefficient vector indexed by ``monomials(n, degree)``:
the exponent tuples of total degree ``degree`` in ``n`` variables, in
lexicographic order (so ``a <= b <= c`` for degrees 2 and 3).  Substitution is
then a dense tensor contraction, and all of it stays exact over F_q.
"""
from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import combinations_with_replacement, permutations

from .finite_field import matmul, np


@lru_cache(None)
def monomials(n, degree):
    """Exponent tuples of total degree ``degree`` in ``n`` variables, sorted."""
    return tuple(combinations_with_replacement(range(n), degree))


def quadratic_matrices(coefficients, n, q):
    """Symmetric matrices of the forms whose coefficients are given.

    A monomial ``x_a x_b`` with ``a != b`` appears twice in ``x^T A x``, so the
    off-diagonal entries carry a factor 1/2.
    """
    C = np.asarray(coefficients, dtype=np.int64)
    A = np.zeros(C.shape[:-1] + (n, n), dtype=np.int64)
    half = pow(2, -1, q)
    for k, (a, b) in enumerate(monomials(n, 2)):
        A[..., a, b] = C[..., k] if a == b else C[..., k] * half % q
        A[..., b, a] = A[..., a, b]
    return A


def quadratic_coefficients(matrices, q):
    """Inverse of ``quadratic_matrices``: symmetric matrices back to coefficients."""
    A = np.asarray(matrices, dtype=np.int64)
    n = A.shape[-1]
    return np.stack([A[..., a, b] if a == b else (A[..., a, b] + A[..., b, a]) % q
                     for a, b in monomials(n, 2)], axis=-1)


def substitute_quadratics(C, B, q):
    """Coefficients of the quadratic map ``C`` in variables ``x_old = B x_new``."""
    N, n = B.shape
    A = quadratic_matrices(C, N, q)
    batch = A.shape[:-2]
    tmp = matmul(A.reshape(-1, N), B, q).reshape(-1, N, n)
    out = np.stack([matmul(B.T, t, q) for t in tmp])
    return quadratic_coefficients(out.reshape(batch + (n, n)), q)


def derivatives(P, N, q):
    """First derivatives: one row per ``(output, variable)``, one column per quadric."""
    m = len(P)
    pairs = monomials(N, 2)
    index = {pair: j for j, pair in enumerate(pairs)}
    D = np.zeros((m * N, len(pairs)), dtype=np.int64)
    rows = np.arange(m) * N
    for col, monomial in enumerate(monomials(N, 3)):
        for v, multiplicity in Counter(monomial).items():
            rest = list(monomial)
            rest.remove(v)
            D[rows + v, index[tuple(rest)]] = (P[:, col] * multiplicity) % q
    return D


def annihilator_matrices(K, N, q):
    """Symmetric matrices of the trace annihilator, from right-kernel coordinates.

    The pairing is ``sum_ij B_ij M_ij``, so the kernel coordinates are the
    symmetric entries directly - no factor 1/2, unlike ``quadratic_matrices``.
    """
    A = np.zeros((len(K), N, N), dtype=np.int64)
    for col, (a, b) in enumerate(monomials(N, 2)):
        A[:, a, b] = K[:, col]
        A[:, b, a] = K[:, col]
    return A % q


def variables_times_quadratics(B, q):
    """Coefficients of ``x_j * q_k(x)``, a linear map times a quadratic map."""
    n = B.shape[-2]
    triples = monomials(n, 3)
    index = {t: i for i, t in enumerate(triples)}
    out = np.zeros(B.shape[:-2] + (len(triples),), dtype=np.int64)
    for j in range(n):
        for k, pair in enumerate(monomials(n, 2)):
            target = index[tuple(sorted(pair + (j,)))]
            out[..., target] = (out[..., target] + B[..., j, k]) % q
    return out


def substitute_cubics(P, B, q):
    """Coefficients of a cubic map in variables ``z_old = B z_new``.

    The symmetric tensor is filled by monomial orbit, then contracted one axis at
    a time.  Separating the monomials needs characteristic other than 2 or 3.
    """
    N, n = B.shape
    if q in (2, 3):
        raise ValueError("Cubic tensor substitution requires q > 3")
    tensor = np.zeros((len(P), N, N, N), dtype=np.int64)
    for k, t in enumerate(monomials(N, 3)):
        orbit = set(permutations(t))
        value = P[:, k] * pow(len(orbit), -1, q) % q
        for a, b, c in orbit:
            tensor[:, a, b, c] = value
    # Contract one axis at a time, reducing modulo q after every multiplication.
    tensor = matmul(tensor.reshape(-1, N), B, q).reshape(len(P), N, N, n)
    tensor = matmul(tensor.transpose(0, 1, 3, 2).reshape(-1, N), B, q).reshape(len(P), N, n, n)
    tensor = matmul(tensor.transpose(0, 2, 3, 1).reshape(-1, N), B, q).reshape(len(P), n, n, n)
    out = np.empty((len(P), len(monomials(n, 3))), dtype=np.int64)
    for k, t in enumerate(monomials(n, 3)):
        out[:, k] = tensor[:, t[0], t[1], t[2]] * len(set(permutations(t))) % q
    return out


def split_bidegrees(P, n):
    """Split a public cubic into its ``x^3``, ``x^2 y``, ``x y^2`` and ``y^3`` blocks.

    Variables are split as ``z = (x, y)`` with ``x`` first; each returned array
    carries the coefficients of one bidegree.
    """
    index = {t: i for i, t in enumerate(monomials(2 * n, 3))}
    pairs = monomials(n, 2)
    triples = monomials(n, 3)
    m = len(P)
    C30 = np.stack([P[:, index[t]] for t in triples], axis=1)
    C21 = np.empty((m, len(pairs), n), dtype=np.int64)
    C12 = np.empty((m, n, len(pairs)), dtype=np.int64)
    for k, (a, b) in enumerate(pairs):
        for j in range(n):
            C21[:, k, j] = P[:, index[(a, b, n + j)]]
            C12[:, j, k] = P[:, index[(j, n + a, n + b)]]
    C03 = np.stack([P[:, index[tuple(n + a for a in t)]] for t in triples], axis=1)
    return C30, C21, C12, C03


def join_bidegrees(C21, C03, n):
    """Rebuild a public cubic from its ``x^2 y`` and ``y^3`` blocks."""
    index = {t: i for i, t in enumerate(monomials(2 * n, 3))}
    P = np.zeros((len(C21), len(index)), dtype=np.int64)
    for k, (a, b) in enumerate(monomials(n, 2)):
        for j in range(n):
            P[:, index[(a, b, n + j)]] = C21[:, k, j]
    for k, t in enumerate(monomials(n, 3)):
        P[:, index[tuple(n + a for a in t)]] = C03[:, k]
    return P


def evaluate(P, points, q):
    """Evaluate the cubic map ``P`` at one point or a stack of points."""
    points = np.asarray(points, dtype=np.int64)
    if points.ndim == 1:
        points = points[None, :]
    mon = np.stack([(points[:, a] * points[:, b] % q) * points[:, c] % q
                    for a, b, c in monomials(points.shape[1], 3)], axis=0)
    return matmul(P, mon, q).T
