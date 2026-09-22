"""Shared helpers over F_q.

Plain linear algebra is delegated to the FLINT backend in `finite_field`
(`rank`, `kernel`, `rref`, `solve`, `inverse`, `det`).  This module adds only
what the attack needs on top of it: scalar inverses that may fail, quadratic
roots, the trace annihilator of a matrix space, and the triangular-family
certificate.
"""
from __future__ import annotations

import itertools

import numpy as np

from .finite_field import det, inverse, kernel, solve  # noqa: F401


def inv(value, q):
    """Inverse of a scalar in F_q."""
    return pow(int(value) % q, -1, q)


def try_inverse(matrix, q):
    """Inverse of a square matrix, or None when it is singular."""
    try:
        return inverse(matrix, q)
    except ValueError:
        return None


def try_solve(matrix, rhs, q):
    """One solution of ``matrix x = rhs``, or None when there is none."""
    try:
        return solve(matrix, rhs, q)
    except ValueError:
        return None


def monomials_up_to(count, degree):
    """All exponent tuples in ``count`` variables of total degree <= degree."""
    return [m for d in range(degree + 1)
            for m in itertools.combinations_with_replacement(range(count), d)]


def evaluate_monomials(monomials, point, q):
    """Value of every monomial at ``point``."""
    values = []
    for monomial in monomials:
        value = 1
        for index in monomial:
            value = value * int(point[index]) % q
        values.append(value)
    return values


def solve_quadratic(a, b, c, q):
    """Roots in F_q of ``a x^2 + b x + c`` (the protocol uses q = 3 mod 4)."""
    a %= q
    b %= q
    c %= q
    if a == 0:
        if b == 0:
            return None if c == 0 else []
        return [(-c) * inv(b, q) % q]
    discriminant = (b * b - 4 * a * c) % q
    inverse_2a = inv(2 * a, q)
    if discriminant == 0:
        return [(-b) * inverse_2a % q]
    if pow(discriminant, (q - 1) // 2, q) != 1:
        return []
    root = pow(discriminant, (q + 1) // 4, q)
    return sorted({(-b + root) * inverse_2a % q, (-b - root) * inverse_2a % q})


def roots(coefficients, q):
    """All roots in F_q of a polynomial given low-to-high coefficients."""
    degree = len(coefficients) - 1
    while degree > 0 and int(coefficients[degree]) % q == 0:
        degree -= 1
    if degree <= 0:
        return []                      # constant or zero polynomial
    values = np.arange(q, dtype=np.int64)
    result = np.zeros(q, dtype=np.int64)
    for coefficient in reversed([int(c) % q for c in coefficients[:degree + 1]]):
        result = (result * values + coefficient) % q
    return [int(v) for v in np.flatnonzero(result == 0)]


def trace_annihilator(forms, q):
    """Trace-orthogonal complement of span(forms), as symmetric matrices.

    These are exactly the quadratic conditions whose common zeros are the
    rank-one points of the span: ``l l^T`` lies in the span if and only if
    ``l^T B l = 0`` for every returned ``B``.
    """
    n = forms[0].shape[0]
    columns = [(i, i) for i in range(n)] + [(i, j) for i in range(n)
                                            for j in range(i + 1, n)]
    pairing = np.array([[(int(M[i, j]) % q) * (1 if i == j else 2) % q
                         for (i, j) in columns] for M in forms], dtype=np.int64)
    basis = kernel(pairing, q)
    out = np.zeros((len(basis), n, n), dtype=np.int64)
    for row, coefficients in enumerate(basis):
        for value, (i, j) in zip(coefficients, columns):
            out[row, i, j] = out[row, j, i] = int(value) % q
    return out


def extract_triangular_forms(forms, flag, q):
    """Triangular basis of span(forms) in the flag given by the rows of ``flag``.

    Certificate used by the flag recovery: level ``i`` must contribute a form
    supported on the first ``i+1`` coordinates with a nonzero square there.
    """
    size = flag.shape[0]
    inverse_flag = try_inverse(flag, q)
    if inverse_flag is None:
        return None, False
    dual = np.stack([(inverse_flag.T @ M @ inverse_flag) % q for M in forms])
    family = []
    for level in range(size):
        constraints = np.array([[dual[k, row, col] for k in range(size)]
                                for row in range(size) for col in range(size)
                                if row > level or col > level],
                               dtype=np.int64).reshape(-1, size)
        candidates = kernel(constraints, q) if constraints.size else np.eye(size, dtype=np.int64)
        for coefficients in candidates:
            matrix = np.einsum("k,kij->ij", coefficients.astype(np.int64), dual) % q
            if matrix[level, level] % q:
                family.append(matrix * inv(matrix[level, level], q) % q)
                break
        else:
            return None, False
    return family, True
