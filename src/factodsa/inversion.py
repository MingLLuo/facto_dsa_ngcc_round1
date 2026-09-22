"""Inversion from the recovered form.

Recovered identity for z = frame @ (x, y):

    P(z)_i = Gamma_i(q(x) + r(y), y),   Gamma_i(a, y) = sum_{p,q} G[i,p,q] a_p y_q

So a target h is met by any pair (a, y) with Gamma(a, y) = h and a - r(y) in
the image of the quadratic map q.  This module solves the bilinear system for
(a, y), checks the image condition, and inverts q through the triangular flag.
"""
from __future__ import annotations

import numpy as np

from . import fieldtools as FX
from .polynomials import evaluate, monomials as mono, quadratic_coefficients


def bilinear_residuals(G, h, point, q, n):
    """Residuals of Gamma(a, y) = h at ``point = a + y``, for checking a solve."""
    a = np.asarray(point[:n], dtype=np.int64) % q
    y = np.asarray(point[n:2 * n], dtype=np.int64) % q
    value = np.einsum("p,ips,s->i", a, np.asarray(G, dtype=np.int64) % q, y) % q
    return (value - np.asarray(h, dtype=np.int64)) % q


def modifier_at(modifier, y, q, n):
    """r(y) for the recovered quadratic correction."""
    y = np.asarray(y, dtype=np.int64) % q
    quadratic = np.array([y[p] * y[s] % q for p, s in mono(n, 2)], dtype=np.int64)
    return (np.asarray(modifier, dtype=np.int64) % q @ quadratic) % q


def q_forms_in_dual_coordinates(Qspace, L, q, n):
    """Coefficients of q_i(L^{-1} t) in the monomial basis of t."""
    inverse = FX.inverse(L, q)
    index = {m: i for i, m in enumerate(mono(n, 2))}
    out = np.zeros((n, len(index)), dtype=np.int64)
    for i in range(n):
        full = np.zeros((n, n), dtype=np.int64)
        for k, (p, s) in enumerate(mono(n, 2)):
            coefficient = int(Qspace[i, k]) % q
            if coefficient:
                full = (full + coefficient * np.outer(inverse[p], inverse[s])) % q
        for u in range(n):
            out[i, index[(u, u)]] = int(full[u, u]) % q
            for v in range(u + 1, n):
                out[i, index[(u, v)]] = (int(full[u, v]) + int(full[v, u])) % q
    return out


def solve_linear_combination(basis_rows, target_row, q):
    """Coefficients c with sum_i c_i * basis_rows[i] = target_row, or None."""
    return FX.try_solve(np.asarray(basis_rows, dtype=np.int64).T,
                        np.asarray(target_row, dtype=np.int64), q)


def triangular_change(Qspace, L, Q_tri, q, n):
    """C with Q_tri(t) = C @ q(x) where t = L x."""
    dual = q_forms_in_dual_coordinates(Qspace, L, q, n)
    triangular = np.stack([quadratic_coefficients([m], q)[0] for m in Q_tri]) % q
    C = np.zeros((n, n), dtype=np.int64)
    for j in range(n):
        coefficients = solve_linear_combination(dual, triangular[j], q)
        if coefficients is None:
            return None
        C[j] = coefficients
    return C


def solve_map(Q_tri, C, L, w, q, n, limit=4096):
    """Solve q(x) = w using the triangular family: Q_tri(t) = C w, x = L^{-1} t."""
    target = (C @ np.array(w, dtype=np.int64)) % q
    found = []

    def rec(j, current):
        if len(found) >= limit:
            return
        if j == n:
            found.append(list(current))
            return
        a = int(Q_tri[j][j, j]) % q
        b = 0
        c = 0
        for i in range(j):
            b = (b + 2 * int(Q_tri[j][j, i]) * current[i]) % q
        for i in range(j):
            for k in range(i, j):
                coefficient = int(Q_tri[j][i, k]) % q
                if i == k:
                    c = (c + coefficient * current[i] * current[i]) % q
                else:
                    c = (c + 2 * coefficient * current[i] * current[k]) % q
        for value in FX.solve_quadratic(a, b, (c - int(target[j])) % q, q):
            rec(j + 1, current + [value])

    rec(0, [])
    inverse = FX.inverse(L, q)
    return [(inverse @ np.array(t, dtype=np.int64)) % q for t in found]


def verify_public(P, z, h, q):
    """Check every public cubic coefficient of P(z) against h."""
    value = evaluate(P, np.asarray(z, dtype=np.int64), q)[0] % q
    return bool(np.array_equal(value, np.asarray(h, dtype=np.int64) % q))
