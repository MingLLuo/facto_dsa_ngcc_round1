"""Bezout-aware reduction of the recovered bilinear system.

The naive route slices the 20 unknowns with 7 generic linear forms, which gives
13 quadratics in 13 variables (Bezout 2^13) and a Macaulay blow-up.  Slicing only
the y-block instead keeps the block structure:

    y = Y0 + sum_j s_j Yj          (7 slices, s in F_q^{m-n})
    M(s) a = h                     (m x n, affine in s)

Solving the first n rows for a and substituting into the remaining m-n rows
leaves m-n equations of degree <= n in m-n unknowns.  For Facto-DSA-128 that is
3 equations of degree 10 in 3 variables: Bezout <= 1000, and the Macaulay size
is a few thousand monomials instead of 2e7.

The residual polynomials are recovered by evaluation + interpolation over F_q,
so only numeric determinants and linear solves are needed.
"""
from __future__ import annotations

import random

import numpy as np

from . import fieldtools as FX
from .msolve_iface import (format_polynomial, monomial_string, parse_msolve_points,
                           run_msolve)


def reduce_bilinear(G, h, q, n, m, *, seed=0):
    """Return a solver for the reduced system plus the y-parametrisation."""
    extra = 2 * n - m          # number of y slices
    free = n - extra           # = m - n, dimension of the reduced system
    if free <= 0:
        return None
    h = np.asarray(h, dtype=np.int64) % q
    rng = random.Random(seed)
    forms = np.array([[rng.randrange(1, q) for _ in range(n)] for _ in range(extra)],
                     dtype=np.int64)
    target = np.array([rng.randrange(q) for _ in range(extra)], dtype=np.int64)
    solution = FX.try_solve(forms, target, q)
    if solution is None:
        return None
    base, directions = solution, FX.kernel(forms, q)

    def y_of(s):
        return (base + np.asarray(s, dtype=np.int64) @ directions) % q

    # M(s)[i][p] = sum_b G[i,p,b] y(s)_b  -> affine in s
    g_flat = np.array(G, dtype=np.int64) % q
    base_vector = np.einsum("ipb,b->ip", g_flat, base) % q
    direction_vectors = np.einsum("ipb,jb->jip", g_flat, directions) % q

    def matrix_of(s):
        return (base_vector + np.einsum("j,jip->ip", np.asarray(s, dtype=np.int64),
                                       direction_vectors)) % q

    head = list(range(n))            # rows used to solve for a
    tail = list(range(n, m))         # rows used as residual equations

    def residual(s):
        M = matrix_of(s)
        top = M[head]
        inverse = FX.try_inverse(top, q)
        if inverse is None:
            return None
        a = (inverse @ h[head]) % q
        # clear the denominator det(M_head) so the residual is polynomial
        return (FX.det(top, q) * (M[tail] @ a - h[tail])) % q, a

    def det_of(s):
        return FX.det(matrix_of(s)[head], q)

    return dict(free=free, extra=extra, base=base, directions=directions,
                y_of=y_of, matrix_of=matrix_of, head=head, tail=tail,
                residual=residual, det_of=det_of)


def _interpolate(reduction, q, sample, seed, budget):
    """Recover one polynomial from a sampling function by evaluation + solving."""
    free = reduction["free"]
    monomials = FX.monomials_up_to(free, reduction["degree"])
    count = len(monomials)
    rng = random.Random(seed)
    rows, values = [], []
    attempts = 0
    while len(rows) < count + 5 and attempts < budget * count:
        attempts += 1
        point = [rng.randrange(q) for _ in range(free)]
        value = sample(point)
        if value is None:
            continue
        rows.append(FX.evaluate_monomials(monomials, point, q))
        values.append(np.asarray(value, dtype=np.int64) % q)
    if len(rows) < count:
        return None
    coefficients = FX.try_solve(np.asarray(rows), np.asarray(values), q)
    if coefficients is None:
        return None
    return dict(monomials=monomials, coefficients=coefficients,
                degree=reduction["degree"])


def interpolate_residuals(reduction, h, q, n, m, *, seed=0, degree=None):
    """Recover the m-n residual polynomials by evaluation + interpolation."""
    reduction = dict(reduction, degree=degree if degree is not None else n)

    def sample(point):
        result = reduction["residual"](point)
        return None if result is None else result[0]

    return _interpolate(reduction, q, sample, seed + 7, 40)


def interpolate_determinant(reduction, q, n, *, seed=0):
    """Recover det(M_head(s)) as a polynomial of total degree <= n."""
    reduction = dict(reduction, degree=n)

    def sample(point):
        determinant = int(reduction["det_of"](point)) % q
        return None if determinant == 0 else [determinant]

    result = _interpolate(reduction, q, sample, seed + 31, 60)
    if result is None:
        return None
    return dict(monomials=result["monomials"],
                coefficients=result["coefficients"].reshape(-1))


def _system_texts(polynomials, determinant_poly, q):
    """msolve input for the saturated reduced system: residuals + z*det - 1."""
    monomials = polynomials["monomials"]
    coefficients = polynomials["coefficients"]
    texts = [format_polynomial(
        [(int(coefficients[k, index]) % q, monomial_string(mono, "s") if mono else None)
         for k, mono in enumerate(monomials)])
        for index in range(coefficients.shape[1])]
    terms = []
    for k, monomial in enumerate(determinant_poly["monomials"]):
        value = int(determinant_poly["coefficients"][k]) % q
        if value:
            name = monomial_string(monomial, "s")
            terms.append((value, f"z*{name}" if name else "z"))
    terms.append((-1, None))
    return texts + [format_polynomial(terms)]


def solve_reduced(reduction, polynomials, q, n, m, *, timeout=300, seed=0):
    """Solve the reduced system with msolve, and lift the points it returns.

    Every residual carries a factor det(M_head), so the raw system also contains
    the spurious surface det = 0.  Saturation via the Rabinowitsch variable z
    with z*det - 1 = 0 keeps exactly the genuine solutions.

    msolve's own output for this system already *is* a univariate
    representation - an elimination polynomial in a separating variable plus
    the other coordinates as polynomials in it - so it is read directly.  There
    is no second solver: if msolve does not answer in that chart, the slice is
    reported as unsolved (`"none"`) and the caller moves on.
    """
    free = reduction["free"]
    determinant_poly = interpolate_determinant(reduction, q, n, seed=seed)
    if determinant_poly is None:
        return [], "none"
    variables = [f"s{j}" for j in range(free)] + ["z"]
    texts = _system_texts(polynomials, determinant_poly, q)

    def residuals(point):
        s = point[:free]
        values = np.array(FX.evaluate_monomials(polynomials["monomials"], s, q),
                          dtype=np.int64)
        determinant = int(np.array(
            FX.evaluate_monomials(determinant_poly["monomials"], s, q), dtype=np.int64)
            @ determinant_poly["coefficients"])
        return np.append((polynomials["coefficients"].T @ values) % q,
                         (int(point[free]) * determinant - 1) % q)

    report = {}
    output = run_msolve(variables, q, texts, timeout=timeout)
    points = parse_msolve_points(output, variables, q, residuals, report=report)
    method = "parametrization" if report.get("decoded") else "none"

    solutions = []
    for point in points:
        s = [int(v) % q for v in point[:free]]
        if reduction["det_of"](s) == 0:
            continue
        result = reduction["residual"](s)
        if result is None:
            continue
        _, a = result
        solutions.append((tuple(int(v) % q for v in a),
                          tuple(int(v) % q for v in reduction["y_of"](s))))
    return solutions, method
