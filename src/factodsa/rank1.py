"""Rank-one points of a matrix space U = span(M_0, ..., M_{k-1}).

A point is a vector l with l l^T in U, equivalently l^T B l = 0 for every B in
the trace annihilator of U.  That homogeneous system is made zero-dimensional by
one random affine slice and handed to msolve; every parsed point is re-verified
against the system before it is returned.
"""
from __future__ import annotations

import random

import numpy as np

from .fieldtools import inv, trace_annihilator
from .msolve_iface import format_polynomial, monomial_string, parse_msolve_points, run_msolve


def form_terms(matrix, n, q):
    """(coefficient, monomial) terms of the quadratic form ``l^T matrix l``."""
    terms = []
    for i in range(n):
        for j in range(i, n):
            value = int(matrix[i, j]) % q
            if value:
                terms.append((value if i == j else 2 * value % q,
                              monomial_string((i, j), "l")))
    return terms


def _normalise(point, q):
    """Rescale a projective point so its first nonzero coordinate is one."""
    for value in point:
        if value % q:
            scale = inv(value, q)
            return tuple(int(v) * scale % q for v in point)
    return tuple(point)


def rank_one_points(M, n, q, *, seed=0, timeout=120, slices=2):
    """All projective rank-one points of span(M), as normalised tuples over F_q."""
    annihilator = trace_annihilator(M, q)
    variables = [f"l{i}" for i in range(n)]
    polynomials = [format_polynomial(form_terms(B, n, q)) for B in annihilator]

    def evaluate(point):
        vector = np.asarray(point, dtype=np.int64) % q
        return [int(vector @ B @ vector) % q for B in annihilator]

    rng = random.Random(seed)
    found = []
    for _ in range(slices):
        slice_terms = [(rng.randrange(1, q), f"l{i}") for i in range(n)]
        slice_terms.append((-1, None))
        system = polynomials + [format_polynomial(slice_terms)]
        output = run_msolve(variables, q, system, timeout=timeout)
        for point in parse_msolve_points(output, variables, q, evaluate):
            normalised = _normalise(point, q)
            if normalised not in found:
                found.append(normalised)
    return found
