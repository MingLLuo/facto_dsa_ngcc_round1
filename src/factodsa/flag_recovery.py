"""Recover the triangular flag of Q from the quadratic space U_Q.

The rank-one points come from msolve (rank1).  Peeling one rank-one direction at
a time, projecting onto the kernel of the corresponding linear form, and
backtracking over the (usually two) candidates per level recovers the true
triangular flag.  Success is decided by the certificate: extracting a triangular
family in the candidate flag must produce nonzero leading squares.
"""
from __future__ import annotations

import numpy as np

from . import fieldtools as FX
from .rank1 import rank_one_points


def project_kernel(M, direction, q):
    """Restrict every form to the hyperplane {x : direction . x = 0}."""
    kernel = FX.kernel(np.asarray(direction, dtype=np.int64).reshape(1, -1) % q, q)
    return [(kernel @ m @ kernel.T) % q for m in M], kernel


def recover_flag(M, q, *, seed=0, timeout=120, max_slices=2, verbose=False):
    """Return (L, Q_triangular, ok) with rows of L the recovered dual flag."""
    n = M[0].shape[0]

    def reconstruct(directions, kernels):
        rows = [directions[0]]
        accumulated = np.eye(n, dtype=np.int64)
        for level in range(1, n):
            kernel = kernels[level - 1]
            pinv = (kernel.T @ FX.inverse(kernel @ kernel.T, q)) % q
            accumulated = (accumulated @ pinv) % q
            rows.append((accumulated @ directions[level]) % q)
        return np.stack(rows, axis=0) % q

    def rec(current, level, directions, kernels):
        k = n - level
        for direction in rank_one_points(current, k, q, seed=seed + level,
                                         timeout=timeout, slices=max_slices):
            direction = np.array(direction, dtype=np.int64) % q
            if k == 1:
                L = reconstruct(directions + [direction], kernels)
                if FX.try_inverse(L, q) is None:
                    continue
                forms, ok = FX.extract_triangular_forms(M, L, q)
                if ok and all(int(forms[i][i, i]) % q for i in range(n)):
                    return L, forms
                continue
            next_forms, kernel = project_kernel(current, direction, q)
            result = rec(next_forms, level + 1, directions + [direction], kernels + [kernel])
            if result is not None:
                return result
        return None

    result = rec(M, 0, [], [])
    if result is None:
        return None, None, False
    if verbose:
        print(f"  flag recovered, diagonal = {[int(result[1][i][i, i]) % q for i in range(n)]}")
    return result[0], result[1], True
