"""Recover the quadratic correction r with ``Gamma(r(y), y) = C03(y)``.

The input describes ``Gamma(a, y)_i = sum_{a,b} G[i,a,b] a_a y_b``; only ``G``
and the pure-``y`` cubic ``C03`` are used.  The correction is expressed in these
public factor coordinates, not in any original key basis.

The points ``e_i`` and ``e_i + e_j`` (``i < j``) form a unisolvent stencil for
the quadratic monomials, and at any point where ``M_G(y)`` has full column rank
the equation ``M_G(y) r(y) = C03(y)`` is an ordinary least-squares-free linear
solve.  A rank-deficient stencil point is replaced by random full-rank points
from a deterministic stream.  One full-rank point certifies uniqueness of the
correction and the final check rebuilds every public cubic coefficient.
"""
from __future__ import annotations

from .finite_field import matmul, np, rank, rref, solve
from .polynomials import evaluate, monomials, variables_times_quadratics


def _field_array(value, q, name):
    a = np.asarray(value)
    if a.dtype.kind not in "iu":
        raise TypeError(f"{name} must contain integer field representatives")
    modulus = np.uint64(q) if a.dtype.kind == "u" else q
    return np.ascontiguousarray(a % modulus, dtype=np.int64)


def _quadratic_rows(points, pair_a, pair_b, q):
    """Values of every quadratic monomial at every point."""
    return (points[:, pair_a] * points[:, pair_b]) % q


def _independent_extension(row, echelon, pivots, q):
    """Reduce one row against pivot-normalized independent rows."""
    reduced = row.copy()
    if pivots:
        reduced = (reduced - matmul(row[pivots][None, :], echelon, q)[0]) % q
    nonzero = np.flatnonzero(reduced)
    if not len(nonzero):
        return None
    pivot = int(nonzero[0])
    reduced = reduced * pow(int(reduced[pivot]), -1, q) % q
    return reduced, pivot


def recover_modifier(G, C03, q, *, max_trials=10000):
    """Return a status dictionary; only ``recovered`` carries the correction.

    ``G`` has shape ``(m, n, n)`` and ``C03`` has shape
    ``(m, len(monomials(n, 3)))`` in the package's monomial order.  Modulus ``q``
    must be prime and within the backend bound ``2 <= q <= 65519``.
    ``max_trials`` bounds the random draws only, so there are at most
    ``n(n+1)/2 + max_trials`` point proposals.

    Status ``inconclusive`` means the budget ran out, not that no correction
    exists; ``incompatible`` means a point equation or the final coefficient
    identity failed, which certifies that none does.
    """
    if isinstance(q, (bool, np.bool_)) or not isinstance(q, (int, np.integer)):
        raise TypeError("q must be an integer prime modulus")
    q = int(q)
    if not 2 <= q <= 65519:
        raise ValueError("the finite-field backend requires 2 <= q <= 65519")
    if isinstance(max_trials, (bool, np.bool_)) or not isinstance(max_trials, (int, np.integer)):
        raise TypeError("max_trials must be a nonnegative integer")
    max_trials = int(max_trials)
    if max_trials < 0:
        raise ValueError("max_trials must be nonnegative")

    g = _field_array(G, q, "G")
    c03 = _field_array(C03, q, "C03")
    if g.ndim != 3 or g.shape[0] < 1 or g.shape[1] < 1 or g.shape[1] != g.shape[2]:
        raise ValueError("G must have shape (m,n,n), with m,n >= 1")
    m, n, _ = g.shape
    pairs = monomials(n, 2)
    s = len(pairs)
    if c03.shape != (m, len(monomials(n, 3))):
        raise ValueError("C03 has the wrong shape for the cubic monomial order")
    pair_a = np.asarray([a for a, _ in pairs], dtype=np.int64)
    pair_b = np.asarray([b for _, b in pairs], dtype=np.int64)
    counts = {
        "stencil_proposals": 0,
        "random_draws": 0,
        "rank_tests": 0,
        "full_column_rank": 0,
        "rank_deficient": 0,
        "dependent_quadratic_rows": 0,
        "consistent_point_solves": 0,
        "accepted_interpolation_points": 0,
    }
    witness = None
    method = "unisolvent_stencil"

    def finish(status, reason, correction=None, coefficient_identity=None,
               interpolation_rank=None, mismatches=None):
        return {
            "status": status,
            "reason": reason,
            "r": correction,
            "q": q,
            "shape": (m, n, n),
            "method": method,
            "max_trials": max_trials,
            "random_seed": 0,
            "probe_counts": dict(counts),
            "uniqueness_witness": witness,
            "checks": {
                "point_equations_exact": counts["consistent_point_solves"]
                    == counts["full_column_rank"],
                "interpolation_rank": interpolation_rank,
                "required_interpolation_rank": s,
                "coefficient_identity": bool(coefficient_identity),
                "full_cubic_coefficient_identity": coefficient_identity,
                "cubic_coefficient_mismatches": mismatches,
                "input_full_column_rank_probe": witness is not None,
                "unique_polynomial_solution": bool(
                    status == "recovered" and witness is not None
                    and coefficient_identity
                ),
            },
        }

    if m < n:
        return finish("inconclusive", "m < n prevents a full-column-rank point witness")

    eye = np.eye(n, dtype=np.int64)
    cross_pairs = [(a, b) for a, b in pairs if a != b]
    stencil = np.asarray(
        [eye[a] for a in range(n)] + [eye[a] + eye[b] for a, b in cross_pairs],
        dtype=np.int64,
    )
    stencil_rows = _quadratic_rows(stencil, pair_a, pair_b, q)
    g_flat = g.reshape(m * n, n)
    stencil_matrices = matmul(g_flat, stencil.T, q).reshape(m, n, s).transpose(2, 0, 1)
    stencil_targets = evaluate(c03, stencil, q)
    accepted_rows, accepted_values = [], []
    all_stencil_values = []

    def point_value(point, matrix, target):
        """Solve one point equation; (value, inconsistent) with value None if unusable."""
        nonlocal witness
        counts["rank_tests"] += 1
        if rank(matrix, q) != n:
            counts["rank_deficient"] += 1
            return None, False
        counts["full_column_rank"] += 1
        if witness is None:
            witness = {"point": point.tolist(), "column_rank": n}
        try:
            value = solve(matrix, target, q)
        except ValueError as error:
            if "Inconsistent" not in str(error):
                raise
            return None, True
        if not np.array_equal(matmul(matrix, value[:, None], q)[:, 0], target):
            raise ArithmeticError("point solve failed its exact public equation")
        counts["consistent_point_solves"] += 1
        return value, False

    # Try the unisolvent stencil first; it succeeds whenever no point is rank deficient.
    for index, point in enumerate(stencil):
        counts["stencil_proposals"] += 1
        value, inconsistent = point_value(point, stencil_matrices[index], stencil_targets[index])
        if inconsistent:
            return finish("incompatible", "C03(y) is outside the column image of M_G(y)",
                          coefficient_identity=False)
        all_stencil_values.append(value)
        if value is not None:
            accepted_rows.append(stencil_rows[index])
            accepted_values.append(value)
            counts["accepted_interpolation_points"] += 1

    if len(accepted_rows) == s:
        # The stencil is unisolvent, so the correction is read off directly:
        # the cross values are quadratic, not linear, in the probe directions.
        correction = np.empty((n, s), dtype=np.int64)
        diagonal_values = all_stencil_values[:n]
        cross_index = {pair: n + i for i, pair in enumerate(cross_pairs)}
        for col, (a, b) in enumerate(pairs):
            correction[:, col] = (
                diagonal_values[a] if a == b else
                (all_stencil_values[cross_index[(a, b)]]
                 - diagonal_values[a] - diagonal_values[b]) % q
            )
        interpolation_rank = s
    else:
        method = "random_full_rank_interpolation"
        if accepted_rows:
            echelon, pivots = rref(np.asarray(accepted_rows, dtype=np.int64), q)
            echelon = echelon[:len(pivots)].copy()
            if len(pivots) != len(accepted_rows):
                raise ArithmeticError("a subset of the unisolvent stencil lost row rank")
        else:
            echelon = np.empty((0, s), dtype=np.int64)
            pivots = []
        rng = np.random.default_rng(0)
        for _ in range(max_trials):
            if len(accepted_rows) == s:
                break
            counts["random_draws"] += 1
            point = rng.integers(q, size=n, dtype=np.int64)
            row = _quadratic_rows(point[None, :], pair_a, pair_b, q)[0]
            extension = _independent_extension(row, echelon, pivots, q)
            if extension is None:
                counts["dependent_quadratic_rows"] += 1
                continue
            matrix = matmul(g_flat, point[:, None], q).reshape(m, n)
            target = evaluate(c03, point, q)[0]
            value, inconsistent = point_value(point, matrix, target)
            if inconsistent:
                return finish("incompatible", "C03(y) is outside the column image of M_G(y)",
                              coefficient_identity=False, interpolation_rank=len(pivots))
            if value is None:
                continue
            normalized, pivot = extension
            # Keep identity columns at every stored pivot, so a future row can be
            # reduced in one modular matrix product.
            echelon = (echelon - echelon[:, pivot, None] * normalized[None, :]) % q
            echelon = np.concatenate([echelon, normalized[None, :]], axis=0)
            pivots.append(pivot)
            accepted_rows.append(row)
            accepted_values.append(value)
            counts["accepted_interpolation_points"] += 1
        if len(accepted_rows) != s:
            return finish(
                "inconclusive",
                f"random budget exhausted: only {len(accepted_rows)} of {s} "
                "independent full-rank points",
                interpolation_rank=len(accepted_rows),
            )
        interpolation_matrix = np.asarray(accepted_rows, dtype=np.int64)
        values = np.asarray(accepted_values, dtype=np.int64)
        interpolation_rank = int(rank(interpolation_matrix, q))
        if interpolation_rank != s:
            raise ArithmeticError("interpolation point certificate lost full rank")
        correction = solve(interpolation_matrix, values, q).T

    variable_quadratics = matmul(
        g.transpose(0, 2, 1).reshape(m * n, n), correction, q
    ).reshape(m, n, s)
    mismatches = int(np.count_nonzero(variables_times_quadratics(variable_quadratics, q) != c03))
    if mismatches:
        return finish("incompatible",
                      "interpolated values fail the complete cubic coefficient identity",
                      coefficient_identity=False, interpolation_rank=interpolation_rank,
                      mismatches=mismatches)
    return finish("recovered",
                  "complete cubic identity and a full-rank point certify a unique correction",
                  correction=correction, coefficient_identity=True,
                  interpolation_rank=interpolation_rank, mismatches=0)


__all__ = ["recover_modifier"]
