"""Synthetic Facto-DSA instances, for testing a run end to end.

``generate`` returns a public coefficient matrix and the secret objects behind
it.  The returned oracle is used **only** to pick a target that is known to be
reachable, so that a run can be scored; the attack itself never sees it.
"""
from __future__ import annotations

from .finite_field import inverse, matmul, np, rank
from .polynomials import monomials, substitute_quadratics, variables_times_quadratics


def invertible(rng, n, q):
    """A uniformly random invertible ``n x n`` matrix over F_q."""
    while True:
        A = rng.integers(q, size=(n, n), dtype=np.int64)
        if rank(A, q) == n:
            return A


def generate(n, m, q, seed):
    """Build the public cubics of one instance, plus the secret map behind them.

    ``Q`` is triangular (``Q_i`` uses only ``x_i`` and lower), ``R`` is general
    quadratic, ``S`` splits ``z = (x, y)`` and ``T`` is the first ``m`` rows of
    an invertible output map.  The returned public map is
    ``P(z) = T * (Q(x) + R(y)) * y`` expanded in the public monomial basis.
    """
    if not (n >= 2 and n < m <= 2 * n - 1 and q > 3):
        raise ValueError("Expected the cubic model parameter range and q > 3")
    rng = np.random.default_rng(seed)
    pairs = monomials(n, 2)
    index = {pair: i for i, pair in enumerate(pairs)}

    Q = np.zeros((n, len(pairs)), dtype=np.int64)
    for i in range(n):
        for col, (a, b) in enumerate(pairs):
            if b <= i:
                Q[i, col] = rng.integers(q)
        Q[i, index[i, i]] = rng.integers(1, q)          # nonzero square
    R = rng.integers(q, size=Q.shape, dtype=np.int64)
    while True:
        T = rng.integers(q, size=(m, 2 * n - 1), dtype=np.int64)
        if rank(T, q) == m:
            break

    S = invertible(rng, 2 * n, q)
    Sinv = inverse(S, q)
    # Expand in public coordinates before anything is handed to the extractor.
    quadratic = (substitute_quadratics(Q, S[:n], q)
                 + substitute_quadratics(R, S[n:], q)) % q
    H = T[:, np.arange(n)[:, None] + np.arange(n)[None, :]]
    linear = matmul(H.reshape(m * n, n), S[n:], q).reshape(m, n, 2 * n)
    terms = matmul(linear.transpose(0, 2, 1).reshape(m * 2 * n, n), quadratic, q)
    P = variables_times_quadratics(terms.reshape(m, 2 * n, -1), q)
    oracle = dict(S=S, Sinv=Sinv, Q=Q, R=R, T=T, W=Sinv[:, :n], C=Sinv[:, n:])
    return P, oracle
