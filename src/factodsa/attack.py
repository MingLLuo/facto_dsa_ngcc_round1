"""End-to-end attack driver.

Public coefficients in, preimage of a given target out.  The chain is

    structural recovery -> triangular flag -> reduced bilinear system ->
    Bezout-aware solve -> image condition + rescaling -> triangular inversion

Nothing here reads a secret key; ``run_case`` builds a reachable target from the
model oracle only to be able to score a run.
"""
from __future__ import annotations

import time

import numpy as np

from .flag_recovery import recover_flag
from .inversion import modifier_at, solve_map, triangular_change, verify_public
from .model import generate
from .polynomials import monomials as mono, quadratic_matrices
from .recovery import extract_structure
from .structural_solve import interpolate_residuals, reduce_bilinear, solve_reduced


def recover_quadratic_space(Qspace, q, n):
    """The recovered quadratic space U_Q as symmetric matrices in x coordinates."""
    return [quadratic_matrices(Qspace[i], n, q) for i in range(n)]


def recover_structure(n, m, q, P, *, flag_timeout=120):
    """Recover everything that does not depend on the target.

    Returns a dictionary that is either ``status="recovered"`` - carrying the
    separated form, the triangular flag and the change of basis - or a failure
    status with the timings so far.  The result can be reused for many targets;
    only :func:`forge_target` depends on the target.
    """
    started = time.perf_counter()
    t_recover = time.perf_counter()
    arrays, info = extract_structure(P, q)
    recover_seconds = time.perf_counter() - t_recover

    G = arrays["G"] % q
    frame = arrays["frame"] % q
    modifier = arrays["modifier"] % q
    Qspace = arrays["Qspace"] % q

    t_flag = time.perf_counter()
    quadratic_space = recover_quadratic_space(Qspace, q, n)
    L, Q_tri, flag_ok = recover_flag(quadratic_space, q, timeout=flag_timeout)
    flag_seconds = time.perf_counter() - t_flag
    if not flag_ok:
        return dict(status="flag_recovery_failed", recover_seconds=recover_seconds,
                    flag_seconds=flag_seconds)

    C = triangular_change(Qspace, L, Q_tri, q, n)
    if C is None:
        return dict(status="change_of_basis_failed", recover_seconds=recover_seconds,
                    flag_seconds=flag_seconds)
    return dict(status="recovered", n=n, m=m, q=q, P=P, G=G, frame=frame,
                modifier=modifier, Qspace=Qspace, L=L, Q_tri=Q_tri, C=C, info=info,
                recover_seconds=recover_seconds, flag_seconds=flag_seconds,
                structure_seconds=time.perf_counter() - started)


def forge_target(structure, h, *, solve_seed=0, timeout=180, rescalings=24, slices=8):
    """Hit one target with an already recovered structure.

    Returns ``status="forged"`` with the preimage ``z``, or
    ``status="no_preimage_from_bilinear_solutions"`` when this target was not
    reached.  Different ``solve_seed`` values explore different slices and
    rescalings of the same target.
    """
    n, m, q, P = structure["n"], structure["m"], structure["q"], structure["P"]
    G, frame, modifier = structure["G"], structure["frame"], structure["modifier"]
    L, Q_tri, C = structure["L"], structure["Q_tri"], structure["C"]
    started = time.perf_counter()

    t_solve = time.perf_counter()
    pairs = []
    solve_method = "none"
    for attempt in range(max(1, slices)):
        seed_attempt = solve_seed + attempt
        reduction = reduce_bilinear(G, h, q, n, m, seed=seed_attempt)
        if reduction is None:
            continue
        polynomials = interpolate_residuals(reduction, h, q, n, m, seed=seed_attempt)
        if polynomials is None:
            continue
        found, solve_method = solve_reduced(reduction, polynomials, q, n, m,
                                            timeout=timeout, seed=seed_attempt)
        if found:
            pairs.extend(found)
            break
    solve_seconds = time.perf_counter() - t_solve

    # Gamma(lam*a, lam^-1*y) = Gamma(a, y): the signer's rescaling, which also
    # moves a - r(y) through the image of q.
    rng_scalings = np.random.default_rng(solve_seed + 991)
    scalings = [1] + [int(v) for v in rng_scalings.integers(1, q, size=rescalings)]
    tried = 0
    for a, y in pairs:
        for lam in scalings:
            tried += 1
            inverse = pow(int(lam), -1, q)
            scaled_a = (int(lam) * np.array(a, dtype=np.int64)) % q
            scaled_y = (inverse * np.array(y, dtype=np.int64)) % q
            w = (scaled_a - modifier_at(modifier, scaled_y, q, n)) % q
            for x in solve_map(Q_tri, C, L, w, q, n):
                z = (frame @ np.concatenate([x, scaled_y])) % q
                if verify_public(P, z, h, q):
                    return dict(status="forged", solve_seconds=solve_seconds,
                                solve_method=solve_method,
                                target_seconds=time.perf_counter() - started,
                                bilinear_solutions=len(pairs), solutions_tried=tried,
                                z=[int(v) for v in z])
    return dict(status="no_preimage_from_bilinear_solutions", solve_method=solve_method,
                solve_seconds=solve_seconds, bilinear_solutions=len(pairs),
                target_seconds=time.perf_counter() - started)


def forge_from_public(n, m, q, P, h, *, solve_seed=0, timeout=180, flag_timeout=120,
                      rescalings=24, slices=8):
    """Recover structure from public coefficients and hit the given target h."""
    started = time.perf_counter()
    structure = recover_structure(n, m, q, P, flag_timeout=flag_timeout)
    if structure["status"] != "recovered":
        return dict(n=n, m=m, total_seconds=time.perf_counter() - started, **structure)
    forged = forge_target(structure, h, solve_seed=solve_seed, timeout=timeout,
                          rescalings=rescalings, slices=slices)
    return dict(n=n, m=m, recover_seconds=structure["recover_seconds"],
                flag_seconds=structure["flag_seconds"],
                total_seconds=time.perf_counter() - started, **forged)


def run_case(n, m, q, seed, *, target_seed=0, solve_seed=0, timeout=180,
             flag_timeout=120, rescalings=24, slices=8):
    """Synthetic driver: build a reachable target from the oracle, then forge."""
    P, oracle = generate(n, m, q, seed)
    rng = np.random.default_rng(target_seed)
    x0 = rng.integers(q, size=n)
    y0 = rng.integers(q, size=n)
    z0 = (oracle["Sinv"] @ np.concatenate([x0, y0])) % q
    triples = mono(2 * n, 3)
    h = np.array([sum(int(P[i, k]) * int(z0[triples[k][0]]) % q
                      * int(z0[triples[k][1]]) % q * int(z0[triples[k][2]])
                      for k in range(len(triples))) % q for i in range(m)])
    return forge_from_public(n, m, q, P, h, solve_seed=solve_seed, timeout=timeout,
                             flag_timeout=flag_timeout, rescalings=rescalings,
                             slices=slices)
