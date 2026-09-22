# Complexity of the implemented chain

Write `f = m - n` (the free dimension left by the y-only reduction) and let `D`
be the degree of the reduced ideal, `D <= n^f`.

| stage | cost | measured |
|---|---|---|
| structural recovery | `O(n^7)` field ops, `O(n^5)` memory (dense bound) | `n=10`: 0.06 s / ~50 MiB; `n=14`: 0.17 s / ~100 MiB; `n=17`: 0.57 s / ~155 MiB; `n=32`: ~19 s / ~1.4 GiB |
| triangular flag | worst case `Theta(2^n)` from the backtracking; each level solves a `k(k-1)/2`-equation rank-one system | `n=10`: 4-76 s; `n=14`: 114 s |
| reduced solve | msolve's F4 plus its parametrisation on a system of degree `D <= n^f`; our side is `O(D)` per root | `f=1`: 0.02 s; `f=3`: 0.40 s |

One rank-one call inside the flag stage costs ~0.14 s at `n = 10`, of which
`msolve` itself is ~0.02 s: the rest is decoding its parametrisation and
re-verifying every candidate point against the system.  So this stage is bound
by our own verification loop rather than by the solver, and the recursion is
serial.  Both are concrete optimisation targets for the larger `f = 1` sets.

## The parameter that decides everything

The four points in the table above are one run each
(`scripts/measure_structure_scaling.py`, records in
`results/structure-scaling.json`).  The times are reproducible to a few per
cent up to `n = 17`; at `n = 32` they moved between 19 s and 24 s across runs,
and the peak-memory readings vary by roughly 20 % on their own.

The official sets differ in `f`:

| set | n | m | f = m-n | Bezout of the reduced system |
|---|---|---|---|---|
| 128 | 10 | 13 | 3 | `10^3` |
| 256 | 17 | 32 | 15 | `17^15 ~ 2^61` |
| 512 | 32 | 62 | 30 | `32^30 ~ 2^150` |

So the official 256- and 512-bit sets are out of reach of this route for a
reason that no amount of hardware changes: the reduced system itself is too
large.  The 128-bit set is reachable only because `m - n` is small there.

## Higher targets under the authors' criteria

The authors' parameter filter (direct algebraic, one-vector, preimage) admits
sets with `m = n + 1`, i.e. `f = 1`:

| target | n | m | f | note |
|---|---|---|---|---|
| 160 | 14 | 15 | 1 | driven end to end: flag 114 s, solve 0.02 s |
| 192 | 17 | 18 | 1 | structure recovery 0.5 s; flag recovery pending |
| 256 | 22 | 23 | 1 | not yet tested |

For `f = 1` the solve stage degenerates to finding the roots of a single
univariate polynomial of degree `<= n`, which costs milliseconds.  The binding
stage then becomes the rank-one/peeling step, whose measured growth is a factor
of about **1.5 to 2.5 per increment of `n`**.

The recorded points are 1.4 s at `n = 6`, 4.9 s at `n = 8`, 25.9 s at `n = 10`
(an earlier revision of the code) and 4-76 s at `n = 10`, 113-294 s at `n = 14`
here.  The 3.5x and 5.3x ratios of the first series are *per two* steps, i.e.
1.9x and 2.3x per step; `n = 10 -> 14` gives 1.45-2.9x per step.  Single runs are
noisy, so extrapolate to an order of magnitude only: from the 114 s baseline
this puts the `f = 1` sets at roughly 5-26 min for `n = 17` and 0.5-35 h for
`n = 22`.

At `f = 3` the stage no longer depends on the solver growing: `msolve` builds
the parametrisation in 0.17 s and what remains is the interpolation of the
residuals from the recovered structure plus reading the roots.  The `256` and
`512` sets fail before any of this, because `D = n^f` itself is `2^61` and
`2^150`.

`figures/bezout-vs-f.png` plots this: the attack cost is governed by `f`, and
the authors' criteria do not constrain `f`.

For the same analysis compared with the public `ngcc-harness` reproduction of
Facto-DSA (`sign-10`), including the measured cost of its XL-based stages, see
`docs/positioning.md`.

## Reachability of other parameter sets

The three stages can be evaluated for any `(n, m)`:

* structure recovery: measured about `0.06 s * (n/10)^5`, bounded by `O(n^7)`;
  never the binding stage (19 s and 1.3 GiB even at `n = 32`);
* flag recovery: `t(14) * r^(n-14)` with `r` between 1.5 and 2.5 (above);
* reduced solve: F4 on `f` equations of degree `<= n`.  What decides its cost is
  the Macaulay scale rather than `D` itself; a usable proxy is the number of
  monomials up to the degree of regularity, `C(f*n, f)`.  Measured calibration:
  the 128-bit instance has `C(30,3) = 4.1e3` and solves in 0.17 s, while
  anything at or above about `1e7` is out of comfortable reach.

| `C(f*n, f)` | `n=10` | `n=17` | `n=22` | `n=32` |
|---|---|---|---|---|
| `f=1` | 10 | 17 | 22 | 32 |
| `f=2` | 1.9e2 | 5.6e2 | 9.5e2 | 2.0e3 |
| `f=3` | 4.1e3 | 2.1e4 | 4.6e4 | 1.4e5 |
| `f=4` | 9.1e4 | 8.1e5 | 2.3e6 | 1.1e7 |
| `f=5` | 2.1e6 | 3.3e7 | 1.2e8 | 8.2e8 |
| `f=6` | 5.0e7 | 1.3e9 | 6.6e9 | 6.4e10 |
| `f=8` | 2.9e10 | 2.4e12 | 1.9e13 | 4.1e14 |

| parameter set | n | m | f | flag (extrapolated) | reduced solve | verdict |
|---|---|---|---|---|---|---|
| official 128 | 10 | 13 | 3 | 4-76 s (measured) | `C(30,3) = 4.1e3`, 0.17 s | **forged**, 13 s |
| `f = 1`, 160-bit target | 14 | 15 | 1 | 113-294 s (measured) | one univariate polynomial | **forged** |
| `f = 1`, 192-bit target | 17 | 18 | 1 | 5-26 min | one univariate polynomial | ~minutes |
| `f = 1`, 256-bit target | 22 | 23 | 1 | 0.5-35 h | one univariate polynomial | ~hours |
| official 256 | 17 | 32 | 15 | 5-26 min | `C(255,15) ~ 1.5e24` | out of reach |
| official 512 | 32 | 62 | 30 | hours to years | `C(960,30) ~ 1.8e57` | out of reach |

Read off the design rule: `f <= 3` keeps the solve stage trivial and makes the
flag stage the binding cost, which is why the `f = 1` sets are the ones within
reach; pushing `m` out to `n + 5` or `n + 6` raises the proxy to `3.3e7` and
`1.3e9` at `n = 17`, past the reach of this route.  The flag column is
order-of-magnitude only: it starts from a single 114 s point at `n = 14`, and the
per-step factor is itself uncertain over 1.45-2.9.
