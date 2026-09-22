# Positioning against the public NGCC reproduction harness

The `ngcc-harness` repository (`github.com/ngcc-dev/ngcc-harness`, commit
`26c8d47`) publishes reproductions for NGCC Round 1 findings.  Facto-DSA is
candidate `sign-10` there, and its whole record is **one** entry in
`security/vulnerabilities.csv`:

```
sign-10-1;Lead;review;The hidden zero subspace has an unpriced algebraic recovery path
```

Status `Lead`, verification `review`: a code/spec-review lead without a runtime
witness.  The repository also ships `sign-10/cryptanalysis/`, a self-contained
reimplementation of Facto-DSA (their Makefile states plainly that no submitter
source is included) with three programs: `selftest`, `attack1` and `attack3`.

## What that harness actually demonstrates

Measured by running its own targets on this machine:

| program | run | result | time |
|---|---|---|---|
| `selftest 10 13 1` | spec-faithful reimplementation at the standard 128-bit set | signing/verification round trip, tampered signature rejected, K2 dimension 10 as expected | instant |
| `attack1 5 8 1 10` | K2 = ker L2 from the public key | K2 recovered, equals the true `ker L2` | ~0 s |
| `attack1 6 10 1 10` | same at n=6 | K2 recovered, solving degree 7, quotient dimension 65 | 2.5 s |
| `attack3 6 10 1 12` | full chain at n=6 | K2, the (x,y) separation, the quadratic space Q and its symmetriser, and an equivalent **triangular central map**; their Algorithm-7 inversion succeeds on 27/200 random targets | ~4 s |

Its own summary is explicit about the stopping point:

> Remaining for a universal forger: the rational-normal-curve identification of
> the Hankel container (spec 3.2.2 Attack A step 3 / Attack D step 7, 'standard
> GRS recovery').

So the published state of the art it encodes is: recover K2 and an equivalent
triangular central map, **then still identify the Hankel container before a
forgery is possible**, demonstrated at `n = 5, 6` with an XL solving-degree
search, and without any check against the submitter's own verifier.

## What this repository adds

| dimension | ngcc-harness `sign-10` | this repository |
|---|---|---|
| scale | `n = 5, 6` (`m = 8, 10`) | **standard 128-bit set `n = 10, m = 13`** |
| structural recovery | K2 from a restricted system, XL solving-degree search | same target, but from the derivative space and its trace-annihilator, with exact coefficient certificates (every public coefficient rebuilt) |
| triangular flag | determinant filtration, `n = 6` | rank-one peeling with a backtracking certificate, `n = 10` (seconds to a minute) |
| step after the flag | needs the Hankel container / GRS recovery | **not needed**: the target equation is solved directly |
| forgery | not reached (their summary says so) | reached: `P(z) = h` verified coefficient by coefficient |
| verification against the authors' code | none (their tree excludes sign-10 sources) | **the authors' own `sig_verify` accepts the forged signature** for a KAT public key, with a flipped-bit negative control |
| complexity | solving degree observed 6-7 at `n = 5, 6` | `O(n^7)` structure + flag + `O(D^3)` with `D <= n^(m-n)`; plus the `f = m - n` analysis |
| parameter insight | none recorded | `f = m - n` decides reachability; the official 256/512 sets have `f = 15, 30`, while the authors' criteria admit `f = 1` sets |

## The one substantive difference

The harness's chain ends at "an equivalent triangular central map, invertible by
the submission's own Algorithm 7".  That is not yet a forger, because the
signer's flow needs the completed target to factor as `A * Y`, and recovering
that requires identifying the Hankel container (their words: standard GRS
recovery).

The route here never needs the Hankel container.  Once the separated form

```
P(B(x, y)) = Gamma(q(x) + r(y), y),      Gamma(a, y)_i = a^T G_i y
```

is certified, the target equation is posed directly: set `a = q(x) + r(y)`, get
a bilinear system in `(a, y)`, slice **only the y-block** and eliminate `a`,
leaving `f = m - n` equations of degree `<= n` in `f` unknowns.  With
`f = 3` for the standard set that is a 3-variable, degree-10 system - small
enough for `msolve` to return a rational parametrisation of it outright, and
the remaining `a - r(y) in image(q)` condition is met with the signer's own
rescaling.  The factorisation machinery that the harness says is still missing
is simply not on the path.

## Attack flow, in the harness's staging vocabulary

| stage | harness stage | what this repository does | cost |
|---|---|---|---|
| 1 | `attack1`: K2 from a restricted XL system | derivative space, trace-annihilator, support of the annihilator images | `O(n^7)` field ops, `O(n^5)` memory; 0.06 s at `n=10`, ~20 s at `n=32` |
| 2 | `attack3` stage 6: determinant filtration | rank-one peeling with a backtracking certificate | measured 4-76 s at `n=10`, 113-294 s at `n=14`; worst case `Theta(2^n)` |
| 3 | not present | y-only slice, eliminate `a`, `f = m - n` residual equations of degree `<= n` | `O(D^3)`, `D <= n^f`; 0.02 s when `f = 1`, 0.40 s for `f = 3` |
| 4 | not present | image condition plus rescaling | a handful of scalar multiplications per candidate |
| 5 | not present (they cannot test it) | encode `z`, call the authors' `sig_verify` | milliseconds |

## Reproducing the comparison

```sh
git clone --depth 1 https://github.com/ngcc-dev/ngcc-harness
cd ngcc-harness/sign-10/cryptanalysis && make all
./build/selftest 10 13 1
./build/attack1 6 10 1 10
./build/attack3 6 10 1 12          # stops at the triangular map
```

On this side, `scripts/run_standard_128.py` (`n = 10`, `m = 13`) and
`scripts/run_acceptance.py --source kat` produce the numbers quoted above.
