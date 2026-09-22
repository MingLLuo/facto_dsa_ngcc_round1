# Attack chain, stage by stage

For why each of these stages must succeed - the invariants that make the scheme
attackable, and what would have to change to remove them - see `route.md`.

## Stage 1 - structural recovery (`src/factodsa/recovery.py`)

Input: the public coefficient matrix `P` (`m x C(2n+2,3)`).

1. Build the first-derivative space `D` (rows: output index times variable) and
   its right kernel.  For a well-formed instance the rank is exactly
   `d* = 3n(n+1)/2`, and the kernel dimension is `n(n-1)/2`.
2. Stack the kernel vectors as symmetric matrices (trace pairing, no division by
   two) and take the span of their images: that span is the hidden `X`-subspace
   `W = ker L2`.
3. Complete `W` to a basis `B0 = [W | C0]`, substitute `z = B0 (x, y)` and split
   the result by bidegree.  The pure `X^3` block must vanish.
4. Solve the `n` small linear systems that remove the mixed `x y^2` block; this
   yields a unique complement `C` and the frame `B = [W | C]`.
5. Extract the quadratic factor space `U_Q` from the `x^2 y` coefficients, the
   bilinear tensor `G` from `U_Q`, and the pure-`Y` cubic `C03`.  A simultaneous
   right symmetriser of `G` and a quadratic modifier `r` complete the picture.

Everything is certified by rebuilding each public coefficient exactly.

## Stage 2 - triangular flag (`src/factodsa/flag_recovery.py`, `rank1.py`)

`U_Q` is an `n`-dimensional space of symmetric matrices that is triangular in
some flag.  The attack:

1. computes the **rank-one points** `l` with `l l^T` in `U_Q`.  These are the
   common zeros of `l^T B l` over a basis of the trace-annihilator; one random
   affine slice makes that system zero-dimensional and `msolve` returns the
   points.  Every returned point is re-verified against the system;
2. peels one point at a time: project all forms onto the kernel of the
   corresponding linear form and repeat on the smaller space;
3. backtracks over the (usually two) candidates per level until the certificate
   `dim(U_Q cap Sym^2 F_i) = i+1` with nonzero leading square holds for every
   level.

## Stage 3 - target inversion (`src/factodsa/structural_solve.py`)

With the recovered frame, `P(B(x,y)) = Gamma(q(x) + r(y), y)` where
`Gamma(a,y)_i = a^T G_i y`.  Setting `a = q(x) + r(y)` turns the target equation
into a bilinear system of `m` equations in `2n` unknowns.

The slicing used here is deliberately asymmetric:

* slice **only the y-block** with `2n - m` generic linear forms, leaving
  `f = m - n` free parameters;
* substituting makes the system linear in `a`; solving the first `n` rows for
  `a` and substituting into the remaining rows leaves `f` residual equations of
  degree `<= n` in `f` unknowns.

For the standard 128-bit set `f = 3`, so this is 3 equations of degree `<= 10`
in 3 variables, Bezout `<= 1000` - small enough that `msolve` returns a
rational parametrisation of the solution set directly: one elimination
polynomial in a separating variable plus the other coordinates as polynomials
in it, which is exactly the univariate representation we want.
`src/factodsa/structural_solve.py` reads that output; nothing else about the
system changes.  This replaces the earlier route of asking `msolve` for a
reduced Groebner basis (`msolve -g 2`) and rebuilding the parametrisation here
with multiplication matrices: the solver was never the cost, our own extraction
was, and that stage went from 8.5 s to 0.40 s.  The recorded `solve_method`
field is `parametrization` when `msolve` answered in this chart and `none`
otherwise; there is no second solver.

Two corrections are essential:

* `msolve`'s parser mishandles repeated factors: writing `s0*s0` is read as
  `s0`.  Squares must be written `s0^2`.  Before this was found, every
  parametrisation looked "wrong" because a different system was being solved.
* every residual carries a factor `det(M_head)`, so the raw system also contains
  a spurious positive-dimensional `det = 0` component.  Adding the Rabinowitsch
  variable `z` with `z*det - 1 = 0` saturates it away.

## Stage 4 - image condition and rescaling (`src/factodsa/inversion.py`)

`a - r(y)` must be a value of `q`.  In the triangular coordinates this is a
coordinate-by-coordinate quadratic solve, which succeeds with moderate
probability; the bilinear system is invariant under
`(a, y) -> (lambda a, lambda^-1 y)` (the signer's rescaling), so a handful of
`lambda` values per sampled solution is enough to land inside the image.

## Stage 5 - acceptance (`src/factodsa/acceptance.py`)

`z` is encoded as `2n` little-endian `uint16` and passed to the authors'
`sig_verify`, compiled from the vendored sources through a helper that
`#include`s `SIG_AlgorithmInstance.c` so its static hash routines are reachable.
The same helper computes the protocol target `h` from the public key and the
message, so the forged signature is checked by the authors' own code path.
