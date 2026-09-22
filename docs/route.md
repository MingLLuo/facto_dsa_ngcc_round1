# A new attack route: why Facto-DSA is attackable at all

The scheme hides a structured map behind two secret changes of coordinates.  It
fails because the *structure itself* is measurable from the public key: the
derivative space is too small, the missing directions name the hidden variable
split, and once the variables are split the target equation collapses to three
unknowns.  This document is the *why*; the procedure is in `attack.md` and the
numbers are in `complexity.md`.

## The scheme, in one page

Standard 128-bit set: `q = 65519`, `n = 10`, `r = 2n = 20`, `D = 2n-1 = 19`,
`m = 13`, `s = D - m = 6`; signature size `2r = 40` bytes; public key
`2 * m * C(r+2,3) = 40040` bytes.

* `S` in `GL_2n(F_q)` splits the inputs into `(X, Y) = S z`.
* `U` in `GL_D(F_q)`; `T` is its first `m` rows.
* `Q` is a homogeneous **triangular** quadratic map,
  `Q_i(X) = a_i X_i^2 + b_i(X_<i) X_i + c_i(X_<i)` with `a_i != 0`.
* `R` is a general homogeneous quadratic map.
* The public map is `P(z) = T * mu(Q(X) + R(Y), Y)`, where `mu` is polynomial
  multiplication truncated to degree `< D`.  It is a cubic map
  `F_q^20 -> F_q^13`.

**Signing.** The signer completes the target `h` to a degree-`<D` polynomial
`f = A * Y` with `deg A, deg Y < n`, computes `W = A - R(Y)` and solves
`Q(X) = W` coordinate by coordinate: triangularity makes each step a quadratic
in one unknown, with branching and a rescaling `lambda` that makes the
triangular preimage almost certain.

**Verification.** `sig_verify(pk, sigma, M)` decodes `z` from `sigma`, derives
the target `h = HashToField(FDSA-H || H_pk(pk) || M)`, evaluates the public cubic
at `z` and compares with `h`.  It is deterministic and uses public data only.

## Glossary

| term | meaning here |
|---|---|
| trace-annihilator | the set of symmetric matrices `B` that pair to zero with every derivative form; a purely linear-algebra object computed from the public coefficients |
| `W`, `ker L2` | the hidden subspace of inputs where the `Y`-part vanishes; recovering it is the same as recovering the variable split |
| bidegree | how many `x` factors and how many `y` factors a term has; the scheme only produces `x^2 y` and `y^3` |
| triangular flag | the nested subspaces in which the quadratic map becomes triangular; it is what makes the submitter's inversion routine work |
| `f = m - n` | the number of free parameters the reduction leaves; it decides the attack cost |

## The precise route

Everything here is stated for the submitted construction

```
P(z) = T * mu( Q(X) + R(Y), Y ),        (X, Y) = S z,     z in F_q^(2n)
```

with `Q` homogeneous triangular (`Q_i = a_i X_i^2 + b_i(X_<i) X_i + c_i(X_<i)`,
`a_i != 0`), `R` a general homogeneous quadratic map, `S` an invertible input
map, and `T` the first `m` rows of an invertible output map.

The design hides `S, T, Q, R`.  The claim of this route is that the composition
has three structural invariants that no choice of keys can remove, and that they
chain into an end-to-end break.

## 1. Invariant one: the derivative space is permanently defective

Differentiate the hidden expression.  Writing `U_Q = span(Q_0, ..., Q_{n-1})`:

```
d_X P  in  X* (x) Y*
d_Y P  in  U_Q  +  Sym^2(Y*)
```

so the whole first-derivative space `D` sits inside

```
E = U_Q  +  (X* (x) Y*)  +  Sym^2(Y*),
dim E = 3n(n+1)/2  =: d*.
```

The ambient space of quadratic forms in `2n` variables has dimension
`n(2n+1)`.  The difference

```
n(2n+1) - 3n(n+1)/2 = n(n-1)/2
```

is not a property of the key: it is a property of *writing a cubic as a
quadratic times `Y`*.  Any such construction is deficient by exactly that
much, and the deficiency lives in the `X (x) X` block.

**Why that exposes the input splitting.**  Let `A` be the trace-orthogonal
complement of `D` (a publicly computable linear-algebra object), and take the
sum of the images of the symmetric matrices in `A`:

> **Proposition (proved, not assumed).**  If the publicly checkable condition
> `dim D = d*` holds, then
> `sum_{B in A} im B = S^{-1}(X + 0) = ker L2`.

The proof is in `../docs/attack.md`'s reference material
(`research notes`, section 3, and `complement-theorem.md` in the wider project)
and it never assumes `Q` is generic: it uses `a_i != 0` to construct, for each
`i`, an explicit element of `A` whose image is a prescribed 2-plane inside `X`.
Consequently the images cover all of `X`, and conjugating by `S` gives the
hidden `X`-subspace.

**Consequence.**  The "input separation" is not a secret.  It is recovered by
linear algebra from the public coefficients - no root finding, no search.

**Why the design's own countermeasure does not help.**  The specification
introduces the `R(Y) * Y` term so that the derivative space is not of the pure
`X^2 + XY` type, which does defeat the Kipnis-Shamir style derivative chain.
But `R` contributes only `Y^2`, `XY` and mixed terms; *nothing except `Q`
produces `X^2`*.  The `n(n-1)/2`-dimensional `X^2` defect therefore survives
the countermeasure, and the trace-annihilator argument above does not need the
totally isotropic structure that the countermeasure removes.

## 2. Invariant two: the second half is fixed by *linear* equations

With `W = ker L2` known, complete it by any section and write
`z = B_X x + C_0 y`.  Because `P` vanishes on `W`, the expansion is

```
P(B_X x + C_0 y) = C21(x,y) + C12(x,y) + C03(y),
```

with the subscripts the bidegree in `(x, y)`.  Every other section using the
same quotient coordinates has the form `C_0 + B_X H`; the induced substitution
`x -> x + H y` sends the unwanted `x y^2` block to

```
C12(x,y) + D_x C21(x,y)[H y].
```

This is **linear in the entries of `H`**, and the hidden `Y`-space guarantees a
solution exists.  Uniqueness is a polynomial-identity argument: if two
solutions differ by `K`, then `D_x C21[K y] = 0` identically, which forces every
`(K y)_a y_b` to vanish and hence `K = 0`.

So the second half of the splitting is determined by solving one linear system,
and the resulting complement is exactly `ker L1 = S^{-1}(0 + Y)`.

**Consequence.**  After stages 1 and 2 - both pure linear algebra over `F_q` -
the public map is *literally* bidegree-separated:

```
P_i(x, y) = sum_{a,b} G_{i,ab} q_a(x) y_b + C03_i(y).
```

That is the whole secret shape of the central map, obtained without any
algebraic solving.

## 3. Invariant three: triangularity is a property of a space, not of a key

Stage 2 also hands over `U_Q = span(Q)`: it is the space spanned by the `x^2 y`
coefficients, and its dimension is `n` (a publicly checkable condition).

The signer's advantage is *not* "knows the map `Q`".  It is "knows a flag in
which `Q` is triangular", because that is what makes Algorithm 7 a
coordinate-by-coordinate quadratic solve.  Triangularity of `Q` with respect to
a flag `F_1 < F_2 < ... < F_n` is equivalent to the *internal* conditions

```
dim( U_Q  cap  Sym^2(F_i) ) = i + 1      for every i,
```

with a nonzero new square coefficient at each level.  These conditions mention
only `U_Q` and the flag, never the hidden basis, so they are invariants of the
subspace `U_Q` up to the `GL_n` freedom that remains.

The route computes the flag from those conditions:

1. rank-one points of `U_Q` - the `l` with `l l^T in U_Q`.  They are the common
   zeros of `l^T B l` over a basis of the trace-annihilator; one affine slice
   makes that a zero-dimensional system for `msolve`, and every returned point is
   re-verified;
2. peel: project onto the kernel of the corresponding linear form and repeat on
   the smaller space;
3. backtrack over the (usually two) candidates per level until the certificate
   above holds at every level.

**Consequence.**  An equivalent triangular central map is recovered, and the
submission's own Algorithm 7 inverts it.  The trapdoor is a *recoverable
property of a public space* rather than a secret function.

This is also the honest weak point of the implementation: the backtracking has
worst-case cost `Theta(2^n)`.  The certificate makes a *wrong* flag fail loudly,
which is why small backtracking suffices in practice, but no polynomial bound is
claimed for this stage.

## 4. Invariant four: the target equation collapses to `f = m - n` unknowns

With the separated form, set `a = q(x) + r(y)` (using the recovered `q` and the
recovered modifier `r`).  The public target equation becomes

```
Gamma(a, y)_i = a^T G_i y = h_i ,      i = 1..m,
```

a system that is **bilinear** in `(a, y) in F_q^n x F_q^n`: `m` equations in
`2n` unknowns.  Its solution set is `2n - m` dimensional.  Two facts make the
inversion cheap.

**(a) Slice only the y-block.**  Cutting `y` by `2n - m` generic linear forms
leaves `f := m - n` free parameters and keeps the system linear in `a`.  Solving
the first `n` rows for `a` and substituting into the remaining rows gives `f`
residual equations in `f` unknowns, and their degree is at most `n`: the
residual is `(M_tail adj(M_head) h_head)/det - h_tail`, and clearing the
`det` denominator (degree `n`) leaves a numerator of degree at most
`1 + (n-1) = n`.  Hence

```
Bezout of the reduced system  <=  n^f.
```

For the standard 128-bit set `f = 3`, i.e. three equations of degree at most 10
in three unknowns.  Compare the naive formulation of the same target: `m = 13`
cubic equations in `2n = 20` unknowns, whose Macaulay cost is what the
specification itself prices at "129 bits" - the reduced system is not a constant
factor better, it is a different problem.

**(b) The image condition is not an obstruction.**  A solution of the bilinear
system is only usable when `a - r(y)` lies in the image of `q`.  But the map
`(x, y) -> (q(x) + r(y), y)` covers the valid part of a 7-dimensional solution
set, and the bilinear system is invariant under

```
(a, y) -> (lambda a, lambda^-1 y),
```

which is exactly the rescaling the signer uses.  It moves `a - r(y)` through
the image of `q` with constant probability, so a few scalars per sampled
solution suffice.

## 5. Why the specification's own defences miss this

* **The derivative countermeasure is aimed at the wrong object.**  Section 3.2.3
  of the specification argues that `R(Y)Y` removes the common totally isotropic
  subspace and so blocks the Kipnis-Shamir start and the `Q(X)Y` chain.  That is
  correct for those chains, but the `X^2` dimension defect is invariant under the
  construction, and the trace-annihilator argument does not use isotropy.

* **The security criteria never price `f = m - n`.**  The three filters used to
  choose parameters constrain `n` and `m` through the direct-algebraic estimate,
  the one-vector search and the preimage search.  All of those are functions of
  `(n, m)` *and of the solving model*, while the attack cost here is governed by
  `f`.  At a fixed target level `f` can be as small as 1 and as large as `n-1`;
  the official 128-bit set happens to have `f = 3`, while the official 256- and
  512-bit sets have `f = 15` and `f = 30`, which is why they are out of reach.
  The same criteria admit `f = 1` sets at 160, 192 and 256 bits, and there the
  entire inversion stage degenerates to the roots of one polynomial of degree
  `<= n`.

* **"Direct algebraic solving" is not what is happening.**  The specification's
  estimate `d_solv = 2m+1` describes F4/F5 on the raw system.  The route above
  never computes a Groebner basis of the raw system; it removes the structural
  degrees of freedom first and leaves a system small enough that a Groebner
  basis (or any equivalent) is cheap.

## 6. What would actually remove the attack

These are the design consequences of the route, stated as countermeasures
rather than as claims:

1. **Make `f = m - n` large.**  Choosing `m` close to `2n-1` (large `s = D-m`
   is not the same thing - `f` is what matters) pushes the reduced system to
   `n^(n-1)`, which is far out of reach.  The official 256- and 512-bit sets
   already do this incidentally; the 128-bit set does not.
2. **Break the bidegree separation.**  The attack's stage 2 is exactly the
   recovery of the `X/Y` splitting; a construction whose public form is not
   `quadratic x Y` would not admit it.  But the separation is precisely what
   makes the signer's completion-plus-factorisation work, so this cannot be
   changed without replacing the trapdoor.
3. **Make the derivative space full-dimensional.**  Impossible while the public
   map is a product of a quadratic map with `Y`: the defect is exactly
   `n(n-1)/2`.

## 7. Scope of the claims

Stages 1-2 and the degree bound in stage 5 are proved statements under the
publicly checkable conditions stated with them, and the implementation checks
those conditions on the instance at hand.  Stage 4 is an exact algorithm with a
heuristic cost profile: the certificate decides correctness, and no polynomial
bound is claimed for it.  The 256- and 512-bit sets are **not** broken by this
route, because `f = 15` and `30` push the reduced system out of reach, as
explained above.

The complete list of what is and is not claimed is in
[`limitations.md`](limitations.md); the runtime evidence is in
[`acceptance.md`](acceptance.md).
