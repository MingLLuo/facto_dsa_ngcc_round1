# Facto-DSA, explained in one page

## What the scheme does

Signing is easy if you know the shape of a certain secret map: a quadratic map
whose first coordinate depends only on the first variable, the second on the
first two, and so on - a "triangular" map.  So the designers build that shape and
then hide it behind two secret changes of coordinates, one on the inputs (`S`)
and one linear projection on the outputs (`T`):

```
public map:   P(z) = T * mu( Q(X) + R(Y), Y ),      (X, Y) = S z
```

`Q` is the triangular quadratic map, `R` a second quadratic map, `mu` polynomial
multiplication.  Verification just checks `P(z) = h` for a message-derived
target `h`.

## The phenomenon: the map has a fingerprint

The last step before the projection is *multiply by `Y`*.  That single fact
leaves a mark.

Look at all first derivatives of the public map.  They are quadratic forms, and
there are 210 of them in principle (in 20 variables).  But they only ever span
**165** dimensions.  The **45** that are missing are exactly the "`X`-squared"
directions:

```
210  all quadratic forms                    (everything possible)
165  what the derivatives actually reach
 45  missing = n(n-1)/2                     (only Q can make X^2, and Q is multiplied by Y)
```

That count belongs to the construction, not to the key: it can be measured on
any public key.

## Why that is a vulnerability

A space that is "too small" has a dual space that is "too big", and that dual
space is computable from the public key with ordinary linear algebra.  Its images
span exactly the hidden `X`-half of the variables:

> `S` is not secret in the way that matters.  What it hides is *which
> coordinates are which*, and that is precisely what the 45 missing dimensions
> reveal.

Once the split is known the rest follows:

1. the public map visibly separates into "quadratic in `x`, times `y`", plus
   "cubic in `y`";
2. that quadratic part *is* the triangular map, up to a change of basis;
3. triangularity is testable from the inside, so an equivalent triangular basis
   can be recovered - and the submitter's own inversion routine only needs *a*
   triangular basis, not the original one;
4. the equation a forger must solve collapses from 13 cubic equations in 20
   unknowns to **3 equations of degree at most 10 in 3 unknowns**.

![Why it breaks](figures/why-it-breaks.png)

## How we implemented it

| step | what it does | cost |
|---|---|---|
| 1 | derivatives plus their dual space, to recover the hidden `X`-half | linear algebra, 0.06 s |
| 2 | one linear system, to recover the hidden `Y`-half | linear algebra |
| 3 | rank-one peeling with a certificate, to recover a triangular basis | 13 s |
| 4 | slice only the `y` variables, eliminate, read msolve's parametrisation | 0.40 s |
| 5 | encode the preimage and hand it to the authors' verifier | instant |

## The result

On the standard 128-bit parameters (`q = 65519`, `n = 10`, `m = 13`) the whole
run takes **14 seconds and 67 MiB**, and the signature it produces is
**accepted by the authors' own `sig_verify`** (return code 0).  Flipping one bit
of that signature makes it reject (code -1).  The small-parameter regression
passes 15/15.

## What this does *not* say

* The 256- and 512-bit sets are **not** broken: their reduction leaves 15 and 30
  unknowns instead of 3, which is out of reach.  The deciding number is `m - n`,
  and the submitter's parameter criteria never mention it.
* We recover an *equivalent* triangular map, not the authors' original secret
  basis.
* The triangular-basis step is correct by certificate, but its running time is
  not provably polynomial.

## Reproduce

```sh
python3 scripts/check_environment.py          # prints READY
python3 scripts/run_standard_128.py           # forged, ~13 s
python3 scripts/run_acceptance.py --source kat # the authors' sig_verify accepts
```
