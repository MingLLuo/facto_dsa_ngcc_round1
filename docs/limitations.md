# Limitations and honest scope

* **What is broken.** The standard 128-bit set, end to end: the authors'
  `sig_verify` accepts signatures produced from the public key alone.  The same
  chain reproduces 15/15 forgeries on small parameters.

* **What is not broken.** The official 256- and 512-bit sets.  After the y-only
  reduction they have `f = m - n = 15` and `30` free parameters, i.e. Bezout
  numbers near `2^61` and `2^150`.  This is a property of the parameter choice,
  not of the available hardware.

* **Flag recovery is the weak link in the implementation, not in the theory.**
  Its worst case is exponential in `n`; measured times are 4-76 s at `n = 10`
  and 114 s at `n = 14`.  Whether `n = 17` or `n = 22` (the `f = 1` sets for
  192- and 256-bit targets) stay practical has not been established here.

* **Recovered structure is not the original secret basis.** The chain certifies
  a separated representation and inverts it; it does not claim to output the
  authors' `S`, `T`, `Q`, `R` in their original coordinates.

* **The structural stage has public preconditions.** The reduction needs
  `dim D = d* = 3n(n+1)/2` and an annihilator support of dimension `n`, both of
  which are checked on the instance at hand.  An instance failing either is
  *rejected* (`StructureInconclusive`), not attacked; no claim is made about
  such keys.

* **The inversion depends on msolve's output convention.** The target system is
  solved by reading msolve's rational parametrisation, which needs msolve to
  pick a separating linear form in a single variable - what it did on every
  instance measured here.  There is no second solver for the other case: such a
  slice is reported as unsolved (`solve_method: none`) and the run moves on.

* **Scale of the end-to-end check.** Acceptance is against the authors' own
  `sig_verify` compiled from the vendored 128-bit reference tree, on the
  vendored 128-bit KAT public key.  The KAT file also carries 256-bit keys, but
  this route cannot reach them for the reason above, so there is no
  acceptance-level check at 256 or 512 bits.

* **Synthetic runs need an oracle to be scored.** `model.generate` builds
  instances together with their secret map, but it is used only to choose a
  target that is known to be reachable.  The attack path itself takes public
  coefficients and a target vector, and the acceptance path takes a real
  protocol target derived from the public key and the message.

* **Probabilistic steps.** The y-slice is random and may miss rational points
  (a few slices are retried), and the image condition is met through rescaling.
  Reported timings include that retrying.

* **No constant-factor claim.** Stage timings are single-machine measurements
  with msolve's Groebner computations in the loop; they are evidence for
  feasibility, not a benchmark.  Repeat runs of structural recovery agree to a
  few per cent up to `n = 17` but move by about 20 % at `n = 32`, and the peak
  memory readings vary by roughly the same amount.

* **Materials.** The vendored specification, reference implementation and KAT
  files are the authors' submission and are used for research; the KAT file
  contains test secrets that this code never reads (only the `PK` field is
  parsed).
