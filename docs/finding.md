# Facto-DSA (`sign-10`): signatures forged from the public key alone

**Result, in one sentence.** Using only a public key, we recover an equivalent
internal description of the scheme and produce a signature that the
submitter's own verification code accepts - at the standard 128-bit parameters.

## Why this matters

Facto-DSA's security rests on the public map being hard to invert.  We show the
map is not hard to *understand*: its shape can be recovered step by step, and
once recovered, the "hard" inversion collapses into a small algebraic problem.
Nothing secret is used at any point.

![Left panel: the derivative space is short by 45 dimensions](../figures/why-it-breaks.png)

The three steps of the attack - the 45-dimensional defect in the derivative
space, the recovery of a triangular basis, and the collapse of the target
equation to 3 unknowns - are set out in the repository
[`README.md`](../README.md) and in one page of plain language in
[`EXPLAINED.md`](../EXPLAINED.md); the right panel of the figure above shows the
reduction of the target system.

## Abstract for the site (about 120 words)

Facto-DSA hides a structured map behind two secret changes of coordinates.  But
the structure is not hidden in a measurable sense.  The space spanned by the
public map's first derivatives is too small by exactly `n(n-1)/2` dimensions,
and the missing directions reveal the split between the two groups of
variables; a linear-algebra computation on the public coefficients recovers it.
In those coordinates the secret shape becomes visible, and the target equation
collapses from 13 cubic equations in 20 unknowns to **3 equations of degree at
most 10 in 3 unknowns**.  A signature forged this way is accepted by the
authors' own verification code at the standard 128-bit parameters.

## What we verified

| witness | result |
|---|---|
| forged signature checked by the authors' `sig_verify`, KAT public key | **accepted** (return code 0), 40-byte signature from a 40040-byte key |
| same signature with one bit flipped | rejected (return code -1) |
| forged signature and structural claims re-checked with SageMath | all checks pass (`acceptance.md`) |
| standard 128-bit instance, end to end | forged in 13.8 s, peak 67 MiB |
| small-parameter regression (three seeds each, `n = 2..6`) | 15/15 forged, 0.07-1.32 s each |
| higher-target set admitted by the submitter's criteria (160-bit, `n=14,m=15`) | flag 114 s, solve 0.02 s, 3/3 valid |
| structural recovery at `n = 32` (the 512-bit set) | ~19 s, ~1.4 GiB, single run (19-24 s across runs) |
| official 256- and 512-bit sets | **not** broken: `f = m - n = 15` and `30` |

## Figures

| file | caption suggestion |
|---|---|
| `figures/why-it-breaks.png` | Left: the public derivative space is 45 dimensions short of the ambient 210, and the missing block names the hidden variable split.  Right: the same construction reduces the target equation to three unknowns. |
| `figures/pipeline.png` | The attack chain: public key -> structural recovery -> triangular flag -> reduced bilinear system -> parametrisation solve -> image condition -> triangular inversion -> verification. |
| `figures/stage-timings-128.png` | Where the time goes at the standard 128-bit parameters. |
| `figures/structure-scaling.png` | Measured structural-recovery time and memory for `n = 10, 14, 17, 32`. |
| `figures/bezout-vs-f.png` | Cost of the reduced system against `f = m - n`; the three `f = 1` sets admitted by the submitter's own criteria are marked. |

Start with `why-it-breaks.png`: it is the whole argument in two panels.  Then
`pipeline.png` shows the order of operations, and `stage-timings-128.png` shows
that no single stage is expensive in a surprising way.  `bezout-vs-f.png` is the
design lesson: the cost is set by `f = m - n`, which the submitter's parameter
criteria do not constrain.

## Reproduce

```sh
python3 scripts/check_environment.py          # prints READY
python3 scripts/run_standard_128.py           # the 128-bit instance, ~13 s
python3 scripts/run_acceptance.py --source kat # the authors' sig_verify accepts
python3 scripts/visualize.py                  # the figures above
```

## Inventory row

Matching the public harness format:

```
sign-10-2;Confirmed;runtime+static;The public key alone yields an equivalent triangular trapdoor and signatures accepted by the submitter's verifier
```

The existing `sign-10-1` entry ("The hidden zero subspace has an unpriced
algebraic recovery path") is a `Lead` at review level and stops at the
triangular map; `positioning.md` states the difference.

## Checklist before publishing

1. **No license is granted.** The first-party code and documentation ship
   without a license file, so no rights are granted beyond reading them here;
   the vendored material keeps its own terms (`NOTICE.md`).
2. **Attribution.** Keep `NOTICE.md` and `vendor/MANIFEST.sha256` with any
   published copy; they record exactly which submitted files were used.
3. **Claims to state carefully.** The 256- and 512-bit sets are *not* broken;
   the flag-recovery stage has a heuristic cost; the recovered object is an
   equivalent trapdoor, not the original secret basis.  Details in
   `limitations.md`.
4. **Cross-link the prior record.** `sign-10-1` in the public harness is a
   `Lead`/`review` for the same construction; `positioning.md` states the
   difference (scale, certificates, and the fact that this route does not need
   the Hankel container).
5. **Machine-dependence.** All timings come from one machine (8 cores, msolve
   0.10.1); quote them as measurements, not as benchmarks.
6. **Provenance.** Acknowledge the earlier material this grew out of and the
   assistants used while building it.  `NOTICE.md` carries the wording, and
   `README.md` has the longer version.

## Where to read more

| document | content |
|---|---|
| [`route.md`](route.md) | the scheme, and why the construction cannot hide its structure |
| [`attack.md`](attack.md) | the same chain stage by stage, with the code |
| [`complexity.md`](complexity.md) | costs, and the `f = m - n` parameter analysis |
| [`acceptance.md`](acceptance.md) | the verification against the authors' code, including the SageMath re-check |
| [`positioning.md`](positioning.md) | side-by-side against the public `ngcc-harness` reproduction |
| [`limitations.md`](limitations.md) | the complete list of what is and is not claimed |
