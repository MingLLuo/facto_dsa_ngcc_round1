# Recorded results

| file | what it contains |
|---|---|
| `small-validation.json` | 3/3 forgeries at `n = 2,3,4,5,6` with `m = 2n-1` (15/15) |
| `standard-128-seed0.json` | standard 128-bit instance, forged, with per-stage timings and the solver path used |
| `acceptance/acceptance.json` | the authors' `sig_verify` accepting a forged signature for a KAT public key |
| `higher-targets.json` | the 160-bit `f = 1` set (`n = 14`, `m = 15`) driven end to end |
| `structure-scaling.json` | structural recovery alone at the four sizes the docs quote (128, 160, 256, 512-bit sets) |

`acceptance/` also keeps the artifacts of that run: `pk.bin`, `message.bin`,
`signature.bin`, `signature_flipped.bin` (negative control) and the compiled
helper.

Numbers quoted in `README.md` and `docs/*` come from these files; they were
produced by the scripts in `scripts/` on the machine described in the README.

The `solve_method` field records how the reduced system was solved:
`parametrization` for msolve's own univariate representation, `none` when msolve
did not answer in that chart.  There is no second solver.

`scripts/verify_with_sage.py` (`make sage-check`) re-checks the forged
signature in `acceptance/` and the structural claims with SageMath, without
importing the attack code.  It needs Sage and the acceptance files above.
