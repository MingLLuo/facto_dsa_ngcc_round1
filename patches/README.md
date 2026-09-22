# Patches against the authors' reference implementation

## `0001-reference-functional-harness.patch`

Adds three files to a `Reference_Implementation` tree (the authors' own layout,
with `Facto-DSA-128/`, `-256/`, `-512/` subdirectories):

* `functional_test.c` - drives the authors' own DRNG plus `sig_keygen`,
  `sig_sign`, `sig_verify` for three message lengths (0, 1, 56 bytes) per
  parameter set, and checks that a valid signature is accepted while an
  appended message byte and a truncated signature are rejected;
* `check_functional.py` - compiles and runs that harness against every parameter
  set present in the tree (missing sets are skipped, so the vendored
  128-bit-only copy works), records return codes, lengths, source hashes and
  peak memory;
* `FUNCTIONAL_CHECK.md` - usage notes.

It does not modify any cryptographic source; it only adds a test driver.  The
patched files were verified to apply cleanly to a copy of `vendor/reference-128`
and to report `3/3 ordinary API cases passed` there.

```sh
cd <Reference_Implementation>
patch --dry-run -p1 < /path/to/0001-reference-functional-harness.patch
patch -p1 < /path/to/0001-reference-functional-harness.patch
python3 check_functional.py --source . --output /tmp/fresh-results
```

## Acceptance helper (not a patch)

`src/factodsa/verify_helper.c` plays the same role for the attack: it `#include`s
`SIG_AlgorithmInstance.c` so the static `sm3_digest`/`xof_field` are reachable,
and exposes the protocol target plus the authors' `sig_verify`.  It is compiled
by `src/factodsa/acceptance.py` against `vendor/reference-128` and needs no
modification of the submitted tree.

