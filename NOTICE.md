# Provenance and intended use

## Authors' material (vendored, unmodified)

`vendor/spec/algorithm-specification-facto-dsa.pdf`,
`vendor/reference-128/*` and `vendor/kat/KAT_SIG_Facto-DSA-128.txt` are the
submitted Facto-DSA materials (ICCS / NGCC submission, June 2026).  They are
vendored here unchanged so that the reproduction is self-contained; the source
paths and SHA-256 digests of the vendored files, and of the first-party code
that ships alongside them (`src/factodsa/`, `patches/`), are recorded in
`vendor/MANIFEST.sha256`.

The KAT file is the authors' own test-vector file.  It contains test secrets as
well as public keys; only the `PK` field is ever read by the code here
(`src/factodsa/acceptance.py`).  Nothing in this repository uses a secret field.

## First-party code and provenance

`src/factodsa/*`, `scripts/*` and `docs/*` are the research implementation and
documentation produced for this reproduction.  They operate on public keys only:
the attack path takes public coefficients plus a target vector, and the
acceptance path additionally uses the authors' verification function.

The work grew out of the same working environment, which already held an
analysis of a previous version of Facto-DSA and a set of implementation
sketches; that earlier material pointed at the attack angle pursued here, and it
is preserved unchanged alongside this folder.  Parts of the mathematical route
were shaped by ideas from **GPT-6 Astra**.  The code and documentation in this
repository were produced with **DeepSeek Flash**, an agentic coding model,
working from that earlier material.  This is unofficial work, not connected to
the submitters of Facto-DSA beyond the public materials vendored above.

Every reported result is a recorded run of the scripts in `scripts/`, and every
claim is decided by an exact certificate or by the authors' own verification
code.

## Intended use

This is cryptanalysis of a published algorithm submission.  It is meant to
document that the standard 128-bit instance is forgeable, to characterise the
parameter families that are and are not reachable, and to support the authors'
parameter selection.  It is not a deployed-system exploit and no real key is
involved.
