"""Reference-implementation acceptance test.

Forge a signature for the *authors' own* verification code and check that the
original ``sig_verify`` accepts it.  Only the public key is used.

    python3 -m factodsa.acceptance --source kat      # a submission KAT public key
    python3 -m factodsa.acceptance --source model    # a synthetic 128-bit key
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import struct
import subprocess
import time
from pathlib import Path

import numpy as np

from .attack import forge_from_public
from .model import generate
from .polynomials import monomials as mono

REPO = Path(__file__).resolve().parents[2]
DEFAULT_REFERENCE = REPO / "vendor" / "reference-128"
DEFAULT_KAT = REPO / "vendor" / "kat" / "KAT_SIG_Facto-DSA-128.txt"
REFERENCE = DEFAULT_REFERENCE


def build_helper(output: Path) -> Path:
    compiler = os.environ.get("CC", "cc")
    if shutil.which(compiler) is None:
        raise SystemExit(f"C compiler '{compiler}' not found; set CC to one that exists")
    helper = output / "verify_helper"
    command = [compiler, "-std=c99", "-O2", "-I", str(REFERENCE),
               str(Path(__file__).resolve().parent / "verify_helper.c"),
               str(REFERENCE / "auxfunc.c"), str(REFERENCE / "drng.c"),
               "-o", str(helper)]
    subprocess.run(command, check=True, capture_output=True, text=True)
    return helper


def first_kat_public_key(path: Path) -> bytes:
    """Extract only the first PK field from a KAT text file."""
    pattern = re.compile(rb"^PK\s*=\s*([0-9A-Fa-f]+)\s*$", re.MULTILINE)
    with path.open("rb") as stream:
        head = stream.read(4 << 20)
    matches = pattern.findall(head)
    if not matches:
        raise SystemExit(f"no PK field found in {path}")
    return bytes.fromhex(matches[0].decode("ascii"))


def encode_signature(z, q):
    return b"".join(struct.pack("<H", int(v) % q) for v in z)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=["model", "kat", "pkfile"], default="kat")
    parser.add_argument("--n", type=int, default=10)
    parser.add_argument("--m", type=int, default=13)
    parser.add_argument("--q", type=int, default=65519)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--message", default="Facto-DSA acceptance test")
    parser.add_argument("--pk-file", type=Path)
    parser.add_argument("--kat", type=Path, default=DEFAULT_KAT)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--output", type=Path, default=REPO / "results" / "acceptance")
    parser.add_argument("--timeout", type=int, default=1800)
    args = parser.parse_args()

    global REFERENCE
    REFERENCE = args.reference
    args.output.mkdir(parents=True, exist_ok=True)
    q, n, m = args.q, args.n, args.m
    terms = len(mono(2 * n, 3))

    if args.source == "model":
        P, _oracle = generate(n, m, q, args.seed)
        P = np.asarray(P, dtype=np.int64) % q
        pk = b"".join(struct.pack("<H", int(v) % q) for v in P.reshape(-1))
    else:
        pk = args.pk_file.read_bytes() if args.source == "pkfile" \
            else first_kat_public_key(args.kat)
        if len(pk) != 2 * m * terms:
            raise SystemExit(f"public key length {len(pk)} != {2 * m * terms}")
        P = np.frombuffer(pk, dtype="<u2").astype(np.int64).reshape(m, terms) % q

    (args.output / "pk.bin").write_bytes(pk)
    (args.output / "message.bin").write_bytes(args.message.encode())
    helper = build_helper(args.output)

    started = time.perf_counter()
    target_text = subprocess.run(
        [str(helper), "h", str(args.output / "pk.bin"),
         str(args.output / "message.bin")],
        check=True, capture_output=True, text=True).stdout.split()
    h = np.array([int(v) for v in target_text], dtype=np.int64) % q
    report = {"source": args.source, "n": n, "m": m, "q": q, "message": args.message,
              "public_key_bytes": len(pk),
              "target_seconds": round(time.perf_counter() - started, 3),
              "target": [int(v) for v in h]}

    started = time.perf_counter()
    forged = forge_from_public(n, m, q, P, h, solve_seed=args.seed + 1,
                               timeout=args.timeout, flag_timeout=args.timeout,
                               slices=6, rescalings=48)
    report["forge_seconds"] = round(time.perf_counter() - started, 3)
    report["forge_status"] = forged["status"]
    report["forge_report"] = {k: v for k, v in forged.items() if k != "z"}
    if forged["status"] != "forged":
        report["acceptance"] = "forge failed"
        (args.output / "acceptance.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
        raise SystemExit(1)

    signature = encode_signature(forged["z"], q)
    (args.output / "signature.bin").write_bytes(signature)
    good = subprocess.run([str(helper), "verify", str(args.output / "pk.bin"),
                           str(args.output / "message.bin"),
                           str(args.output / "signature.bin")],
                          capture_output=True, text=True)
    flipped = bytearray(signature)
    flipped[0] ^= 1
    (args.output / "signature_flipped.bin").write_bytes(bytes(flipped))
    bad = subprocess.run([str(helper), "verify", str(args.output / "pk.bin"),
                          str(args.output / "message.bin"),
                          str(args.output / "signature_flipped.bin")],
                         capture_output=True, text=True)
    report["signature_bytes"] = len(signature)
    report["sig_verify_return_code"] = good.stdout.strip()
    report["negative_control_flipped_bit"] = bad.stdout.strip()
    report["acceptance"] = ("ACCEPTED by the authors' sig_verify"
                            if good.returncode == 0 else "REJECTED by sig_verify")
    (args.output / "acceptance.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "forge_report"}, indent=2))
    print(json.dumps(report["forge_report"], indent=2))
    raise SystemExit(0 if good.returncode == 0 else 1)


if __name__ == "__main__":
    main()
