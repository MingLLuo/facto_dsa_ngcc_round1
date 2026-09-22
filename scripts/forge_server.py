"""Local interactive forger: your message in, verified signatures out.

    python3 scripts/forge_server.py            # open the URL it prints
    python3 scripts/forge_server.py --count 3 --port 8765

The structural recovery and the triangular flag are target-independent, so they
run once at start-up (about 13 s for the vendored 128-bit KAT key).  After that
every request is cheap: the target comes from the submitters' own SM3/XOF code
and each candidate preimage is checked with their own ``sig_verify``, both
compiled from ``vendor/reference-128``.  Ask for several signatures and you get
several distinct preimages of the same target.

The server binds to 127.0.0.1 and talks to nothing outside this machine.
"""
from __future__ import annotations

import argparse
import json
import re
import struct
import subprocess
import sys
import tempfile
import time
from math import comb
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

from factodsa.acceptance import (DEFAULT_KAT, DEFAULT_REFERENCE, build_helper,  # noqa: E402
                                 encode_signature, first_kat_public_key)
from factodsa.attack import forge_target, recover_structure  # noqa: E402
from factodsa.polynomials import evaluate  # noqa: E402

UI = Path(__file__).resolve().parent / "forge_ui.html"
STATE = {}


def load_public_key(args):
    """Return (public key bytes, coefficient matrix), from KAT, file or model."""
    if args.source == "model":
        from factodsa.model import generate
        P, _ = generate(args.n, args.m, args.q, args.seed)
        P = np.asarray(P, dtype=np.int64) % args.q
        return b"".join(struct.pack("<H", int(v) % args.q) for v in P.reshape(-1)), P
    public_key = (args.pk_file.read_bytes() if args.source == "pkfile"
                  else first_kat_public_key(args.kat))
    P = np.frombuffer(public_key, dtype="<u2").astype(np.int64).reshape(args.m, -1) % args.q
    expected = comb(2 * args.n + 2, 3)
    if P.shape[1] != expected:
        raise SystemExit(f"public key has {P.shape[1]} coefficients per row, "
                         f"expected C(2n+2,3) = {expected} for n = {args.n}")
    return public_key, P


def parse_target(text, q, m):
    """A target is ``m`` decimal values, or ``2m`` bytes of hex (little-endian)."""
    cleaned = re.sub(r"[,\s]+", " ", text.strip()).strip()
    if not cleaned:
        raise ValueError("empty target")
    if re.fullmatch(r"\d+(?: \d+)*", cleaned):          # decimal list wins
        values = [int(v) for v in cleaned.split()]
        if len(values) != m:
            raise ValueError(f"expected {m} numbers, got {len(values)}")
        return [v % q for v in values]
    try:
        raw = bytes.fromhex(re.sub(r"[^0-9a-fA-F]", "", text))
    except ValueError as error:
        raise ValueError("target is neither a decimal list nor hex") from error
    if len(raw) != 2 * m:
        raise ValueError(f"expected {m} decimal numbers or {2 * m} bytes of hex")
    return list(struct.unpack(f"<{m}H", raw))
    raise ValueError(f"expected {m} numbers or {2 * m} bytes of hex")


def helper_target(message: str) -> list:
    """Target for a message, from the submitters' own sm3/xof code."""
    with tempfile.TemporaryDirectory(prefix="forge-") as tmp:
        path = Path(tmp) / "m.bin"
        path.write_bytes(message.encode())
        out = STATE["run"](["h", str(STATE["pkfile"]), str(path)])
        return [int(v) for v in out.split()]


def helper_verify(message: str, signature: bytes) -> int:
    """The submitters' sig_verify return code for this message and signature."""
    with tempfile.TemporaryDirectory(prefix="forge-") as tmp:
        message_path = Path(tmp) / "m.bin"
        signature_path = Path(tmp) / "s.bin"
        message_path.write_bytes(message.encode())
        signature_path.write_bytes(signature)
        return int(STATE["run"](["verify", str(STATE["pkfile"]), str(message_path),
                                 str(signature_path)]))


def forge(request):
    """Forge `count` distinct signatures for one message or target."""
    q, n, m = STATE["q"], STATE["n"], STATE["m"]
    structure = STATE["structure"]
    count = max(1, min(int(request.get("count", 1)), 20))
    mode = request.get("mode", "message")
    started = time.perf_counter()

    if mode == "target":
        message = None
        target = parse_target(str(request.get("target", "")), q, m)
    else:
        message = str(request.get("message", ""))
        target = helper_target(message)

    h = np.array(target, dtype=np.int64) % q
    results, seen = [], set()
    for attempt in range(count * 8):
        if len(results) >= count:
            break
        seed = int(request.get("seed", 0)) + attempt
        begun = time.perf_counter()
        forged = forge_target(structure, h, solve_seed=seed,
                              timeout=STATE["timeout"], slices=8)
        z = forged.get("z")
        if not z or tuple(z) in seen:
            continue
        seen.add(tuple(z))
        signature = encode_signature(z, q)
        checked = evaluate(STATE["P"], np.array(z, dtype=np.int64), q)[0] % q
        results.append(dict(
            seed=seed,
            z=[int(v) for v in z],
            hex=signature.hex(" "),
            seconds=round(time.perf_counter() - begun, 2),
            coefficients_match=bool(np.array_equal(checked, h)),
            sig_verify_return=(helper_verify(message, signature) if message is not None
                               else None),
        ))
    return dict(target=[int(v) for v in h], mode=mode, message=message,
                count=len(results), requested=count,
                seconds=round(time.perf_counter() - started, 2),
                signatures=results)


class Handler(BaseHTTPRequestHandler):
    def _send(self, body: bytes, status=200, kind="application/json"):
        self.send_response(status)
        self.send_header("Content-Type", kind)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path.startswith("/info"):
            self._send(json.dumps(STATE["info"]).encode())
        else:
            self._send(UI.read_bytes(), kind="text/html; charset=utf-8")

    def do_POST(self):
        length = int(self.headers.get("Content-Length", 0))
        try:
            request = json.loads(self.rfile.read(length) or b"{}")
            self._send(json.dumps(forge(request)).encode())
        except Exception as error:                     # surfaced in the page
            self._send(json.dumps({"error": f"{type(error).__name__}: {error}"}).encode(), status=400)

    def log_message(self, fmt, *args):
        print("  " + fmt % args, flush=True)


def prepare(args):
    """Load the key, build the helper, recover the structure once."""
    workdir = Path(tempfile.mkdtemp(prefix="forge-run-"))
    helper = build_helper(workdir)
    public_key, P = load_public_key(args)
    public_key_path = workdir / "pk.bin"
    public_key_path.write_bytes(public_key)
    STATE["pkfile"] = public_key_path

    def run(argv):
        # `sig_verify` exits 1 on a rejected signature, which is a result, not
        # an error; anything above that is a real failure.
        out = subprocess.run([str(helper), *argv], capture_output=True, text=True)
        if out.returncode > 1:
            raise RuntimeError(out.stderr.strip() or f"helper exited {out.returncode}")
        return out.stdout.strip()

    STATE.update(run=run, pk=public_key, P=P, q=args.q, n=args.n, m=args.m,
                 timeout=args.timeout)
    print(f"public key: {len(public_key)} bytes, {P.shape[0]} x {P.shape[1]} coefficients")
    print("recovering the structure once (target-independent)...", flush=True)
    started = time.perf_counter()
    structure = recover_structure(args.n, args.m, args.q, P,
                                  flag_timeout=args.flag_timeout)
    if structure["status"] != "recovered":
        raise SystemExit(f"recovery failed: {structure['status']}")
    STATE["structure"] = structure
    seconds = time.perf_counter() - started
    STATE["info"] = dict(source=args.source, key_bytes=len(public_key), n=args.n,
                         m=args.m, q=args.q, recovery_seconds=round(seconds, 1),
                         count_default=args.count)
    return seconds


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=["kat", "model", "pkfile"], default="kat")
    parser.add_argument("--pk-file", type=Path)
    parser.add_argument("--kat", type=Path, default=DEFAULT_KAT)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--n", type=int, default=10)
    parser.add_argument("--m", type=int, default=13)
    parser.add_argument("--q", type=int, default=65519)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--count", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--flag-timeout", type=int, default=1800)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    seconds = prepare(args)
    print(f"structure recovered in {seconds:.1f} s "
          f"(flag {STATE['structure']['flag_seconds']:.1f} s); serving on "
          f"http://{args.host}:{args.port}/", flush=True)
    ThreadingHTTPServer((args.host, args.port), Handler).serve_forever()


if __name__ == "__main__":
    main()
