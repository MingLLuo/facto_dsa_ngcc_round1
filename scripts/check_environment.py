"""Verify that everything needed to reproduce the results is present.

Run this first.  It exits non-zero if a required component is missing, so a
failure here means the repository is not ready rather than the attack failing.
"""
from __future__ import annotations

import importlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

REQUIRED = ["numpy"]
OPTIONAL = ["matplotlib"]
REFERENCE = REPO / "vendor" / "reference-128"
KAT = REPO / "vendor" / "kat" / "KAT_SIG_Facto-DSA-128.txt"
SPEC = REPO / "vendor" / "spec" / "algorithm-specification-facto-dsa.pdf"
MIN_MSOLVE = (0, 10)


def report(label, ok, detail=""):
    print(f"[{'ok  ' if ok else 'FAIL'}] {label}{'  ' + detail if detail else ''}")
    return ok


def main():
    ok = True
    print("Facto-DSA reproduction: environment check")
    print(f"python {sys.version.split()[0]}  ({sys.executable})")
    ok &= report("python >= 3.10", sys.version_info >= (3, 10))

    for name in REQUIRED:
        try:
            module = importlib.import_module(name)
            ok &= report(f"python package {name}", True, getattr(module, "__version__", ""))
        except Exception as exc:  # pragma: no cover - diagnostic path
            ok &= report(f"python package {name}", False, type(exc).__name__)
    for name in OPTIONAL:
        try:
            module = importlib.import_module(name)
            report(f"optional package {name}", True, getattr(module, "__version__", ""))
        except Exception:
            report(f"optional package {name}", False, "figures will be skipped")

    compiler = shutil.which(os.environ.get("CC", "cc"))
    ok &= report("C compiler", compiler is not None, compiler or "")

    msolve = shutil.which("msolve")
    if msolve:
        try:
            text = subprocess.run([msolve, "-h"], capture_output=True, text=True,
                                  timeout=30).stdout
            version = ""
            for line in text.splitlines():
                if "version" in line:
                    version = line.strip().split()[-1]
                    break
            parsed = tuple(int(part) for part in version.split(".")[:2]) if version else (0, 0)
            ok &= report("msolve >= 0.10", parsed >= MIN_MSOLVE, version)
        except Exception as exc:
            ok &= report("msolve runs", False, type(exc).__name__)
    else:
        ok &= report("msolve on PATH", False, "install msolve >= 0.10")

    for name in ["SIG_AlgorithmInstance.c", "SIG_AlgorithmInstance.h", "auxfunc.c",
                 "auxfunc.h", "drng.c", "drng.h"]:
        ok &= report(f"vendored {name}", (REFERENCE / name).is_file())
    report("vendored KAT public keys", KAT.is_file(), str(KAT.name))
    report("vendored specification", SPEC.is_file(), str(SPEC.name))

    def usable_backend():
        """The library name that actually loaded and answered, or an exception."""
        from factodsa.finite_field import _lib, rank
        import numpy as np
        if rank(np.eye(3, dtype=np.int64), 65519) != 3:
            raise RuntimeError("backend returned a wrong rank")
        return Path(_lib._name).name

    try:
        name = usable_backend()
    except Exception:
        print("[info] FLINT backend missing or unusable; trying to build it")
        shell = shutil.which("sh")
        try:
            if shell is None:
                raise RuntimeError("no 'sh' to run build_backend.sh")
            subprocess.run([shell, str(REPO / "scripts" / "build_backend.sh")],
                           check=True)
            name = usable_backend()
        except Exception as exc:
            name = None
            ok &= report("FLINT backend", False, f"{type(exc).__name__}: {exc}")
    if name:
        ok &= report("FLINT backend", True, name)

    try:
        from factodsa import attack, rank1  # noqa: F401
        ok &= report("attack package imports", True)
    except Exception as exc:
        ok &= report("attack package imports", False, f"{type(exc).__name__}: {exc}")

    print()
    print("READY" if ok else "NOT READY - fix the FAIL lines above")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
