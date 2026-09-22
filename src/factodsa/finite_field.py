"""Exact small-prime matrix operations through the installed FLINT library."""
import ctypes
from pathlib import Path

try:
    import numpy as np
except ImportError as error:                      # pragma: no cover - setup path
    raise ImportError(
        "the attack needs numpy; install it (pip install numpy) or run inside "
        "an environment that has it"
    ) from error

_root = Path(__file__).resolve().parent
_SYMBOLS = ("ff_rank", "ff_rref", "ff_kernel", "ff_det", "ff_solve", "ff_mul")
_lib = None
for _candidate in ("libfinite_field.dylib", "libfinite_field.so", "finite_field.dylib"):
    _path = _root / _candidate
    if _path.exists():
        try:
            _library = ctypes.CDLL(str(_path))
        except OSError:
            continue
        if all(hasattr(_library, _symbol) for _symbol in _SYMBOLS):
            _lib = _library
            break
if _lib is None:
    raise OSError(
        "finite_field native backend missing or out of date (needs "
        f"{', '.join(_SYMBOLS)}); run build_backend.sh first"
    )
_ptr = ctypes.POINTER(ctypes.c_uint64)
_base = [_ptr, ctypes.c_long, ctypes.c_long, ctypes.c_ulong]
_lib.ff_rank.argtypes = _base
_lib.ff_rank.restype = ctypes.c_long
for _name in ("ff_kernel", "ff_rref"):
    getattr(_lib, _name).argtypes = _base + [_ptr]
    getattr(_lib, _name).restype = ctypes.c_long
_lib.ff_det.argtypes = [_ptr, ctypes.c_long, ctypes.c_ulong]
_lib.ff_det.restype = ctypes.c_uint64
_lib.ff_solve.argtypes = [_ptr, ctypes.c_long, ctypes.c_long, _ptr, ctypes.c_long, ctypes.c_ulong, _ptr]
_lib.ff_solve.restype = ctypes.c_int
_lib.ff_mul.argtypes = _lib.ff_solve.argtypes
_lib.ff_mul.restype = None

def _array(A, q):
    a = np.asarray(A)
    if a.ndim != 2 or a.dtype.kind not in "iu":
        raise ValueError("Expected a two-dimensional integer matrix")
    if q < 2 or q > 65519:
        raise ValueError("This exact-integer backend is bounded to q <= 65519")
    return np.ascontiguousarray(a % q, dtype=np.uint64)

def rank(A, q):
    a = _array(A, q)
    return int(_lib.ff_rank(a.ctypes.data_as(_ptr), *a.shape, q))

def kernel(A, q):
    """Return a row basis of the right nullspace."""
    a = _array(A, q)
    out = np.zeros((a.shape[1], a.shape[1]), dtype=np.uint64)
    k = _lib.ff_kernel(a.ctypes.data_as(_ptr), *a.shape, q, out.ctypes.data_as(_ptr))
    return out[:k].astype(np.int64)

def rref(A, q):
    """Return full reduced matrix and its pivot columns."""
    a = _array(A, q); out = np.zeros_like(a)
    k = _lib.ff_rref(a.ctypes.data_as(_ptr), *a.shape, q, out.ctypes.data_as(_ptr))
    pivots = [int(np.flatnonzero(out[i])[0]) for i in range(k)]
    return out.astype(np.int64), pivots

def solve(A, B, q):
    """One solution of A X=B, or raise for an inconsistent system."""
    a = _array(A, q)
    b = np.asarray(B, dtype=np.int64)
    vector = b.ndim == 1
    if vector: b = b[:, None]
    b = _array(b, q)
    if b.shape[0] != a.shape[0]: raise ValueError("Row mismatch")
    out = np.zeros((a.shape[1], b.shape[1]), dtype=np.uint64)
    ok = _lib.ff_solve(a.ctypes.data_as(_ptr), *a.shape, b.ctypes.data_as(_ptr),
                      b.shape[1], q, out.ctypes.data_as(_ptr))
    if not ok: raise ValueError("Inconsistent finite-field linear system")
    result = out.astype(np.int64)
    if not np.array_equal((a.astype(np.int64) @ result) % q, b):
        raise ArithmeticError("Exact solve certificate failed")
    return result[:, 0] if vector else result

def det(A, q):
    """Determinant of a square matrix over F_q."""
    a = _array(A, q)
    if a.shape[0] != a.shape[1]:
        raise ValueError("Expected a square matrix")
    return int(_lib.ff_det(a.ctypes.data_as(_ptr), a.shape[0], q)) % q

def inverse(A, q):
    a = _array(A, q)
    if a.shape[0] != a.shape[1] or rank(a, q) != len(a):
        raise ValueError("Expected an invertible square matrix")
    return solve(a, np.eye(len(a), dtype=np.int64), q)

def column_basis(A, q):
    R, pivots = rref(np.asarray(A).T, q)
    return R[:len(pivots)].T

def matmul(A, B, q):
    a = _array(A, q); b = _array(B, q)
    if a.shape[1] != b.shape[0]: raise ValueError("Matrix product shape mismatch")
    out = np.zeros((a.shape[0], b.shape[1]), dtype=np.uint64)
    _lib.ff_mul(a.ctypes.data_as(_ptr), *a.shape, b.ctypes.data_as(_ptr),
                b.shape[1], q, out.ctypes.data_as(_ptr))
    return out.astype(np.int64)
