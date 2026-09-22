"""Small helpers shared by the runner scripts."""
from __future__ import annotations

import sys


def peak_mib():
    """Peak resident memory of this process in MiB, or ``None`` if unavailable.

    ``resource`` is POSIX-only, and the two platforms report different units
    (bytes on macOS, kibibytes on Linux).
    """
    try:
        import resource
    except ImportError:                                  # pragma: no cover
        return None
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return round(value / (1024 ** 2 if sys.platform == "darwin" else 1024), 1)
