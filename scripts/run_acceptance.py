"""Forge a signature and hand it to the authors' own verifier.

Wrapper around ``factodsa.acceptance`` so that every entry point is a script
under ``scripts/`` and no PYTHONPATH has to be set by hand.  Options are the
same, e.g. ``--source kat`` (default), ``--source model``, ``--source pkfile``.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from factodsa.acceptance import main  # noqa: E402

if __name__ == "__main__":
    main()
