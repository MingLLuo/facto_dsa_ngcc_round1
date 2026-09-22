"""Thin msolve interface: format a system, run it, parse the parametrisation.

Only used to solve small structured systems (rank-one point finding and the
recovered bilinear system).  Every parsed point is verified against the input
equations before it is returned, so a wrong reading of msolve's sign or
ordering convention shows up as "no points" instead of silent garbage.
"""
from __future__ import annotations

import ast
import itertools
import subprocess
import tempfile
from pathlib import Path

from .fieldtools import roots


def format_polynomial(terms):
    """terms: iterable of (coefficient, monomial) with monomial None for a constant."""
    pieces = []
    for coefficient, monomial in terms:
        if coefficient == 0:
            continue
        negative = coefficient < 0
        magnitude = abs(coefficient)
        if monomial is None:
            body = str(magnitude)
        else:
            body = monomial if magnitude == 1 else f"{magnitude}*{monomial}"
        if not pieces:
            pieces.append(("-" if negative else "") + body)
        else:
            pieces.append((" - " if negative else " + ") + body)
    return "".join(pieces) if pieces else "0"


def monomial_string(indices, prefix, start=0):
    """msolve's parser mishandles repeated factors (x*x), so group exponents."""
    counts = {}
    order = []
    for index in indices:
        index = int(index) + start
        if index not in counts:
            counts[index] = 0
            order.append(index)
        counts[index] += 1
    parts = []
    for index in order:
        power = counts[index]
        parts.append(f"{prefix}{index}" if power == 1 else f"{prefix}{index}^{power}")
    return "*".join(parts)


def run_msolve(variables, q, polynomials, timeout=120, threads=8, extra_args=()):
    """Run msolve on the given system and return its raw output text."""
    text = ",".join(variables) + "\n" + str(q) + "\n" + ",\n".join(polynomials) + "\n"
    with tempfile.TemporaryDirectory(prefix="msolve-run-") as tmp:
        inp = Path(tmp) / "system.ms"
        out = Path(tmp) / "out.txt"
        inp.write_text(text)
        try:
            subprocess.run(["msolve", "-f", str(inp), "-o", str(out), "-t", str(threads),
                            *extra_args],
                           capture_output=True, text=True, timeout=timeout, check=False)
        except subprocess.TimeoutExpired:
            return None
        except FileNotFoundError as error:
            raise RuntimeError(
                "msolve not found on PATH; install msolve >= 0.10 and make sure "
                "the 'msolve' binary is reachable"
            ) from error
        return out.read_text() if out.exists() else None


def _split_top_level(text):
    """Split msolve output into top-level bracketed blocks."""
    depth = 0
    start = None
    blocks = []
    for index, char in enumerate(text):
        if char == "[":
            if depth == 0:
                start = index
            depth += 1
        elif char == "]":
            depth -= 1
            if depth == 0 and start is not None:
                blocks.append(text[start:index + 1])
                start = None
    return blocks


def parse_msolve_points(output, variables, q, evaluate, report=None):
    """Parse msolve output into explicit F_q points.

    ``evaluate`` maps a candidate point (tuple) to a list of residuals; only
    points with all-zero residuals are kept, so a wrong reading of msolve's
    conventions shows up as "no points" rather than as silent garbage.
    ``report``, when given, is updated with ``decoded=True`` as soon as a block
    about this variable list is found, so a caller can tell "msolve answered,
    the slice has no F_q point" from "msolve declined this chart".
    """
    if not output:
        return []
    blocks = _split_top_level(output)
    points = []
    for block in blocks:
        try:
            data = ast.literal_eval(block)
        except (ValueError, SyntaxError):
            continue
        if not isinstance(data, list) or len(data) < 2:
            continue
        live = data[0]
        body = data[1]
        if not isinstance(body, list) or len(body) < 6:
            continue
        try:
            nvars = int(body[1])
            variables_in = body[3]
            linear_form = body[4]
            parametrisation = body[5]
        except (TypeError, IndexError, ValueError):
            continue
        if nvars != len(variables) or list(variables_in) != list(variables):
            continue
        if report is not None:
            report["decoded"] = True
        if isinstance(live, int) and live != 0:
            continue  # msolve reports -1 (no solutions) or 1 (positive dimensional)
        # msolve wraps the parametrisation as [ <flag>, [elim, den, coords] ]
        if (isinstance(parametrisation, list) and len(parametrisation) == 2
                and isinstance(parametrisation[1], list)):
            parametrisation = parametrisation[1]
        if not isinstance(parametrisation, list) or len(parametrisation) < 3:
            continue
        nonzero = [i for i, c in enumerate(linear_form) if int(c) % q]
        if len(nonzero) != 1:
            continue  # only the single-parameter chart is decoded
        slot = nonzero[0]
        elim_degree, elim_coefficients = parametrisation[0]
        _, den_coefficients = parametrisation[1]
        coordinate_blocks = parametrisation[2]
        others = [i for i in range(nvars) if i != slot]
        if len(coordinate_blocks) != len(others):
            continue
        for root in roots(list(elim_coefficients), q):
            root = root * pow(int(linear_form[slot]) % q, -1, q) % q
            denominator = 0
            for coefficient in reversed(list(den_coefficients)):
                denominator = (denominator * root + coefficient) % q
            if denominator == 0:
                continue
            inverse = pow(denominator, -1, q)
            values = []
            for entry in coordinate_blocks:
                # msolve nests each coordinate as [[deg, [coeffs]]]
                while isinstance(entry, list) and len(entry) == 1 and isinstance(entry[0], list):
                    entry = entry[0]
                if not isinstance(entry, list) or len(entry) != 2:
                    values = None
                    break
                _, coefficients = entry
                acc = 0
                for coefficient in reversed(list(coefficients)):
                    acc = (acc * root + coefficient) % q
                values.append(acc * inverse % q)
            if values is None:
                continue
            # msolve fixes neither the coordinate order nor the per-coordinate
            # sign for us, so enumerate both orders and every sign pattern and
            # keep exactly the points the input system verifies.
            orders = [values, list(reversed(values))]
            width = len(others)
            # Above a dozen coordinates the sign patterns are not enumerated
            # (2^width would explode); the all-plus and all-minus conventions
            # still cover what msolve has produced in practice.
            patterns = ([(1,) * width, (-1,) * width] if width > 12
                        else list(itertools.product((1, -1), repeat=width)))
            for ordered in orders:
                for pattern in patterns:
                    point = [0] * nvars
                    point[slot] = root
                    for index, value, sign in zip(others, ordered, pattern):
                        point[index] = sign * value % q
                    if all(r % q == 0 for r in evaluate(tuple(point))):
                        points.append(tuple(point))
    unique = []
    for point in points:
        if point not in unique:
            unique.append(point)
    return unique
