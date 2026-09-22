#!/bin/sh
# Build the exact small-prime linear-algebra backend against the local FLINT.
#
# FLINT is looked for in this order: $FACTO_FLINT_PREFIX, whatever pkg-config
# reports, the Homebrew prefix, the active conda prefix, then /usr and
# /usr/local.  Set CC to choose the compiler.  macOS and Linux are supported;
# anything else stops with a message rather than producing a library that
# cannot be loaded.
set -eu

REPO=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
TARGET="$REPO/src/factodsa"
CC=${CC:-cc}

PREFIX=
for candidate in "${FACTO_FLINT_PREFIX:-}" \
                 "$(pkg-config --variable=prefix flint 2>/dev/null)" \
                 "$(brew --prefix flint 2>/dev/null)" \
                 "${CONDA_PREFIX:-}" \
                 /usr /usr/local; do
  if [ -n "$candidate" ] && [ -f "$candidate/include/flint/nmod_mat.h" ]; then
    PREFIX="$candidate"
    break
  fi
done

if [ -z "$PREFIX" ]; then
  echo "FLINT headers not found.  Install FLINT, or set FACTO_FLINT_PREFIX to a" >&2
  echo "prefix that contains include/flint/nmod_mat.h and lib/libflint.*" >&2
  exit 1
fi
if ! command -v "$CC" >/dev/null 2>&1; then
  echo "C compiler '$CC' not found; set CC to one that exists." >&2
  exit 1
fi

case "$(uname -s)" in
  Darwin) LIB=libfinite_field.dylib
          "$CC" -dynamiclib -O3 -I"$PREFIX/include" -L"$PREFIX/lib" \
                "$TARGET/finite_field.c" -o "$TARGET/$LIB" -lflint ;;
  Linux)  LIB=libfinite_field.so
          "$CC" -shared -fPIC -O3 -I"$PREFIX/include" -L"$PREFIX/lib" \
                "$TARGET/finite_field.c" -o "$TARGET/$LIB" -lflint ;;
  *)      echo "unsupported platform '$(uname -s)': this script builds for macOS" >&2
          echo "and Linux only.  Elsewhere, compile finite_field.c into the" >&2
          echo "library name that src/factodsa/finite_field.py looks for." >&2
          exit 1 ;;
esac

echo "built $TARGET/$LIB (FLINT prefix $PREFIX, compiler $CC)"
