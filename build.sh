#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

if [[ -z "${GENIE:-}" ]]; then
  echo "error: GENIE is not set" >&2
  echo "       Set GENIE to the root of an already-built GENIE installation." >&2
  exit 2
fi
export GENIE

if [[ ! -x "$GENIE/bin/genie-config" ]]; then
  echo "error: could not find executable \$GENIE/bin/genie-config" >&2
  echo "       Set GENIE to the root of an already-built GENIE installation." >&2
  exit 2
fi

if ! command -v root-config >/dev/null 2>&1; then
  echo "error: root-config was not found in PATH" >&2
  exit 2
fi

export PATH="$GENIE/bin:$PATH"
export GXMLPATH="${GXMLPATH:-$GENIE/config}"
export DYLD_LIBRARY_PATH="$GENIE/lib:${DYLD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$GENIE/lib:${LD_LIBRARY_PATH:-}"

CXX="${CXX:-c++}"
SRC="$SCRIPT_DIR/xsec_scan.cxx"
OUT_DIR="$SCRIPT_DIR/out"
EXE="$OUT_DIR/xsec_scan"
mkdir -p "$OUT_DIR"

ROOT_CFLAGS="$(root-config --cflags)"
ROOT_LIBS="$(root-config --glibs) -lMinuit -lGeom -lEG -lGenVector -lMathMore"
ROOT_LIBDIR="$(root-config --libdir 2>/dev/null || true)"

GENIE_LIBS="$("$GENIE/bin/genie-config" --libs | sed 's/-lGPhHadTransp//g')"

EXT_CFLAGS=""
EXT_LIBS=""

if command -v lhapdf-config >/dev/null 2>&1; then
  EXT_CFLAGS+=" $(lhapdf-config --cflags 2>/dev/null || true)"
  EXT_LIBS+=" $(lhapdf-config --libs 2>/dev/null || true)"
fi

if command -v xml2-config >/dev/null 2>&1; then
  EXT_CFLAGS+=" $(xml2-config --cflags 2>/dev/null || true)"
  EXT_LIBS+=" $(xml2-config --libs 2>/dev/null || true)"
fi

if command -v gsl-config >/dev/null 2>&1; then
  EXT_CFLAGS+=" $(gsl-config --cflags 2>/dev/null || true)"
  EXT_LIBS+=" $(gsl-config --libs 2>/dev/null || true)"
fi

if command -v log4cpp-config >/dev/null 2>&1; then
  EXT_CFLAGS+=" $(log4cpp-config --cflags 2>/dev/null || true)"
  EXT_LIBS+=" $(log4cpp-config --libs 2>/dev/null || true)"
elif command -v pkg-config >/dev/null 2>&1 && pkg-config --exists log4cpp; then
  EXT_CFLAGS+=" $(pkg-config --cflags log4cpp)"
  EXT_LIBS+=" $(pkg-config --libs log4cpp)"
fi

PYTHIA6_LIBDIR="${PYTHIA6_LIBDIR:-}"
if [[ -z "$PYTHIA6_LIBDIR" && -n "${PYTHIA6:-}" ]]; then
  if [[ -d "$PYTHIA6/lib" ]]; then
    PYTHIA6_LIBDIR="$PYTHIA6/lib"
  elif [[ -d "$PYTHIA6" ]]; then
    PYTHIA6_LIBDIR="$PYTHIA6"
  fi
fi
if [[ -n "$PYTHIA6_LIBDIR" && -d "$PYTHIA6_LIBDIR" ]]; then
  EXT_LIBS+=" -L$PYTHIA6_LIBDIR"
fi

RPATH_FLAGS="-Wl,-rpath,$GENIE/lib"
if [[ -n "$ROOT_LIBDIR" ]]; then
  RPATH_FLAGS+=" -Wl,-rpath,$ROOT_LIBDIR"
fi
if [[ -n "$PYTHIA6_LIBDIR" && -d "$PYTHIA6_LIBDIR" ]]; then
  RPATH_FLAGS+=" -Wl,-rpath,$PYTHIA6_LIBDIR"
fi

"$CXX" -std=c++17 -O2 -I"$GENIE/src" $ROOT_CFLAGS $EXT_CFLAGS ${GENIE_EXTRA_CFLAGS:-} \
  "$SRC" -o "$EXE" \
  $ROOT_LIBS $EXT_LIBS $GENIE_LIBS ${GENIE_EXTRA_LDFLAGS:-} $RPATH_FLAGS

echo "$EXE"
