#!/bin/bash

# === build.sh ===
# Custom build script for tensor-omics project

set -e

# Load .env configuration if present
if [ -f .env ]; then
  source .env
else
  echo "No .env file found, using default settings."
fi

# Set Fortran compiler
FC=${FC:-gfortran}

# Determine alignment (optional, customize as needed)
ALIGN=${DEFAULT_ALIGNMENT:-32}
if lscpu | grep -q amx; then ALIGN=128; fi

# Create build directory if it doesn't exist
mkdir -p build

# === Build shared library for R ===
echo "Building shared library for R..."
SRC="src/tox_sorting.F90 src/tox_normalization.F90 src/tox_normalization_wrappers.F90"
OUT="build/tox_normalization.so"

$FC -shared -fPIC $SRC -o $OUT

if [ -f "$OUT" ]; then
  echo "Shared library created at $OUT"
else
  echo "ERROR: Failed to create shared library"
  exit 1
fi

# === Build Python module with f2py ===
echo "Building Python module with f2py..."
F2PY_MOD="tox_normalization_py"
F2PY_SRC="src/tox_normalization.F90 src/tox_normalization_wrappers.F90"

python3 -m numpy.f2py -c -m $F2PY_MOD $F2PY_SRC

if ls ${F2PY_MOD}*.so 1> /dev/null 2>&1; then
  mv ${F2PY_MOD}*.so build/
  echo "Python module created at build/${F2PY_MOD}*.so"
else
  echo "ERROR: Failed to create Python module"
  exit 1
fi

echo "Build finished successfully."