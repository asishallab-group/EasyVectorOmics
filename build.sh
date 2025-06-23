#!/bin/bash

# === build.sh ===
# Optimized build script for tensor-omics project according to F42 guidelines (section 7)

set -e # If something fails, exit immediately

# Load .env configuration if present
if [ -f .env ]; then
  source .env
else
  echo "No .env file found, using default settings."
fi

# === Detect alignment based on CPU capabilities ===
ALIGN=32  # Fallback default
if lscpu | grep -q amx; then
  ALIGN=128
elif lscpu | grep -q avx512; then
  ALIGN=64
elif lscpu | grep -q avx2; then
  ALIGN=32
elif lscpu | grep -q sse2; then
  ALIGN=16
fi
export DEFAULT_ALIGNMENT=$ALIGN

# === Detect compiler and set appropriate flags ===
FC=${FC:-gfortran}

if ! command -v $FC &>/dev/null; then
  echo "ERROR: Compiler $FC not found"
  exit 1
fi



echo "Detected Fortran compiler: $FC"
echo "Using alignment: ${DEFAULT_ALIGNMENT}"

if [[ "$FC" == *ifx* || "$FC" == *ifort* ]]; then
  FPM_FLAGS=${FPM_FLAGS_IFX}
else
  FPM_FLAGS=${FPM_FLAGS_GFORTRAN}
fi

echo "Using FPM_FLAGS: $FPM_FLAGS"

# === Set output directory ===
mkdir -p build

# === Build shared library for R ===
echo "Building shared library for R..."
SRC="src/tox_sorting.F90 src/tox_normalization.F90 src/tox_normalization_wrappers.F90"
OUT="build/tox_normalization.so"

$FC -shared -fPIC $SRC -o $OUT $FPM_FLAGS -DDEFAULT_ALIGNMENT=$DEFAULT_ALIGNMENT

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

python3 -m numpy.f2py -c -m $F2PY_MOD $F2PY_SRC --fcompiler=$FC \
  --opt="$FPM_FLAGS -DDEFAULT_ALIGNMENT=$DEFAULT_ALIGNMENT"

if ls ${F2PY_MOD}*.so 1> /dev/null 2>&1; then
  mv ${F2PY_MOD}*.so build/
  echo "Python module created at build/${F2PY_MOD}*.so"
else
  echo "ERROR: Failed to create Python module"
  exit 1
fi

echo "Build complete with compiler: $FC, alignment: $ALIGN bytes, optimization flags: $FPM_FLAGS"

echo "Build finished successfully."
