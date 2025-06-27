#!/bin/bash
# run_tests.sh | Compile and run tests with selected compiler and flags

SOURCE_DIR="src"
TEST_DIR="test"
BUILD_DIR="build"
EXECUTABLE="$BUILD_DIR/run_tests"

# Detect alignment
ALIGN=32
if lscpu | grep -q amx; then
  ALIGN=128
elif lscpu | grep -q avx512; then
  ALIGN=64
elif lscpu | grep -q avx2; then
  ALIGN=32
elif lscpu | grep -q sse2; then
  ALIGN=16
fi

# Detect compiler and flags
if [[ "$FC" == "ifx" || "$FC" == "ifort" ]]; then
  FLAGS="-O3 -qopenmp -xHost -align array64byte -qopt-zmm-usage=high -qopt-prefetch=3 -qopt-matmul -fPIC"
  COMPILER="ifx"
else
  FLAGS="-O3 -march=native -mtune=native -fopenmp -ffast-math -funroll-loops -ftree-vectorize -fassociative-math -fPIC"
  COMPILER="gfortran"
fi

MAX_PERF_FLAG=""
for arg in "$@"; do
  if [[ "$arg" == "--max-performance" ]]; then
    MAX_PERF_FLAG="-DMAX_PERFORMANCE"
  fi
done

mkdir -p $BUILD_DIR

# Compile all source and test files
$COMPILER $FLAGS -DDEFAULT_ALIGNMENT=$ALIGN $MAX_PERF_FLAG -I$SOURCE_DIR -I$TEST_DIR $SOURCE_DIR/*.F90 $TEST_DIR/*.f90 -o $EXECUTABLE

if [ $? -ne 0 ]; then
  echo "Compilation failed."
  exit 1
fi

# Run the executable, passing any arguments (for test selection)
$EXECUTABLE "$@"