#!/bin/bash

SOURCE_DIR="src"
TEST_DIR="test"
BUILD_DIR="build"
EXECUTABLE="$BUILD_DIR/run_tests"

# Detect alignment - FIX: Initialize with default value
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

echo "Detected alignment: $ALIGN"

# Detect compiler and flags
if [[ "$FC" == "ifx" || "$FC" == "ifort" ]]; then
  FLAGS="-O3 -qopenmp -xHost -align array64byte -qopt-zmm-usage=high -qopt-prefetch=3 -qopt-matmul -fPIC"
  MODULE_FLAG="-module $BUILD_DIR"
  COMPILER="ifx"
else
  FLAGS="-O3 -march=native -mtune=native -fopenmp -ffast-math -funroll-loops -ftree-vectorize -fassociative-math -fPIC"
  MODULE_FLAG="-J$BUILD_DIR"
  COMPILER="gfortran"
fi

LIBS="-lzip -lxxhash"
FLAGS="$FLAGS $LIBS"

echo "Using compiler: $COMPILER"

MAX_PERF_FLAG=""
REMOVE_ZIP="true"
for arg in "$@"; do
  if [[ "$arg" == "--max-performance" ]]; then
    MAX_PERF_FLAG="-DMAX_PERFORMANCE"
  fi
  if [[ "$arg" == "--keep-zip" ]]; then
    REMOVE_ZIP="false"
  fi
done

mkdir -p $BUILD_DIR

# Clean up any existing run_tests file/directory
if [ -e "$EXECUTABLE" ]; then
  echo "Removing existing $EXECUTABLE..."
  rm -rf "$EXECUTABLE"
fi

echo "Building main library first..."
bash build.sh 

echo "Compiling test modules..."
# Compile test modules using .mod files from build/
$COMPILER $FLAGS $MODULE_FLAG -DDEFAULT_ALIGNMENT=$ALIGN $MAX_PERF_FLAG \
  -I$BUILD_DIR -I$SOURCE_DIR -I$TEST_DIR \
  -c $TEST_DIR/*.f90

compilation_result=$?
echo "Test compilation exit code: $compilation_result"

# Move object files to build/
mv *.o $BUILD_DIR/ 2>/dev/null || true
mv *.mod $BUILD_DIR/ 2>/dev/null || true

echo "Linking executable..."
# Find all object files including the C object file if it exists
OBJECT_FILES=($BUILD_DIR/*.o)

# Link everything together with proper libraries
$COMPILER $FLAGS -I$BUILD_DIR \
  "${OBJECT_FILES[@]}" -o $EXECUTABLE

linking_result=$?
echo "Linking exit code: $linking_result"

if [ $linking_result -eq 0 ]; then
    echo "Running tests..."
    # Run the executable
    $EXECUTABLE "$@"
    TEST_RESULT=$?
    echo "Test execution completed with exit code: $TEST_RESULT"
else
    echo "Error: Failed to link test executable"
    exit 1
fi

# Cleanup
rm -f test_*.bin 2>/dev/null || true

if [[ $REMOVE_ZIP == "true" ]]; then
  rm -f test_archive_*_*.zip
fi 

exit $TEST_RESULT