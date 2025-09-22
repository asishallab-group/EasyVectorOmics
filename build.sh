#!/bin/bash
# build.sh | Optimized build script for FPM with dynamic alignment and xxHash support
# Build with selected profile and alignment parameter:
# Default fallback alignment for the most likely situation:
ALIGN=32
# Detect capabilities in order of descending priority:
if lscpu | grep -q amx; then
ALIGN=128
elif lscpu | grep -q avx512; then
ALIGN=64
elif lscpu | grep -q avx2; then
ALIGN=32
elif lscpu | grep -q sse2; then
ALIGN=16
fi

# Check if xxHash library is available
if [ -f /usr/lib/libxxhash.so ] || [ -f /usr/lib64/libxxhash.so ]; then
    HAS_XXHASH=1
    # Arch Linux typically places headers in /usr/include and libraries in /usr/lib
    XXHASH_FLAGS="-I/usr/include"
    XXHASH_LIBS="-L/usr/lib -lxxhash"
else
    echo "Warning: xxHash library not found. Trying to compile without it..."
    HAS_XXHASH=0
    XXHASH_FLAGS=""
    XXHASH_LIBS=""
fi
LIBS="-lzip"
# Detect compiler and choose appropriate profile:
if [[ "$FC" == "ifx" || "$FC" == "ifort" ]]; then
  FLAGS="-O3 -qopenmp -xHost -align array64byte -qopt-zmm-usage=high -qopt-prefetch=3 -qopt-matmul -fPIC $XXHASH_FLAGS"
  COMPILER="ifx"
  C_COMPILER="icc"
else
  FLAGS="-O3 -march=native -mtune=native -fopenmp -ffast-math -funroll-loops -ftree-vectorize -fassociative-math -fPIC $XXHASH_FLAGS"
  COMPILER="gfortran"
  C_COMPILER="gcc"
fi

# Detect --max-performance flag
MAX_PERF_FLAG=""
for arg in "$@"; do
  if [[ "$arg" == "--max-performance" ]]; then
    MAX_PERF_FLAG="-DMAX_PERFORMANCE"
  fi
done

# Clean build directory if it exists
if [ -d "build" ]; then
  rm -rf build
fi

# Create build directory
mkdir -p build

# Build with FPM
export FC
fpm build --compiler $COMPILER --flag "$FLAGS" --flag "-DDEFAULT_ALIGNMENT=$ALIGN" --flag "$MAX_PERF_FLAG" --flag $LIBS

# Move .mod, .o and .so files from FPM build directories to root
for compiler_dir in build/${COMPILER}_*; do
  if [ -d "$compiler_dir" ]; then
    echo "Processing FPM directory: $compiler_dir"
    # Move .so files if they exist
    mv "$compiler_dir"/*.so build/ 2>/dev/null || true
    # Move .mod files if they exist
    find "$compiler_dir" -name "*.mod" -exec mv {} build/ \; 2>/dev/null || true
    # Move .o files from subdirectories if they exist
    find "$compiler_dir" -name "*.o" -exec mv {} build/ \; 2>/dev/null || true
    # Remove the compiler-specific directory
    rm -rf "$compiler_dir"
  fi
done

echo "Build complete with compiler: $COMPILER, alignment: $ALIGN bytes"
if [ $HAS_XXHASH -eq 1 ]; then
    echo "xxHash support: Enabled"
else
    echo "xxHash support: Disabled (using fallback hashing)"
fi

# Verify that we have the necessary files for test_runner.sh
mod_count=$(find build -name "*.mod" | wc -l)
obj_count=$(find build -name "*.o" | wc -l) 
so_count=$(find build -name "*.so" | wc -l)

if [ $mod_count -eq 0 ] || [ $obj_count -eq 0 ]; then
  echo "Warning: Missing .mod or .o files needed for test compilation"
fi