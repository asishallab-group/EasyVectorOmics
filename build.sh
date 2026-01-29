#!/bin/bash
# build.sh | Optimized build script for FPM with dynamic alignment
# Build with selected profile and alignment parameter:
# Default fallback alignment for the most likely situation:

source build_utils.sh

COMPILER=$(get_compiler)
FLAGS=$(get_flags)
ALIGN=$(get_alignment)
handle_args "$@"

# Clean build directory if it exists
if [[ -d "build" && -z "$KEEP_OLD_BUILD_DIR" ]]; then
  rm -rf build
fi

# Ensure build directory exists before compiling .f files
mkdir -p build

# Define flags for Fortran 77 and Fortran 90
FFLAGS_F77="-O2 -std=legacy -ffixed-line-length-none -fallow-argument-mismatch -fPIC -w"
FFLAGS_F90="$FLAGS"

# Define source files
SRC_F="external/loess/loessf.f external/loess/supp.f external/loess/d1mach.f external/loess/dqrsl.f external/loess/dsvdc.f external/loess/daxpy.f external/loess/ddot.f external/loess/dnrm2.f external/loess/dscal.f external/loess/dswap.f external/loess/dcopy.f external/loess/drot.f"
SRC_F90="external/loess/drotg.f90"

# Define object files
OBJ_F="$(echo $SRC_F | sed 's/\.f/.o/g')"
OBJ_F90="$(echo $SRC_F90 | sed 's/\.f90/.o/g')"
OBJ="$OBJ_F $OBJ_F90"

# Compile Fortran 77 files
for f_file in $SRC_F; do
  echo "Compiling Fortran 77 file: $f_file"
  $COMPILER $FFLAGS_F77 -c $f_file -o build/$(basename "$f_file" .f).o || {
    echo "Failed to compile $f_file";
    exit 1;
  }
done

# Correct the path for Fortran 90 files in the loop
for f90_file in $SRC_F90; do
  echo "Compiling Fortran 90 file: $f90_file"
  $COMPILER $FFLAGS_F90 -c $f90_file -o build/$(basename "$f90_file" .f90).o || {
    echo "Failed to compile $f90_file";
    exit 1;
  }
done

# Create static library
ar rcs build/libloess_netlib.a build/*.o
ranlib build/libloess_netlib.a

echo "Static library created: build/libloess_netlib.a"

# Build with FPM after compiling .f files
fpm build --compiler $COMPILER --flag "$FLAGS $DIRECTIVES" --flag "-DDEFAULT_ALIGNMENT=$ALIGN" --flag "$MAX_PERF_FLAG" --flag -g

check_exit_code "Build with fpm failed"

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

# Link libloess_netlib.a with tensor-omics.so
if [ -f build/libloess_netlib.a ]; then
  echo "Linking tensor-omics.so with libloess_netlib.a"
  $COMPILER -shared -o build/libtensor-omics.so build/*.o build/libloess_netlib.a || {
    echo "Failed to link tensor-omics.so with libloess_netlib.a";
    exit 1;
  }
else
  echo "libloess_netlib.a not found. Ensure it is built before linking."
  exit 1
fi

echo "Build complete with compiler: $COMPILER, alignment: $ALIGN bytes"
