#!/bin/bash

SOURCE_DIR="src"
TEST_DIR="test"
BUILD_DIR="build"
EXECUTABLE="$BUILD_DIR/run_tests"

source build_utils.sh

COMPILER=$(get_compiler)
FLAGS=$(get_flags)
ALIGN=$(get_alignment)
MODULE_FLAG=$(get_module_flag $BUILD_DIR)

echo "Detected alignment: $ALIGN"

echo "Testing safeguard for C kind mismatches:"
bash <<'EOF'
failed=0
for flag in TEST_KIND_MISMATCH_C_INT TEST_KIND_MISMATCH_C_DOUBLE TEST_KIND_MISMATCH_C_DOUBLE_COMPLEX; do
  echo -en "$flag: "
  if bash build.sh "$@" -D${flag} 1>/dev/null 2>/dev/null; then
    echo "failure"
    failed=1
  else
    echo "success"
  fi
done
exit $failed
EOF
check_exit_code "Kind Mismatch Test failed"


echo "Compiling src/"
bash build.sh $@
check_exit_code "Build failed"

echo "Using compiler: $COMPILER"

handle_args $@

mkdir -p $BUILD_DIR

# Clean up any existing run_tests file/directory
if [ -e "$EXECUTABLE" ]; then
  echo "Removing existing $EXECUTABLE..."
  rm -rf "$EXECUTABLE"
fi

echo "Compiling test modules..."
# Then compile test/ modules using .mod files from build/
$COMPILER $FLAGS $MODULE_FLAG -DDEFAULT_ALIGNMENT=$ALIGN $MAX_PERF_FLAG \
  -I$BUILD_DIR \
  -c $TEST_DIR/*.[fF]90

# Move object files to build/
mv *.o $BUILD_DIR/ 2>/dev/null || true
mv *.mod $BUILD_DIR/ 2>/dev/null || true

check_build "*.o" "*.mod"
check_exit_code "Module compilation failed"


echo "Linking executable..."
# Finally link everything together
$COMPILER $FLAGS -I$BUILD_DIR \
  $BUILD_DIR/*.o -o $EXECUTABLE

check_exit_code "Executable compilation failed"

echo "Running tests..."
# Run the executable
$EXECUTABLE $ARGS

check_exit_code "Tests failed"

rm -f test_*.bin
rm -f test_*.zip