#!/bin/bash
# Run PLUM unit tests

set -e  # Exit on error

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BUILD_DIR="${SCRIPT_DIR}/build"

# Check if build directory exists
if [ ! -d "${BUILD_DIR}" ]; then
    echo "Error: Build directory not found at ${BUILD_DIR}"
    echo "Please run 'mkdir build && cd build && cmake .. && cmake --build .' first"
    exit 1
fi

# Check if unit_tests directory exists
if [ ! -d "${BUILD_DIR}/unit_tests" ]; then
    echo "Error: Unit tests not found in ${BUILD_DIR}/unit_tests"
    echo "Please rebuild with unit tests enabled (BUILD_TESTS=ON)"
    exit 1
fi

echo "Running PLUM unit tests..."
echo "=========================="
echo

cd "${BUILD_DIR}/unit_tests"

# Run all tests with ctest
if command -v ctest &> /dev/null; then
    echo "Running all tests with ctest..."
    ctest --output-on-failure
else
    echo "ctest not found, running tests individually..."

    # Find all test executables
    for test_exe in *_test; do
        if [ -x "$test_exe" ]; then
            echo ""
            echo "Running $test_exe..."
            echo "-------------------"
            ./"$test_exe"
        fi
    done
fi

echo ""
echo "=========================="
echo "All tests completed!"
