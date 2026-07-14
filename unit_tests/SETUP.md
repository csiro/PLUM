# PLUM Unit Tests Setup Guide

## Overview

PLUM now includes comprehensive unit tests using Google Test framework. Tests are organized by class in the `unit_tests/` directory.

## Prerequisites

Before building and running tests, ensure you have:

1. **Lime library**: PLUM depends on the `lime` library (CSIRO internal)
   ```bash
   git submodule update --init --recursive
   ```

2. **Solver libraries** (at least one required):
   - Gurobi: Set `GUROBI_HOME` environment variable
   - HiGHS: Set `HIGHS_HOME` environment variable
   - lp_solve: Set `LP_SOLVE_HOME` environment variable

3. **Build tools**:
   - CMake 3.14+
   - C++20 compatible compiler (g++ 10+)
   - Make or Ninja

## Quick Start

```bash
# 1. Navigate to PLUM directory
cd /path/to/PLUM

# 2. Initialize submodules (if not already done)
git submodule update --init --recursive

# 3. Set up solver environment (example for HiGHS)
export HIGHS_HOME=/path/to/highs

# 4. Create build directory and configure
mkdir build && cd build
cmake ..

# 5. Build PLUM and tests
cmake --build .

# 6. Run all tests
cd unit_tests
ctest

# Or use the convenience script from PLUM root
cd /path/to/PLUM
./run_unit_tests.sh
```

## Build Configuration

### Enable/Disable Tests

Tests are enabled by default. To disable:

```bash
cmake -DBUILD_TESTS=OFF ..
```

### Build Types

```bash
# Debug build (with symbols)
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Release build (optimized, default)
cmake -DCMAKE_BUILD_TYPE=Release ..

# Release with debug info
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
```

### Static Build

```bash
cmake -DSTATIC=ON ..
```

## Test Executables

After building, test executables are located in `build/unit_tests/`:

- `metabolite_test` - Metabolite class tests (15 test cases)
- `reaction_test` - Reaction class tests (33 test cases)
- `scenario_test` - Scenario class tests (13 test cases)
- `solution_test` - Solution class tests (15 test cases)
- `multisol_test` - MultiSol class tests (5 test cases)
- `params_test` - Params class tests (2 test cases)

**Total: 83 test cases**

## Running Tests

### All Tests

```bash
# Using ctest (recommended)
cd build/unit_tests
ctest

# With verbose output
ctest --verbose

# Show output only on failure
ctest --output-on-failure
```

### Individual Test Suites

```bash
cd build/unit_tests

# Run specific test executable
./metabolite_test
./reaction_test
./scenario_test
./solution_test
./multisol_test
./params_test
```

### Filtered Tests

Google Test supports filtering:

```bash
# Run specific test case
./metabolite_test --gtest_filter=MetaboliteTest.ConstructorInitializesCorrectly

# Run all tests in a fixture
./reaction_test --gtest_filter=ReactionTest.*

# Run tests matching pattern
./reaction_test --gtest_filter=*Flux*

# Exclude tests
./reaction_test --gtest_filter=-*ExchangeReaction*

# Multiple filters (OR)
./reaction_test --gtest_filter=*Flux*:*Coeff*
```

### Test Output Options

```bash
# Colorized output
./metabolite_test --gtest_color=yes

# Brief output
./metabolite_test --gtest_brief=1

# Repeat tests (for flakiness detection)
./metabolite_test --gtest_repeat=10

# Shuffle test order
./metabolite_test --gtest_shuffle

# Random seed for reproducible shuffle
./metabolite_test --gtest_shuffle --gtest_random_seed=12345
```

## Test Structure

### Test Organization

```
unit_tests/
├── CMakeLists.txt          # Test build configuration
├── README.md               # Test documentation
├── SETUP.md               # This file
├── metabolite_test.cpp     # Metabolite class tests
├── reaction_test.cpp       # Reaction class tests
├── scenario_test.cpp       # Scenario class tests
├── solution_test.cpp       # Solution class tests
├── multisol_test.cpp       # MultiSol class tests
└── params_test.cpp         # Params class tests
```

### Test Fixture Pattern

Each test file follows this structure:

```cpp
#include <gtest/gtest.h>
#include "mosh/<class>.h"

using namespace mosh;

// Test fixture class
class <Class>Test : public ::testing::Test {
protected:
    void SetUp() override {
        // Create test objects
        object = new <Class>(...);
    }

    void TearDown() override {
        // Clean up
        delete object;
    }

    <Class>* object;
};

// Test cases
TEST_F(<Class>Test, TestCaseName) {
    EXPECT_EQ(expected, actual);
    ASSERT_TRUE(condition);
}
```

## Troubleshooting

### CMake Configuration Fails

**Problem**: `add_subdirectory given source "lime" which is not an existing directory`

**Solution**: Initialize the lime submodule:
```bash
git submodule update --init --recursive
```

**Problem**: `No solver libraries found`

**Solution**: Set at least one solver environment variable:
```bash
export HIGHS_HOME=/path/to/highs
# or
export GUROBI_HOME=/path/to/gurobi
# or
export LP_SOLVE_HOME=/path/to/lpsolve
```

### Build Fails

**Problem**: Compiler errors about C++20 features

**Solution**: Use a modern compiler:
```bash
# Check compiler version
g++ --version  # Should be 10+

# Or specify compiler explicitly
cmake -DCMAKE_CXX_COMPILER=g++-11 ..
```

**Problem**: Linker errors about solver libraries

**Solution**: Ensure solver libraries are accessible:
```bash
# Add solver lib directory to library path
export LD_LIBRARY_PATH=/path/to/solver/lib:$LD_LIBRARY_PATH
```

### Tests Fail to Run

**Problem**: `error while loading shared libraries`

**Solution**: Add solver libraries to LD_LIBRARY_PATH:
```bash
export LD_LIBRARY_PATH=$HIGHS_HOME/build/lib:$LD_LIBRARY_PATH
# or
export LD_LIBRARY_PATH=$GUROBI_HOME/lib:$LD_LIBRARY_PATH
```

**Problem**: Test executable not found

**Solution**: Ensure tests were built:
```bash
cmake -DBUILD_TESTS=ON ..
cmake --build .
```

### Individual Test Failures

**Problem**: Test assertions fail

**Solution**: Run with verbose output to diagnose:
```bash
./test_name --gtest_verbose
./test_name --gtest_filter=FailingTest.*
```

**Problem**: Memory leaks detected

**Solution**: Run with valgrind:
```bash
valgrind --leak-check=full ./metabolite_test
```

## Continuous Integration

### GitHub Actions Example

```yaml
name: Unit Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
        with:
          submodules: recursive

      - name: Install Dependencies
        run: |
          sudo apt-get update
          sudo apt-get install -y cmake g++

      - name: Configure
        run: |
          mkdir build && cd build
          cmake -DBUILD_TESTS=ON ..

      - name: Build
        run: cmake --build build

      - name: Test
        run: |
          cd build/unit_tests
          ctest --output-on-failure
```

### Jenkins Pipeline Example

```groovy
pipeline {
    agent any
    stages {
        stage('Checkout') {
            steps {
                checkout scm
                sh 'git submodule update --init --recursive'
            }
        }
        stage('Build') {
            steps {
                sh '''
                    mkdir -p build && cd build
                    cmake -DBUILD_TESTS=ON ..
                    cmake --build .
                '''
            }
        }
        stage('Test') {
            steps {
                sh 'cd build/unit_tests && ctest --output-on-failure'
            }
        }
    }
}
```

## Adding New Tests

1. **Create test file**: `unit_tests/<class>_test.cpp`

2. **Write test fixture**:
```cpp
class NewClassTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize
    }
    void TearDown() override {
        // Cleanup
    }
};
```

3. **Add test cases**:
```cpp
TEST_F(NewClassTest, TestCaseName) {
    EXPECT_EQ(expected, actual);
}
```

4. **Update CMakeLists.txt**:
```cmake
add_executable(newclass_test newclass_test.cpp)
target_link_libraries(newclass_test GTest::gtest_main plumlib)
gtest_discover_tests(newclass_test)
```

5. **Build and run**:
```bash
cmake --build build
cd build/unit_tests && ./newclass_test
```

## Best Practices

### Test Naming

- Test files: `<class>_test.cpp`
- Test fixtures: `<Class>Test`
- Test cases: Descriptive action/state names (e.g., `ConstructorInitializesCorrectly`)

### Assertions

- **EXPECT_**: Non-fatal, test continues
- **ASSERT_**: Fatal, test stops
- Use ASSERT_ when subsequent code depends on the condition

### Test Independence

- Each test should be independent
- Use SetUp/TearDown for initialization/cleanup
- Don't rely on test execution order

### Test Coverage

- Test normal cases
- Test boundary conditions
- Test error conditions
- Test edge cases (empty, null, large values)

## Performance

### Build Time

Approximate build times (Release mode):
- PLUM library: ~30-60 seconds
- Each test executable: ~10-20 seconds
- Total (with 6 test suites): ~3-5 minutes

### Test Execution Time

- Individual test suite: < 1 second
- All tests (83 cases): < 5 seconds

## Future Enhancements

Potential additions:

1. **Solver Tests**: Test LP/MILP solver implementations
2. **PathFinder Tests**: Test reachability algorithms
3. **Integration Tests**: Test with toy problems
4. **Performance Tests**: Benchmark solver performance
5. **Memory Tests**: Valgrind integration
6. **Coverage Reports**: gcov/lcov integration
7. **Fuzz Testing**: Random input generation
8. **Property-Based Tests**: RapidCheck integration (for future consideration)

## References

- [Google Test Documentation](https://google.github.io/googletest/)
- [CMake Testing](https://cmake.org/cmake/help/latest/manual/ctest.1.html)
- PLUM Documentation: `CLAUDE.md`, `fileformat.md`, `doc/overview/overview.md`
- UML Diagrams: `doc/uml/`
