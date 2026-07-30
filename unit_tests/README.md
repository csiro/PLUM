# PLUM Unit Tests

This directory contains unit tests for the PLUM metabolic gap-filling solver using Google Test framework.

## Test Organization

Tests are organized by class with each test file following the pattern `<class>_test.cpp`:

- `metabolite_test.cpp` - Tests for the Metabolite class
- `reaction_test.cpp` - Tests for the Reaction class
- `scenario_test.cpp` - Tests for the Scenario class
- `solution_test.cpp` - Tests for the Solution class
- `multisol_test.cpp` - Tests for the MultiSol class
- `params_test.cpp` - Tests for the Params class

## Building Tests

Tests are built automatically when you build PLUM. To disable test building:

```bash
cmake -DBUILD_TESTS=OFF ..
```

### Build Process

```bash
# From PLUM root directory
mkdir build && cd build
cmake ..
cmake --build .
```

The test executables will be created in `build/unit_tests/`.

## Running Tests

### Run All Tests

```bash
cd build/unit_tests
ctest
```

Or run the convenience script from the PLUM root:

```bash
./run_unit_tests.sh
```

### Run Individual Test Suites

```bash
cd build/unit_tests
./metabolite_test
./reaction_test
./scenario_test
./solution_test
./multisol_test
./params_test
```

### Run with Verbose Output

```bash
./metabolite_test --gtest_verbose
```

### Run Specific Tests

```bash
# Run specific test case
./metabolite_test --gtest_filter=MetaboliteTest.ConstructorInitializesCorrectly

# Run all tests in a fixture
./reaction_test --gtest_filter=ReactionTest.*

# Run tests matching a pattern
./reaction_test --gtest_filter=*Flux*
```

## Test Structure

Each test file follows this structure:

```cpp
#include <gtest/gtest.h>
#include "mosh/<class>.h"

using namespace mosh;

class <Class>Test : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize test objects
    }

    void TearDown() override {
        // Clean up test objects
    }

    // Test member variables
};

TEST_F(<Class>Test, TestName) {
    EXPECT_EQ(expected, actual);
    ASSERT_TRUE(condition);
}
```

## Test Coverage

### Metabolite Tests (metabolite_test.cpp)
- Constructor initialization
- Name setters/getters
- Index management
- Boolean flags (source, cycle_met, c_source, dummy)
- Output formatting (write_to)
- Edge cases (empty names, long names, special characters)

### Reaction Tests (reaction_test.cpp)
- Constructor initialization
- Name and index management
- Objective coefficient handling
- Active/selected/dummy flags
- Flux bounds and known flux
- Reduced cost
- Metabolite coefficient management (set_coeff, remove_met)
- Input/output metabolite lists
- Exchange (EX) and demand (DM) reaction detection
- Export reaction detection
- Formula generation and nicification
- Reaction comparison (same_as, reverse_of)
- Output formatting (write_to)

### Scenario Tests (scenario_test.cpp)
- Constructor initialization
- Adding metabolites and reactions
- Index assignment for metabolites and reactions
- Accessor methods
- Empty scenario handling

### Solution Tests (solution_test.cpp)
- Constructor initialization
- Flux setting and getting (by index and reaction pointer)
- uses_react predicate
- Fill operation
- Copy constructor
- Flux value independence
- Edge cases (large values, negative values)

### MultiSol Tests (multisol_test.cpp)
- Constructor initialization
- Solution count tracking
- Accessor methods

### Params Tests (params_test.cpp)
- Constructor initialization
- Object creation

## Google Test Features Used

- **TEST_F**: Test fixtures with SetUp/TearDown
- **EXPECT_EQ**: Non-fatal assertion for equality
- **EXPECT_NE**: Non-fatal assertion for inequality
- **EXPECT_TRUE/FALSE**: Boolean assertions
- **EXPECT_GT/LT**: Comparison assertions
- **ASSERT_***: Fatal assertions (stop test on failure)

## Adding New Tests

1. Create a new test file: `<class>_test.cpp`
2. Add test fixture class inheriting from `::testing::Test`
3. Implement `SetUp()` and `TearDown()` methods
4. Add test cases using `TEST_F(<Fixture>, <TestName>)`
5. Update `CMakeLists.txt` to include the new test executable
6. Rebuild and run tests

## Dependencies

- Google Test (fetched automatically via CMake FetchContent)
- PLUM library (plumlib)
- Lime library
- Solver libraries (same as PLUM)

## Continuous Integration

These tests can be integrated into CI/CD pipelines:

```bash
# Build
cmake -S . -B build
cmake --build build

# Test
cd build/unit_tests && ctest --output-on-failure
```

## Troubleshooting

### Tests Fail to Build

- Ensure solver environment variables are set (GUROBI_HOME, HIGHS_HOME, LP_SOLVE_HOME)
- Check that lime submodule is initialized: `git submodule update --init --recursive`

### Tests Fail to Run

- Verify solver libraries are in LD_LIBRARY_PATH
- Check that plumlib built successfully

### Individual Test Failures

- Run with `--gtest_verbose` for detailed output
- Use `--gtest_filter` to isolate failing tests
- Check test fixture SetUp/TearDown for proper initialization

## Future Enhancements

Potential areas for additional testing:

- Solver tests (IntSolver, IncrSolver, LnsMxSolver)
- PathFinder reachability tests
- Experiment tests
- File I/O tests (reading PLD format)
- Integration tests with toy problems
- Performance benchmarks
- Memory leak detection (valgrind integration)
