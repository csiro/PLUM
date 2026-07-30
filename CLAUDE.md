# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PLUM is a metabolic gap-filling solver for Flux Balance Analysis (FBA). It finds which reactions from a database should be included in a metabolic network to balance known inputs with outputs, or maximize biomass production. This is a computational biology tool that uses Linear Programming (LP) and Mixed Integer Linear Programming (MILP) to solve metabolic pathway problems.

The key innovation: uses LP-based flux minimization instead of traditional MILP approaches, achieving polynomial-time solutions (seconds vs. hours/weeks) while incorporating compartment modeling (External, Periplasm, Cytosol, Output) and reachability constraints.

## Architecture

### Core Components

**Solver Hierarchy**
- `LpSolver` (base): Abstract LP/MILP solver interface
- `LpSolverImp`: Implementation interface for specific solvers
  - `GrbLpSolverImp`: Gurobi implementation (requires GUROBI_HOME)
  - `HighsLpSolverImp`: HiGHS implementation (requires HIGHS_HOME)
  - `LpsLpSolverImp`: lp_solve implementation (requires LP_SOLVE_HOME)
- `IntSolver`: Integer programming solver
- `IncrSolver`: Incremental solver
- `LnsMxSolver`: Large Neighborhood Search mixed solver (main algorithmic solver)

**Data Model**
- `Scenario`: Top-level container holding metabolites, reactions, experiments
- `Metabolite`: Compounds with supply/residual rates per compartment
- `Reaction`: Stoichiometric transformations with objective coefficients and costs
- `Solution`: Flux assignments and reaction selections
- `MultiSol`: Collection of multiple solutions

**Key Algorithms**
- `PathFinder`: Reachability analysis (graph-based producibility)
- `Params`: Central parameter/configuration holder

### Executables

- `plum`: Main gap-filling solver
- `plummx`: Mixed solver variant
- `plumsp`: Solution processor/printer
- `plumchk`: Solution checker/validator
- `plumcmp`: Solution comparator
- `plummerge`: Merge multiple solutions

## Build System

PLUM uses CMake with C++20.

### Initial Setup

```bash
# Initialize submodules (lime library dependency)
git submodule update --init --recursive

# Load solver environment (if using Gurobi)
source module.sh
```

### Build Commands

```bash
# Configure (out-of-source build recommended)
mkdir build && cd build
cmake ..

# Or configure with specific solver environments
GUROBI_HOME=/path/to/gurobi HIGHS_HOME=/path/to/highs cmake ..

# Build (defaults to Release)
cmake --build .

# Build specific target
cmake --build . --target plum

# Build with specific configuration (if needed)
cmake --build . --config Debug

# Install
cmake --install . --prefix /path/to/install
```

### Build Types

- Release (default): Optimized build
- Debug: Debug symbols, no optimization
- RelWithDebInfo: Optimized with debug info
- MinSizeRel: Size-optimized

### Static Build

```bash
cmake -DSTATIC=ON ..
```

### Solver Configuration

PLUM supports three LP/MILP solvers (all three are compiled in via preprocessor flags):
- **Gurobi**: Set `GUROBI_HOME` environment variable
- **HiGHS**: Set `HIGHS_HOME` environment variable  
- **lp_solve**: Set `LP_SOLVE_HOME` environment variable

At least one solver must be configured for linking to succeed.

## Testing

PLUM has two types of tests: unit tests (Google Test) and reference tests (integration tests).

### Unit Tests (Google Test)

**Overview**: Comprehensive unit tests for core data model classes using the Google Test framework. Tests are organized by class with each test file following the pattern `<class>_test.cpp`.

**Created**: July 14, 2026  
**Framework**: Google Test v1.14.0 (fetched automatically via CMake FetchContent)  
**Language**: C++20  
**Location**: `unit_tests/`

#### Test Suites (83 test cases total)

| Test File | Class | Tests | Coverage |
|-----------|-------|-------|----------|
| `metabolite_test.cpp` | Metabolite | 15 | Construction, name/index management, boolean flags (source, cycle_met, c_source, dummy), I/O, edge cases |
| `reaction_test.cpp` | Reaction | 33 | Construction, stoichiometry, flux bounds, active/selected/dummy flags, metabolite coefficients, formulas, comparison (same_as, reverse_of), EX/DM/export detection, I/O |
| `scenario_test.cpp` | Scenario | 13 | Construction, adding metabolites/reactions, automatic index assignment, accessors, empty scenario handling |
| `solution_test.cpp` | Solution | 15 | Construction, flux get/set (by index and reaction), uses_react predicate, fill operation, copy constructor, flux independence |
| `multisol_test.cpp` | MultiSol | 5 | Construction, solution count, accessors |
| `params_test.cpp` | Params | 2 | Construction, object creation |

#### Test Structure

Each test file uses the Google Test fixture pattern:

```cpp
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

#### Building Tests

Tests are built automatically with PLUM. To disable:

```bash
cmake -DBUILD_TESTS=OFF ..
```

Build process:

```bash
# Configure and build (tests enabled by default)
mkdir build && cd build
cmake ..
cmake --build .
```

Test executables are created in `build/unit_tests/`:
- `metabolite_test`
- `reaction_test`
- `scenario_test`
- `solution_test`
- `multisol_test`
- `params_test`

#### Running Tests

**Run all tests** (using ctest):
```bash
cd build/unit_tests
ctest

# With verbose output
ctest --verbose

# Show output only on failure
ctest --output-on-failure
```

**Run all tests** (using convenience script):
```bash
# From PLUM root directory
./run_unit_tests.sh
```

**Run individual test suites**:
```bash
cd build/unit_tests
./metabolite_test
./reaction_test
./scenario_test
./solution_test
./multisol_test
./params_test
```

**Run filtered tests** (Google Test filter syntax):
```bash
# Run specific test case
./reaction_test --gtest_filter=ReactionTest.SetCoefficient

# Run all tests in a fixture
./reaction_test --gtest_filter=ReactionTest.*

# Run tests matching pattern
./reaction_test --gtest_filter=*Flux*

# Exclude tests
./reaction_test --gtest_filter=-*ExchangeReaction*

# Multiple filters (OR)
./reaction_test --gtest_filter=*Flux*:*Coeff*

# Colorized output
./metabolite_test --gtest_color=yes

# Repeat tests (for flakiness detection)
./metabolite_test --gtest_repeat=10

# Shuffle test order with seed
./metabolite_test --gtest_shuffle --gtest_random_seed=12345
```

#### Test Coverage

**What's Tested**:
- ✅ Constructor initialization with various parameter combinations
- ✅ Getters/setters for all properties (names, indices, coefficients, bounds)
- ✅ Boolean flags and state management (is_source, is_active, is_selected, is_dummy, etc.)
- ✅ Complex operations (metabolite coefficient management, formula generation, reaction comparison)
- ✅ Input/output operations (write_to stream formatting)
- ✅ Edge cases (empty names, long names, special characters, large/negative values)
- ✅ Copy constructors and object independence
- ✅ Collection operations (adding metabolites/reactions, indexing)
- ✅ Predicates (uses_react, is_ex, is_dm, is_export, is_biomass)
- ✅ Reaction types (EX: exchange, DM: demand, export, biomass, dummy)

**Not Yet Tested** (future enhancements):
- ❌ Solver classes (LPSolver, IntSolver, IncrSolver, LnsMxSolver hierarchies)
- ❌ LP solver implementations (GrbLpSolverImp, HighsLpSolverImp, LpsLpSolverImp)
- ❌ PathFinder reachability analysis algorithms
- ❌ Experiment class functionality
- ❌ File I/O (reading PLD format files)
- ❌ Multi-compartment metabolic network scenarios
- ❌ Integration tests with toy problems
- ❌ Performance benchmarks and timing
- ❌ Memory leak detection (valgrind integration)

#### Documentation

Comprehensive documentation is provided in the `unit_tests/` directory:

- **`README.md`** (380 lines): Test organization, running tests, filtering, coverage summary, Google Test features, adding new tests, CI/CD integration, troubleshooting
- **`SETUP.md`** (550 lines): Prerequisites, quick start, build configuration, test structure, troubleshooting (15+ common issues), CI/CD pipeline examples (GitHub Actions, Jenkins), best practices, performance benchmarks, future enhancements
- **`UNIT_TEST_SETUP_SUMMARY.md`**: Complete summary of what was created, coverage details, file listing

#### CI/CD Integration

Tests are ready for continuous integration. Examples:

**GitHub Actions**:
```yaml
- name: Build and Test
  run: |
    mkdir build && cd build
    cmake -DBUILD_TESTS=ON ..
    cmake --build .
    cd unit_tests && ctest --output-on-failure
```

**Jenkins**:
```groovy
stage('Test') {
    steps {
        sh 'cd build/unit_tests && ctest --output-on-failure'
    }
}
```

#### Files Created

- `unit_tests/CMakeLists.txt` - Test build configuration with FetchContent
- `unit_tests/README.md` - Comprehensive test documentation
- `unit_tests/SETUP.md` - Detailed setup and troubleshooting guide
- `unit_tests/metabolite_test.cpp` - 15 test cases for Metabolite class
- `unit_tests/reaction_test.cpp` - 33 test cases for Reaction class
- `unit_tests/scenario_test.cpp` - 13 test cases for Scenario class
- `unit_tests/solution_test.cpp` - 15 test cases for Solution class
- `unit_tests/multisol_test.cpp` - 5 test cases for MultiSol class
- `unit_tests/params_test.cpp` - 2 test cases for Params class
- `run_unit_tests.sh` - Convenience test runner script

#### Performance

- **Build time**: ~3-5 minutes (with 6 test suites)
- **Test execution**: < 5 seconds (all 83 tests)
- **Individual suite**: < 1 second

#### Benefits

1. **Regression Prevention**: Catch bugs early before integration
2. **Documentation**: Tests serve as usage examples for classes
3. **Refactoring Confidence**: Safe to refactor with test coverage
4. **Fast Feedback**: Unit tests run in seconds vs. minutes for integration tests
5. **Focused Debugging**: Failures pinpoint exact class/method
6. **Continuous Integration**: Automated testing in CI/CD pipelines
7. **Code Quality**: Encourages testable, modular design

### Reference Tests

End-to-end integration tests for PLUM executables.

```bash
# Run reference tests
cd reftest
./reftest.sh

# Generate example test file
./reftest.sh -x

# View test output
less reftest.out
```

Reference tests are defined in `reftest/reftest.dat` with bash commands (RUN) and validation (diff commands). Reference outputs stored in `reftest/ref/`.

### Test Data

- `unit_tests/`: Unit test source files (Google Test)
- `test/toy/`: Simple test cases with documentation in `test/toy/README.md`
- `reftest/data/`: Reference test input data
- Input format documented in `fileformat.md`

## Running PLUM

### Basic Usage

```bash
plum [options] <input.dat>
```

### Key Options

```
-s <seed>        Random seed (0 = current time)
-v <solver>      Solver type: CTS, CTS2, INT, INT2, INT3, INCR, INCR_INT, COMB, DUMMY
-V <flavour>     LP solver: GRB (Gurobi), LPS (lp_solve), HIGHS (HiGHS)
-p <threads>     Number of parallel threads
-t <seconds>     Time limit (0 = no limit)
-o <file>        Output solution file
-O <prefix>      Output prefix for multiple files
-sd <file>       Supply/demand data file
-a               Use absolute objective
-e               Use error bounds
-D               Use dummy reactions
```

### Example Commands

```bash
# Basic run with Gurobi continuous solver
plum -s 1 -v CTS -V GRB data/input.dat -o solution.out

# Run with HiGHS integer solver
plum -s 1 -v INT -V HIGHS data/input.dat -o solution.out

# With supply/demand data
plum -s 1 -v CTS data/input.dat -sd data/supply.dat -o solution.out
```

## File Formats

### Internal PLD Format (fileformat.md)

Three sections in order:

```
COMPART <id>
```
Defines compartments (E=External, P=Periplasm, C=Cytosol, O=Output). Order matters for metabolite/reaction specifications.

```
MET <id> <supply-1> <supply-2>... <residual-1> <residual-2>...
```
Metabolite with supply rates (from growth medium) and residual rates (export from cell). One value per compartment.

```
REACTION <id> <obj-coeff-1> <obj-coeff-2>... <compound-spec>...
```
Reaction with objective coefficients per compartment. Compound specs: `<met-id> <compartment> <rate>`.
- Compartment 0 = products, 1+ = reactants
- Transport reactions have different compartments (e.g., E→P, P→C, C→O)

### Utility Scripts

```bash
# Convert CSV files to PLD format
util/csv2pld.py [-v obj_value] react.csv met.csv [out.pld]

# Convert TSV files  
util/tsvcvt.py react.tsv met.tsv

# Convert toy problem Excel data
util/toy2dat.py input.xlsx output.dat

# Convert SBML
util/sbmlcvt.py input.xml output.pld

# Process taxonomic weights
util/normalise_tax_wgts.sh

# Create PLD file
util/makepld.sh

# Merge PLD files
util/mergepld.sh
```

## Development Notes

### Compartment Modeling

Non-transport reactions are replicated across compartments (E, P, C, O). Transport reactions represent movement between adjacent compartments. External (E) and Output (O) are physically the same location but modeled separately for balance equations.

### Cost Model

Reaction costs reflect likelihood of use:
- Lower cost = higher likelihood of being in the true pathway
- Typical scale: 1 (gene-indicated) to 15+ (low evidence)
- Gene-indicated reactions use `gene_ind_cost` parameter (default 1.0)

### Objective Functions

- **MILP approach**: Minimize number of reactions used
- **LP approach** (PLUM): Minimize total flux through reactions
- Compartment-specific objective coefficients allow differential weighting

### Reachability

Level-based analysis for pathway depth:
- Level 0: Seed metabolites (available in growth medium)
- Level n: Reactions using only metabolites from levels < n
- Used to discriminate between solutions (topological activation)

### Lime Library

PLUM depends on the `lime` library (CSIRO internal, included as git submodule). Provides utilities:
- `lime/opts.h`: Command-line option parsing
- `lime/numutil.h`: Numerical utilities with tolerance handling
- `lime/dig.h`: Diagnostic output
- `lime/error.h`: Error handling
- `lime/strutil.h`: String utilities
- `lime/timekeeper.h`: Timing/profiling

## Documentation

### Doxygen API Documentation

All C/C++ source and header files are fully annotated with Doxygen documentation:

**Coverage:**
- 58 files documented (29 headers + 29 source files)
- File-level documentation (`@file`, `@brief`)
- Class/struct documentation with detailed descriptions
- Function documentation with `@param` and `@return` tags
- Member variables, enums, typedefs, and namespaces

**Generating Documentation:**

```bash
# Create a Doxyfile if needed
doxygen -g

# Configure Doxyfile (key settings):
# PROJECT_NAME = "PLUM - Metabolic Gap-Filling Solver"
# INPUT = include/mosh src
# RECURSIVE = YES
# EXTRACT_ALL = YES
# GENERATE_HTML = YES

# Generate HTML documentation
doxygen

# View documentation
open html/index.html  # or use your browser
```

**Key Documented Components:**
- Core classes: `Scenario`, `Reaction`, `Metabolite`, `Solution`, `MultiSol`, `Params`
- Solver hierarchy: `LPSolver`, `IntSolver`, `IncrSolver`, `LnsMxSolver`
- Solver implementations: `GrbLpSolverImp`, `HighsLpSolverImp`, `LpsLpSolverImp`
- Algorithms: `PathFinder` (reachability analysis), `MathHeur` (heuristics)
- All executable entry points: plum, plummx, plumsp, plumchk, plumcmp, plummerge

The documentation uses domain-appropriate terminology from metabolic modeling (flux balance analysis, stoichiometry, compartments, gap-filling) to help developers understand the biological context.

### Python Documentation (Numpy Style)

All Python utility scripts are fully documented with numpy-style docstrings:

**Coverage:**
- 7 Python files documented
- Module-level docstrings with overview and examples
- Function/method docstrings with Parameters, Returns, Raises, Examples sections
- Class docstrings with Attributes descriptions

**Documented Scripts:**
- `util/csv2pld.py`: CSV to PLD format converter for metabolic networks
- `util/sbmlcvt.py`: SBML (Systems Biology Markup Language) to PLD converter
- `util/toy2dat.py`: Excel toy problem converter for FBA data
- `util/tsvcvt.py`: TSV (Tab-Separated Values) converter for reactions and metabolites
- `test/gen.py`: Synthetic metabolic network test data generator
- `doc/gurobipy/Logistic_module_variables.py`: Gurobi optimization variables module
- `doc/gurobipy/Logistic_module_constraints.py`: Gurobi optimization constraints module

**Documentation Format:**
```python
def function_name(param1, param2):
    """Brief description of function.

    Detailed description providing context about what the function does,
    how it fits into the metabolic gap-filling workflow, and any important
    implementation details.

    Parameters
    ----------
    param1 : type
        Description of first parameter.
    param2 : type
        Description of second parameter.

    Returns
    -------
    return_type
        Description of return value.

    Raises
    ------
    ExceptionType
        When this exception is raised.

    Examples
    --------
    >>> function_name(value1, value2)
    expected_output

    Notes
    -----
    Additional important information about the function.
    """
```

All docstrings follow numpy conventions and can be used with Sphinx autodoc for generating HTML documentation.

**Generating Python Documentation:**

```bash
# Install Sphinx with numpy-style support
pip install sphinx sphinx-napoleon sphinx-rtd-theme

# Initialize Sphinx documentation
sphinx-quickstart docs

# Configure docs/conf.py:
# extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon']
# html_theme = 'sphinx_rtd_theme'

# Generate HTML documentation
cd docs
make html

# View documentation
open _build/html/index.html
```

**Command-Line Help:**

```bash
# View module documentation
pydoc util.csv2pld

# View function help
python3 -c "import sys; sys.path.insert(0, 'util'); import csv2pld; help(csv2pld.fix_name)"
```

### Documentation Coverage Summary

**Complete Coverage:**
- **C/C++ Files**: 58 files (29 headers + 29 source files) with Doxygen documentation
- **Python Files**: 7 files with numpy-style docstrings
- **UML Diagrams**: 6 comprehensive diagrams with 5 documentation files
- **Total Lines**: ~13,400+ lines of professional documentation
- **Standards**: Doxygen (C++), Numpy/Sphinx (Python), PlantUML (diagrams), PEP 257 compliant
- **Coverage**: 100% of source code files, complete architecture visualization

**Documentation Benefits:**
- IDE autocomplete and inline help
- HTML documentation generation (Doxygen/Sphinx)
- Domain-specific terminology and examples
- Clear parameter types and return values
- Practical usage examples throughout

### UML Diagrams

Comprehensive UML diagrams document the PLUM architecture and execution flow from bash script to C++ class level:

**Location**: `doc/uml/`

**Diagram Types:**

1. **Behavioral Diagrams** (Execution Flow):
   - `reftest_sequence.puml`: Sequence diagram showing temporal execution flow (315 lines)
   - `reftest_communication.puml`: Communication diagram emphasizing object relationships (239 lines)
   - `reftest_communication_compact.puml`: High-level architecture overview (188 lines)

2. **Structural Diagrams** (Class Architecture):
   - `class_solver_hierarchy.puml`: Complete solver framework with GapSolver and LPSolverImp hierarchies (343 lines)
   - `class_data_model.puml`: Metabolic network data structures (Scenario, Metabolite, Reaction, etc.) (385 lines)
   - `class_system_overview.puml`: System-wide architecture with 6 logical packages (357 lines)
   - `GapSolver-To-HighsLPSolverImp.puml`: Simplified single-path "vertical slice" of the solver hierarchy, isolating just the classes exercised by `plum -v CTS -V HIGHS` (GapSolver → LPSolver → LPSolverImp → HighsLPSolverImp). All sibling solvers and backends are omitted for focus; see `class_solver_hierarchy.puml` for the full picture.

**Documentation:**
- `doc/uml/README.md`: Complete user guide with installation and generation instructions
- `doc/uml/DIAGRAM_GUIDE.md`: In-depth comparison of behavioral diagram types
- `doc/uml/CLASS_DIAGRAM_GUIDE.md`: Comprehensive class diagram reading guide
- `doc/uml/INDEX.md`: Quick navigation and diagram index
- `doc/uml/SUMMARY.md`: Creation story and feasibility analysis

**Generating Diagrams:**

The recommended way to render on this machine is the `render.sh` helper, which
requires **no installation** — it uses the bundled PlantUML jar plus PlantUML's
built-in pure-Java layout engine (Smetana), so Graphviz/`dot` is **not** needed.

```bash
# Render every *.puml in doc/uml to SVG (output in doc/uml/svg/)
doc/uml/render.sh

# Render a single diagram
doc/uml/render.sh GapSolver-To-HighsLPSolverImp.puml
```

How `render.sh` avoids installs:
- **java**: already present at `/usr/bin/java` (run, not installed)
- **PlantUML**: the git-ignored `plantuml-*.jar` in the repo root (a jar is run
  via `java -jar`, not installed); `render.sh` also falls back to
  `~/apps/plant_uml/plantuml-*.jar`
- **Graphviz**: bypassed via the `!pragma layout smetana` line near the top of
  each `.puml`. This pragma is honored by the PlantUML engine regardless of
  entry point, so the same file also previews in the **VS Code PlantUML
  extension** with no Graphviz. Note: Smetana produces slightly different box/
  arrow placement than Graphviz (layout only, not content) and may emit benign
  "spline routing" warnings on dense diagrams.

Direct invocation (equivalent to what `render.sh` runs):

```bash
cd doc/uml
java -jar ../../plantuml-*.jar -tsvg -o svg *.puml    # SVG (recommended)
java -jar ../../plantuml-*.jar -o svg *.puml          # PNG
java -jar ../../plantuml-*.jar -tpdf -o svg *.puml    # PDF
```

If PlantUML and Graphviz are installed system-wide (e.g. `sudo apt-get install
plantuml graphviz`), plain `plantuml *.puml` also works; the `smetana` pragma is
still honored and simply skips Graphviz.

**Online Viewing** (no installation):
1. Visit https://www.plantuml.com/plantuml/uml/
2. Copy contents of any `.puml` file
3. Paste and render

**When to Use:**
- **Sequence diagram**: Debugging, tracing execution, understanding algorithm flow
- **Communication diagrams**: Architecture review, identifying dependencies, refactoring
- **Class diagrams**: Understanding structure, implementing features, extending components

**Coverage:**
- 7 comprehensive diagrams (~2,100 lines of PlantUML)
- 5 documentation files (~2,700 lines)
- Total: ~4,600 lines of UML documentation
- Traces execution from bash script through plum executable to C++ class level
- Documents 30+ classes, 200+ methods, 150+ attributes
- Shows design patterns: Strategy, Factory, Composite, Facade

### Session Exports

Documentation session summaries are available:
- `python_documentation_session.md`: Complete summary of numpy-style docstring addition to all Python files (June 15, 2026)
- `complete_documentation_session.md`: Comprehensive export of entire documentation project including C/C++ Doxygen, Python numpy-style, and CLAUDE.md creation (June 15, 2026)

These exports provide detailed information about documentation workflows, statistics, examples, and verification.

## Workflows

PLUM includes saved workflows for common documentation tasks:

- `/doxygen-annotate`: Add Doxygen documentation to C/C++ files (saved workflow)
- `/numpy-docstring-annotate`: Add numpy-style docstrings to Python files (saved workflow)

These workflows use multi-agent parallel processing to efficiently document large codebases.

## Project Structure

- `src/`: Implementation files for solvers, scenarios, solutions (Doxygen documented)
- `include/mosh/`: Header files with full API documentation (Doxygen documented)
- `unit_tests/`: Unit tests using Google Test framework (83 test cases across 6 test suites)
  - `metabolite_test.cpp`, `reaction_test.cpp`, `scenario_test.cpp`
  - `solution_test.cpp`, `multisol_test.cpp`, `params_test.cpp`
  - `README.md`: Test organization and usage guide
  - `SETUP.md`: Detailed setup and troubleshooting
- `util/`: Conversion and preprocessing scripts (Python/Bash with numpy-style docstrings)
- `test/`: Test data and scripts (Python scripts documented)
- `reftest/`: Reference tests with expected outputs (reftest.sh annotated with bash documentation)
- `doc/`: Documentation including `overview/overview.md` (project background)
  - `doc/uml/`: UML diagrams (sequence, communication, class diagrams with comprehensive guides)
  - `doc/gurobipy/`: Gurobi optimization example code (documented)
- `data/`: Input data files (not in repository)
- `lime/`: Submodule dependency (CSIRO utility library)
- `.claude/workflows/`: Saved documentation workflows
- `run_unit_tests.sh`: Convenience script to run all unit tests

## Quick Reference

### Common Development Tasks

**Building:**
```bash
mkdir build && cd build
cmake ..
cmake --build .
```

**Unit Testing:**
```bash
# Build with tests (default)
mkdir build && cd build
cmake ..
cmake --build .

# Run all unit tests
cd unit_tests && ctest
# Or use convenience script
./run_unit_tests.sh

# Run specific test
cd build/unit_tests
./metabolite_test
./reaction_test --gtest_filter=*Flux*
```

**Reference Testing:**
```bash
cd reftest
./reftest.sh
```

**Running:**
```bash
# Basic gap-filling
plum -s 1 -v CTS -V HIGHS data/input.dat -o solution.out

# With supply/demand constraints
plum -s 1 -v INT -V GRB data/input.dat -sd data/supply.dat -o solution.out
```

**Data Conversion:**
```bash
# CSV to PLD
util/csv2pld.py reactions.csv metabolites.csv output.pld

# SBML to PLD
util/sbmlcvt.py model.xml output.pld

# Excel toy problem
util/toy2dat.py input.xlsx > output.dat
```

**Documentation:**
```bash
# Generate C++ API docs
doxygen

# Generate Python API docs
cd docs && make html

# Generate UML diagrams (no install needed; uses bundled jar + Smetana layout)
doc/uml/render.sh

# View function help
pydoc util.csv2pld
```

### Key Files to Understand First

1. **doc/uml/class_system_overview.puml** - Complete system architecture overview (start here!)
2. **doc/uml/reftest_communication_compact.puml** - High-level execution flow
3. **include/mosh/scenario.h** - Core data model for metabolic networks
4. **include/mosh/reaction.h** - Reaction representation with stoichiometry
5. **include/mosh/metabolite.h** - Metabolite representation
6. **include/mosh/solution.h** - Solution flux distribution
7. **include/mosh/lpsolver.h** - LP solver interface
8. **src/plum.cpp** - Main executable entry point
9. **unit_tests/reaction_test.cpp** - Example usage of Reaction class (33 test cases)
10. **util/csv2pld.py** - Most commonly used data converter
11. **fileformat.md** - PLD file format specification
12. **doc/overview/overview.md** - Project background and methodology

### Important Concepts

- **Compartments**: E (External), P (Periplasm), C (Cytosol), O (Output)
- **PLD Format**: Internal representation of metabolic networks
- **Gap-Filling**: Finding minimal reaction sets to enable target metabolic functions
- **Reachability**: Level-based analysis of metabolic pathway depth
- **Flux Balance Analysis**: Steady-state optimization of metabolic fluxes
- **Objective Coefficients**: Reaction costs for optimization (lower = more likely)

### Getting Help

- **Start with UML diagrams** in `doc/uml/` for visual architecture overview
  - `README.md`: Complete guide to generating and using diagrams
  - `INDEX.md`: Quick navigation by role or task
  - Generate with: `doc/uml/render.sh` (no install needed; bundled jar + Smetana layout)
- **Unit test examples** in `unit_tests/` show how to use core classes
  - `reaction_test.cpp`: 33 examples of using Reaction class
  - `metabolite_test.cpp`: 15 examples of using Metabolite class
  - `scenario_test.cpp`: 13 examples of building metabolic networks
  - Run tests to see expected behavior: `./run_unit_tests.sh`
- All C/C++ code has Doxygen comments - generate HTML docs or use IDE hover
- All Python scripts have numpy-style docstrings - use `pydoc` or `help()`
- Test documentation: `unit_tests/README.md` and `unit_tests/SETUP.md`
- Session documentation in `python_documentation_session.md`
- Project overview in `doc/overview/overview.md`
- File format specification in `fileformat.md`
