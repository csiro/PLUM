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

### Reference Tests

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

## Project Structure

- `src/`: Implementation files for solvers, scenarios, solutions (Doxygen documented)
- `include/mosh/`: Header files with full API documentation (Doxygen documented)
- `util/`: Conversion and preprocessing scripts (Python/Bash)
- `test/`: Test data and scripts
- `reftest/`: Reference tests with expected outputs
- `doc/`: Documentation including `overview/overview.md` (project background)
- `data/`: Input data files (not in repository)
- `lime/`: Submodule dependency (CSIRO utility library)
