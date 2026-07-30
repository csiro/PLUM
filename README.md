# PLUM

**PLUM** is a metabolic gap-filling solver for Flux Balance Analysis (FBA). Given a
metabolic network with known inputs and outputs, it determines which reactions from a
reaction database should be included to balance the network — or to maximize biomass
production.

## Repository Status:
 - ONHOLD 2026-07-30: Need to be rebuild in the proper environment 

## Why PLUM?

Traditional gap-filling relies on Mixed Integer Linear Programming (MILP), which can take
hours or weeks on realistic networks. PLUM's key innovation is to use **LP-based flux
minimization** instead, achieving polynomial-time solutions (seconds rather than
hours/weeks) while still incorporating:

- **Compartment modeling** — External (E), Periplasm (P), Cytosol (C), and Output (O)
- **Reachability constraints** — level-based analysis of metabolic pathway depth

## How It Works

PLUM models a metabolic network as a set of metabolites and stoichiometric reactions
distributed across compartments. Non-transport reactions are replicated across
compartments, while transport reactions move metabolites between adjacent ones.

- **MILP approach:** minimize the *number* of reactions used
- **LP approach (PLUM):** minimize *total flux* through reactions

Each reaction carries a cost reflecting the likelihood it belongs to the true pathway
(lower cost = more likely; gene-indicated reactions get the lowest cost). Reachability
analysis assigns levels to metabolites and reactions (level 0 = seed metabolites available
in the growth medium) to discriminate between candidate solutions.

## Architecture

**Data model**
- `Scenario` — top-level container for metabolites, reactions, and experiments
- `Metabolite` — compounds with supply/residual rates per compartment
- `Reaction` — stoichiometric transformations with objective coefficients and costs
- `Solution` / `MultiSol` — flux assignments and reaction selections

**Solver hierarchy**
- `LpSolver` (abstract base) → `LpSolverImp` implementations:
  `GrbLpSolverImp` (Gurobi), `HighsLpSolverImp` (HiGHS), `LpsLpSolverImp` (lp_solve)
- `IntSolver`, `IncrSolver`, `LnsMxSolver` (Large Neighborhood Search — the main
  algorithmic solver)

**Algorithms**
- `PathFinder` — graph-based reachability / producibility analysis
- `Params` — central configuration holder

## Executables

| Tool | Purpose |
|------|---------|
| `plum` | Main gap-filling solver |
| `plummx` | Mixed solver variant |
| `plumsp` | Solution processor / printer |
| `plumchk` | Solution checker / validator |
| `plumcmp` | Solution comparator |
| `plummerge` | Merge multiple solutions |

## Building

PLUM uses **CMake** and **C++20**.

```bash
# Initialize the lime submodule (CSIRO utility library)
git submodule update --init --recursive

# Configure and build (out-of-source)
mkdir build && cd build
cmake ..
cmake --build .
```

PLUM supports three LP/MILP backends, all compiled in via preprocessor flags. At least one
must be configured for linking to succeed:

- **Gurobi** — set `GUROBI_HOME`
- **HiGHS** — set `HIGHS_HOME`
- **lp_solve** — set `LP_SOLVE_HOME`

```bash
GUROBI_HOME=/path/to/gurobi HIGHS_HOME=/path/to/highs cmake ..
```

A static build is available with `cmake -DSTATIC=ON ..`.

See **`HOW-TO-BUILD-PLUM.md`** for a detailed walkthrough.

## Usage

```bash
plum [options] <input.dat>
```

Common options:

```
-s <seed>     Random seed (0 = current time)
-v <solver>   Solver type: CTS, INT, INCR, COMB, ...
-V <flavour>  LP backend: GRB (Gurobi), HIGHS, LPS (lp_solve)
-p <threads>  Parallel threads
-t <seconds>  Time limit (0 = no limit)
-o <file>     Output solution file
-sd <file>    Supply/demand data file
```

Examples:

```bash
# Basic gap-filling with the HiGHS continuous solver
plum -s 1 -v CTS -V HIGHS data/input.dat -o solution.out

# Integer solver via Gurobi, with supply/demand constraints
plum -s 1 -v INT -V GRB data/input.dat -sd data/supply.dat -o solution.out
```

## File Format

PLUM's internal **PLD format** has three ordered sections (full spec in `fileformat.md`):

```
COMPART <id>                              # define compartments (E, P, C, O)
MET <id> <supply...> <residual...>        # one value per compartment
REACTION <id> <obj-coeff...> <specs...>   # spec: <met-id> <compartment> <rate>
```

Converters in `util/` produce PLD from other formats:

```bash
util/csv2pld.py reactions.csv metabolites.csv output.pld   # CSV
util/sbmlcvt.py model.xml output.pld                       # SBML
util/toy2dat.py input.xlsx > output.dat                    # Excel toy problem
```

## Testing

**Unit tests** (Google Test, 83 cases across 6 suites) are built with PLUM by default:

```bash
cd build && cmake .. && cmake --build .
cd unit_tests && ctest --output-on-failure
# or, from the project root:
./run_unit_tests.sh
```

Disable with `cmake -DBUILD_TESTS=OFF ..`.

**Reference (integration) tests** exercise the executables end-to-end:

```bash
cd reftest && ./reftest.sh
```

## Documentation

- **`doc/uml/`** — architecture and execution-flow UML diagrams (start with
  `class_system_overview.puml`). Render with `doc/uml/render.sh` — no install needed;
  it uses the bundled PlantUML jar plus the pure-Java Smetana layout engine.
- **Doxygen** — all C/C++ files are annotated; run `doxygen` to generate HTML.
- **Python** — utility scripts use numpy-style docstrings (Sphinx-ready).
- **`doc/overview/overview.md`** — project background and methodology.

## Project Structure

```
src/          Solver, scenario, and solution implementations
include/mosh/  Public headers (fully Doxygen documented)
unit_tests/   Google Test suites
reftest/      Reference/integration tests and expected outputs
util/         Data conversion scripts (CSV, SBML, TSV, Excel → PLD)
test/         Test data and generators
doc/          Documentation, UML diagrams, overview
data/         Input data files (not tracked in the repository)
lime/         CSIRO utility library (git submodule)
```

## Key Concepts

- **Compartments** — E (External), P (Periplasm), C (Cytosol), O (Output)
- **Gap-filling** — finding a minimal reaction set to enable target metabolic functions
- **Reachability** — level-based analysis of pathway depth
- **Flux Balance Analysis** — steady-state optimization of metabolic fluxes
- **Cost model** — reaction cost encodes likelihood of use (1 = gene-indicated, 15+ = low
  evidence)
