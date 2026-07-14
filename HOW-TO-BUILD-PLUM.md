# HOW TO BUILD PLUM

## Step-by-Step Build Instructions

This guide provides complete instructions for building PLUM (Metabolic Gap-Filling Solver) from source using CMake and C++20.

## Prerequisites

### Required Software

- **C++ Compiler**: GCC 10+ or Clang 11+ with C++20 support  
  -- for kubeflow VSCode Notebook image June 2026; $gcc --version -> 11.4.0 
- **CMake**: Version 3.15 or higher  
  -- for kubeflow VSCode Notebook image June 2026; $cmake --version -> cmake version 3.22.1
- **Git**: For cloning and submodule management  
  -- for kubeflow VSCode Notebook image June 2026; $git --version -> git version 2.34.1

### Required LP/MILP Solver (at least one)

PLUM requires at least one of the following LP/MILP solvers:

- [**Gurobi**: Commercial optimizer requires license](https://www.gurobi.com/why-gurobi)
- [**HiGHS**: Open-source LP/MILP solver](https://ergo-code.github.io/HiGHS/stable/installation/)
- [**lp_solve**: Open-source LP solver](https://lp-solve.github.io/)

## Step 1: Clone the Repository

```bash
# Clone PLUM repository
git clone <repository-url> PLUM
cd PLUM
```

## Step 2: Initialize Submodules

PLUM depends on the `lime` library (CSIRO utility library) included as a git submodule.

```bash
# Initialize and fetch all submodules
git submodule update --init --recursive
```

This will clone the `lime` library into the `lime/` directory.

## Step 3: Configure Solver Environment Variables

Set environment variables for the LP/MILP solvers you have installed. **You must configure at least one solver** for PLUM to build and link successfully.

### Option A: Gurobi

```bash
export GUROBI_HOME=/path/to/gurobi
export LD_LIBRARY_PATH=$GUROBI_HOME/lib:$LD_LIBRARY_PATH
```

Alternatively, if you have a `module.sh` script that loads Gurobi:

```bash
source module.sh
```

### Option B: HiGHS

```bash
export HIGHS_HOME=/path/to/highs
export LD_LIBRARY_PATH=$HIGHS_HOME/lib:$LD_LIBRARY_PATH
```

### Option C: lp_solve

```bash
export LP_SOLVE_HOME=/path/to/lpsolve
export LD_LIBRARY_PATH=$LP_SOLVE_HOME/lib:$LD_LIBRARY_PATH
```

### Configuring Multiple Solvers

You can configure multiple solvers simultaneously. PLUM will compile support for all configured solvers:

```bash
export GUROBI_HOME=/path/to/gurobi
export HIGHS_HOME=/path/to/highs
export LP_SOLVE_HOME=/path/to/lpsolve
export LD_LIBRARY_PATH=$GUROBI_HOME/lib:$HIGHS_HOME/lib:$LP_SOLVE_HOME/lib:$LD_LIBRARY_PATH
```

## Step 4: Create Build Directory

Use an out-of-source build (recommended to keep the source tree clean):

```bash
# Create and enter build directory
mkdir build
cd build
```

## Step 5: Configure with CMake

Run CMake to configure the build system:

```bash
# Basic configuration (uses environment variables)
cmake ..
```

### Configuration Options

**Specify solvers explicitly** (if environment variables are not set):

```bash
cmake -DGUROBI_HOME=/path/to/gurobi -DHIGHS_HOME=/path/to/highs ..
```

**Build type** (default is Release):

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
# Other options: Debug, RelWithDebInfo, MinSizeRel
```

**Static build** (statically link libraries):

```bash
cmake -DSTATIC=ON ..
```

**Custom install prefix**:

```bash
cmake -DCMAKE_INSTALL_PREFIX=/custom/install/path ..
```

**Combined example**:

```bash
cmake -DCMAKE_BUILD_TYPE=Release \
      -DGUROBI_HOME=/opt/gurobi \
      -DHIGHS_HOME=/usr/local/highs \
      -DCMAKE_INSTALL_PREFIX=$HOME/plum-install \
      ..
```

## Step 6: Build PLUM

Compile all executables:

```bash
# Build all targets (uses all available cores)
cmake --build .
```

### Build Options

**Parallel build** (specify number of jobs):

```bash
cmake --build . -j 8
```

**Build specific target** (executable):

```bash
# Build only the main plum executable
cmake --build . --target plum

# Build other executables
cmake --build . --target plummx
cmake --build . --target plumsp
cmake --build . --target plumchk
```

**Build with specific configuration** (if using multi-config generator):

```bash
cmake --build . --config Debug
```

**Verbose build** (see full compiler commands):

```bash
cmake --build . --verbose
# or
make VERBOSE=1
```

## Step 7: Verify Build Success

After a successful build, executables will be in the `build/` directory:

```bash
# List built executables
ls -lh plum plummx plumsp plumchk plumcmp plummerge

# Check executable
./plum --help
```

## Step 8: Install (Optional)

Install PLUM executables to the system or custom prefix:

```bash
# Install to CMAKE_INSTALL_PREFIX (default: /usr/local)
cmake --install .

# Install to custom prefix
cmake --install . --prefix $HOME/plum-install
```

After installation, executables will be in `<prefix>/bin/`.

## Step 9: Run Tests

Verify the build by running reference tests:

```bash
# From build directory, navigate to reftest
cd ../reftest

# Run all reference tests
./reftest.sh

# View test output
less reftest.out
```

Tests compare PLUM output against reference solutions in `reftest/ref/`. All tests should pass.

## Build Types Explained

| Build Type      | Optimization | Debug Info | Use Case                          |
|-----------------|--------------|------------|-----------------------------------|
| Release         | -O3          | No         | Production use (default)          |
| Debug           | -O0          | Yes        | Development, debugging with gdb   |
| RelWithDebInfo  | -O2          | Yes        | Profiling, debugging optimized    |
| MinSizeRel      | -Os          | No         | Minimize binary size              |

## Troubleshooting

### CMake Can't Find Solver

**Error**: `Could not find Gurobi/HiGHS/lp_solve`

**Solution**: Ensure environment variables are set correctly and libraries are installed:

```bash
# Verify paths exist
ls $GUROBI_HOME/lib
ls $HIGHS_HOME/lib

# Set explicitly in CMake command
cmake -DGUROBI_HOME=/path/to/gurobi ..
```

### No Solver Configured

**Error**: `At least one solver must be configured`

**Solution**: Configure at least one of GUROBI_HOME, HIGHS_HOME, or LP_SOLVE_HOME.

### Linking Errors

**Error**: `undefined reference to Gurobi/HiGHS functions`

**Solution**: Ensure `LD_LIBRARY_PATH` includes solver library directories:

```bash
export LD_LIBRARY_PATH=$GUROBI_HOME/lib:$LD_LIBRARY_PATH
```

### Submodule Not Initialized

**Error**: `lime library not found`

**Solution**: Initialize submodules:

```bash
git submodule update --init --recursive
```

### C++20 Not Supported

**Error**: Compiler errors about C++20 features

**Solution**: Update compiler to GCC 10+ or Clang 11+:

```bash
# Check compiler version
g++ --version
clang++ --version

# Use specific compiler
cmake -DCMAKE_CXX_COMPILER=g++-11 ..
```

## Quick Build Script

For convenience, here's a complete build script:

```bash
#!/bin/bash
# build_plum.sh - Quick PLUM build script

set -e  # Exit on error

# 1. Initialize submodules
echo "Initializing submodules..."
git submodule update --init --recursive

# 2. Load solver environment (adjust paths as needed)
export GUROBI_HOME=/opt/gurobi
export HIGHS_HOME=/usr/local/highs
export LD_LIBRARY_PATH=$GUROBI_HOME/lib:$HIGHS_HOME/lib:$LD_LIBRARY_PATH

# 3. Create build directory
echo "Creating build directory..."
mkdir -p build
cd build

# 4. Configure
echo "Configuring with CMake..."
cmake -DCMAKE_BUILD_TYPE=Release ..

# 5. Build
echo "Building PLUM..."
cmake --build . -j $(nproc)

# 6. Verify
echo "Build complete! Executables:"
ls -lh plum plummx plumsp plumchk plumcmp plummerge

echo "Done! Run './plum --help' to get started."
```

Make it executable and run:

```bash
chmod +x build_plum.sh
./build_plum.sh
```

## Next Steps

After building PLUM successfully:

1. **Read file format documentation**: `fileformat.md`
2. **Explore test data**: `test/toy/README.md`
3. **Run example**: 
   ```bash
   ./build/plum -s 1 -v CTS -V HIGHS test/toy/input.dat -o solution.out
   ```
4. **Convert data**: Use utilities in `util/` to convert CSV/SBML to PLD format
5. **View UML diagrams**: `doc/uml/` for architecture overview

## Summary

```bash
# Minimal build commands
git submodule update --init --recursive
export GUROBI_HOME=/path/to/gurobi  # or HIGHS_HOME or LP_SOLVE_HOME
mkdir build && cd build
cmake ..
cmake --build .
./plum --help
```

For detailed usage information, see `CLAUDE.md` and the generated Doxygen documentation.
