#!/bin/bash

################################################################################
# reftest.sh - Reference test runner for PLUM metabolic gap-filling solver
#
# Description:
#   Executes a series of regression tests for PLUM by running commands from
#   a test definition file and validating outputs against reference files.
#   The test file contains bash RUN commands (to execute PLUM) and diff
#   commands (to validate output against expected results in ref/).
#
# Usage:
#   ./reftest.sh [-x] [testfile]
#
# Options:
#   -x          Print an example test file and exit
#   -?, --help  Display usage information
#
# Arguments:
#   testfile    File containing RUN and TEST commands [default: reftest.dat]
#
# Test File Format:
#   - Lines are executed as bash commands via eval
#   - Empty lines create section separators in output
#   - Typical pattern: run PLUM command, then diff against ref/ files
#   - Example:
#       plum -s 1 -v CTS data/input.dat -o out/test.out
#       diff out/test.out ref/test.out
#
# Output:
#   Results written to reftest.out with timestamps and command echoing.
#   View with: less reftest.out
#
# Exit Status:
#   0 on success, 1 on usage error or missing test file
#
# Dependencies:
#   - PLUM executables (plum, plummx, etc.) must be in PATH or ../build/
#   - colour command for terminal output (optional, for colored messages)
#   - Test data files in data/ directory
#   - Reference outputs in ref/ directory
#
# Author: CSIRO
# Project: PLUM - Metabolic Gap-Filling Solver
################################################################################

################################################################################
# Usage - Display help message
#
# Prints usage information and command-line options to stderr, then exits.
#
# Arguments:
#   $@  Optional error message to display before usage
#
# Exit Status:
#   Always exits with status 1
################################################################################
Usage ()
{
    if [ $# -gt 0 ]; then echo $@ >& 2; fi
    echo "Usage: `basename $0` [-x] [reftest.sh]" >& 2
    echo "       Run tests " >&2
    echo "       Input file contains bash commands to do runs "
    echo "          and diff commands to test output" >&2
    echo "  Switches" >&2
    echo "    -x: Output an example reftest.sh" >&2
    echo "  Args" >&2
    echo "    testfile: RUN and TEST commands [reftest.dat]" >&2
    exit 1
}

################################################################################
# Example - Generate example test file
#
# Outputs a sample test definition file showing the expected format for
# reftest.dat. Demonstrates typical test patterns including:
#   - Directory setup/cleanup
#   - PLUM solver invocations with different solver types (CTS, INT)
#   - Output validation using diff commands
#
# Exit Status:
#   Always exits with status 1 (example output only)
################################################################################
Example () {
    cat << __EOF__
# Prep
if [ ! -d out ]; then mkdir out; else rm -v out/*; fi

# Test 1 - CTS solver
plum -s 1 -O out/res data/PA01.dat -sd data/sd-4abut.dat -o out/test1.out
diff out/test1.out ref/test1.out

# Test 2 - INT solver
plum -s 1 -O out/res data/PA01.dat -sd data/sd-4abut.dat -o out/test2.out -v INT
diff out/test2.out ref/test2.out
__EOF__
    exit 1
}

################################################################################
# run - Echo and execute a command
#
# Prints the command to stdout before executing it. Useful for logging
# command execution in test output.
#
# Arguments:
#   $@  Command and arguments to execute
#
# Returns:
#   Exit status of the executed command
################################################################################
function run() {
    echo "$@"
    "$@"
}

################################################################################
# Command-line argument processing
#
# Parse switches before positional arguments. Recognized options:
#   -x       : Generate example test file
#   -?, --help : Display usage information
################################################################################
while [ "${1:0:1}" = "-" ]
do
  arg="$1"
  shift
  case $arg in
      "-x")  Example ;;
      "-?")  Usage ;;
  "--help")  Usage ;;
      *) Usage "Unknown option" $arg ;;
  esac
done

################################################################################
# Initialize test execution
#
# Set up test file path (default: reftest.dat) and output file (reftest.out).
# Validate that the test file exists and is readable.
################################################################################
testfile="${1:-reftest.dat}"
if [ ! -r $testfile ]
then
    echo "Couldn't open file $testfile" >&2
    exit
fi

out="reftest.out"
echo "reftest run on `date`" > $out

################################################################################
# Main test execution loop
#
# Read test file line by line and execute each command:
#   - Empty lines: Add section separator to output
#   - Non-empty lines: Echo command (in yellow), execute via eval, and
#     capture output to both stdout and reftest.out using tee
#
# Commands are executed in the current shell context, allowing for:
#   - Environment variable persistence across commands
#   - Shell built-ins and control structures
#   - Output redirection and piping
################################################################################
while read -r line
do
    if [ "$line" == "" ]
    then
        echo "================================================" >> $out
    else
        colour YELLOW; echo $line; colour RESET
        echo "------------------------------------------------" >> $out
        echo $line >> $out
        eval "$line" | tee -a $out
    fi
done < $testfile

################################################################################
# Test completion message
#
# Display final status and instructions for viewing test results.
################################################################################
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
colour YELLOW
echo "Done. Output in $out. View with"
colour RESET
echo "less $out"
