/**
 * @file lpsutil.cpp
 * @brief Utility functions for lp_solve linear programming solver integration
 *
 * Provides error checking and validation for return values from the lp_solve
 * library, which is used as one of the LP/MILP solver backends for metabolic
 * gap-filling and flux balance analysis computations in PLUM.
 */

#include "lime/error.h"

#include "mosh/lpsutil.h"

using namespace std;
using namespace mosh;

/**
 * @brief Validates lp_solve solver return values and crashes on error
 *
 * Checks if the lp_solve library call succeeded by verifying the return value
 * equals 1 (success). If the return value indicates failure, terminates the
 * program with an error message indicating where the failure occurred.
 *
 * This function is used throughout the LPS solver implementation to ensure
 * that linear programming operations complete successfully during metabolic
 * flux analysis.
 *
 * @param retval The return value from an lp_solve library call (1 = success)
 * @param where Description of the lp_solve operation that was attempted
 *
 * @throws Calls limeCrash() to terminate execution if retval != 1
 */
void mosh::check_lps_return (unsigned char retval, std::string where)
{
    if (retval != 1) 
        limeCrash ("Return val " + to_string(retval) + " from " + where);
}
