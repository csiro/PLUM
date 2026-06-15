/**
 * @file highsutil.cpp
 * @brief Utility functions for HiGHS linear programming solver integration
 *
 * This file provides helper functions for checking and handling status codes
 * returned by the HiGHS solver during flux balance analysis and metabolic
 * gap-filling operations.
 */

#include "lime/error.h"

#include "mosh/highsutil.h"


using namespace std;
using namespace mosh;

/**
 * @brief Check and handle HiGHS solver status codes
 *
 * Evaluates the status returned by HiGHS solver operations and issues
 * appropriate warnings for non-OK statuses. This function provides
 * centralized error handling for all HiGHS solver interactions during
 * linear programming operations for metabolic network analysis.
 *
 * @param status The HighsStatus code returned by a HiGHS operation
 * @param where A descriptive string indicating the operation context
 *              (e.g., "addVar", "changeColCost") for debugging purposes
 */
void mosh::check_highs_status (HighsStatus status, std::string where)
{
    switch (status) {
    case HighsStatus::kOk:
        // Do nothing
        break;
    case HighsStatus::kWarning:
        limeWarning ("HiGHS warning from " + where);
        break;
    case HighsStatus::kError:
        limeWarning ("HiGHS error from " + where);
        break;
    }
}
