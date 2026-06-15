/**
 * @file highsutil.h
 * @brief Utility functions for HiGHS linear programming solver integration
 *
 * This file provides helper utilities for working with the HiGHS LP/MIP solver,
 * including status checking and error handling for optimization operations used
 * in metabolic network gap-filling and flux balance analysis.
 */
#pragma once

#include "highs/Highs.h"

/**
 * @namespace mosh
 * @brief Core namespace for metabolic optimization and simulation framework
 *
 * Contains utilities and functions for metabolic network analysis, gap-filling,
 * and flux balance analysis using linear programming techniques.
 */
namespace mosh
{
    /**
     * @brief Checks the status returned by HiGHS solver and handles errors
     *
     * Validates the HighsStatus code from HiGHS LP/MIP solver operations.
     * If the status indicates an error or warning condition, appropriate
     * error handling is performed. This is used to ensure optimization
     * problems in metabolic gap-filling complete successfully.
     *
     * @param status The HighsStatus code returned by a HiGHS solver operation
     * @param where Descriptive string indicating where the check is being performed,
     *              used for error reporting and debugging context
     *
     * @throws May throw or terminate on critical solver failures
     */
    void check_highs_status (HighsStatus status, std::string where);
}
