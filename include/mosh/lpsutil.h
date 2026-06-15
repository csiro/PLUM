/**
 * @file lpsutil.h
 * @brief Utility functions for LP solver return value checking
 *
 * This file provides utility functions for error handling and validation
 * of linear programming solver operations used in metabolic network analysis
 * and flux balance analysis computations.
 */

#pragma once

/**
 * @brief Main namespace for metabolic optimization and solver handling utilities
 */
namespace mosh
{
    /**
     * @brief Checks the return value from an LP solver operation and handles errors
     *
     * Validates the return code from a linear programming solver call and throws
     * an appropriate exception or logs an error if the operation failed. This function
     * is used throughout the flux balance analysis workflow to ensure solver operations
     * complete successfully.
     *
     * @param retval The return value from the LP solver operation to validate
     * @param where A descriptive string indicating the context or location where the
     *              LP solver was called, used for error reporting and debugging
     *
     * @throws std::runtime_error if the LP solver return value indicates an error condition
     */
    void check_lps_return (unsigned char retval, std::string where);
}
