/**
 * @file constants.h
 * @brief Core constants and enumerations for the MOSH metabolic optimization framework
 *
 * This file defines fundamental constants used in metabolic flux balance analysis,
 * including flux bounds, sentinel values, solver status codes, and supported solver types.
 */

#pragma once

/**
 * @brief Main namespace for the MOSH (Metabolic Optimization Solver Harness) library
 *
 * Contains all constants, types, and utilities for metabolic gap-filling and
 * flux balance analysis operations.
 */
namespace mosh
{
    /** @brief Upper bound for dummy reaction fluxes in gap-filling algorithms
     *
     * Dummy reactions are artificial reactions added during metabolic gap-filling
     * to identify missing pathways. This constant limits their maximum flux value
     * to prevent unrealistic solutions.
     */
    constexpr double DUMMY_FLUX_UB = 1000.0f;
    
    /** @brief Sentinel value indicating a parameter or constraint should be ignored */
    constexpr int IGNORE = 9999;
    /** @brief Sentinel value indicating a parameter value is not constrained */
    constexpr int DONT_CARE = 9999;

    /**
     * @brief Checks if a value matches the DONT_CARE sentinel
     * @param val The value to check
     * @return true if the value indicates a don't-care constraint, false otherwise
     */
    bool isDontCare (double val);
    /**
     * @brief Checks if a value matches the IGNORE sentinel
     * @param val The value to check
     * @return true if the value should be ignored, false otherwise
     */
    bool isIgnoreVal (double val);

    /**
     * @enum StatusEnum
     * @brief Optimization solver status codes
     *
     * Indicates the outcome of a metabolic flux balance analysis or gap-filling optimization.
     *
     * @var StatusEnum::UNSOLVED
     * Problem has not been solved yet
     *
     * @var StatusEnum::CTS_OPTIMAL
     * Continuous optimization found an optimal solution
     *
     * @var StatusEnum::INT_OPTIMAL
     * Integer/mixed-integer optimization found an optimal solution
     *
     * @var StatusEnum::SUB_OPTIMAL
     * Solver found a feasible but suboptimal solution
     *
     * @var StatusEnum::INFEASIBLE_
     * Problem is infeasible (no solution exists)
     *
     * @var StatusEnum::RANDOM
     * Random or undefined status
     */
    enum StatusEnum {
        UNSOLVED, CTS_OPTIMAL, INT_OPTIMAL, SUB_OPTIMAL, INFEASIBLE_, RANDOM
    };
    /**
     * @enum Flavour
     * @brief Supported linear/mixed-integer programming solver backends
     *
     * @var Flavour::GUROBI
     * Gurobi optimization solver
     *
     * @var Flavour::LP_SOLVE
     * lp_solve open-source solver
     *
     * @var Flavour::HIGHS
     * HiGHS high-performance optimization solver
     */
    enum Flavour {GUROBI, LP_SOLVE, HIGHS};
}
