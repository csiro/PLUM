/**
 * @file intsolver.h
 * @brief Integer programming solver for metabolic gap-filling problems
 *
 * This file defines the IntSolver class which extends GapSolver to solve
 * metabolic network gap-filling problems using Mixed Integer Programming (MIP).
 * The solver supports multiple formulations for optimizing reaction addition
 * while satisfying flux balance constraints.
 */
#pragma once

#include <vector>
#include <memory>

#include "mosh/gapsolver.h"
#include "mosh/lpsolverimp.h"

/**
 * @namespace mosh
 * @brief Metabolic Optimization and gap-filling Solver Harness
 */
namespace mosh
{
    /**
     * @class IntSolver
     * @brief Mixed Integer Programming solver for metabolic gap-filling
     *
     * IntSolver extends GapSolver to solve metabolic network gap-filling problems
     * using integer programming techniques. It identifies the minimal set of reactions
     * to add to a metabolic network to achieve target biomass production across
     * multiple experimental conditions.
     *
     * The solver supports two formulations:
     * - Formulation 1: Uses biomass multiplier for objective scaling
     * - Formulation 2: Direct formulation without biomass multiplier
     *
     * The solver can optionally use a continuous (CTS) solution as a warm start
     * and performs rounding to obtain integer feasible solutions.
     */
    class IntSolver : public GapSolver
    {
    public:
        /**
         * @brief Construct a new IntSolver object
         *
         * @param scenario Pointer to the metabolic network scenario containing reactions and constraints
         * @param params Pointer to solver parameters and settings
         * @param seed Random seed for reproducible stochastic behavior
         * @param which_formulation Formulation type (1 or 2) for the MIP problem
         * @param imp Shared pointer to the LP solver implementation backend
         * @param cts_sol Shared pointer to continuous solution for warm-starting (can be nullptr)
         */
        IntSolver(
            Scenario* scenario, Params* params, int seed,
            int which_formulation,
            LPSolverImpPtr imp, MultiSolPtr cts_sol
        ) :
            GapSolver(scenario, params, seed),
            which_formulation_(which_formulation),
            imp_(imp),
            cts_sol_(cts_sol),
            biomass_obj_mult_(0.0f),
            rounding_iters_(0),
            mip_gap_(-1.0f)
        {
        }

        /**
         * @brief Solve the integer gap-filling problem
         *
         * Formulates and solves the MIP problem to find the minimal set of reactions
         * to add to the metabolic network. The method builds the integer formulation,
         * optionally performs rounding iterations, and returns the optimal solution.
         *
         * @return SolutionPtr Shared pointer to the solution containing selected reactions and flux distributions
         */
        SolutionPtr solve() override;

        /**
         * @brief Generate a summary string of the solution
         *
         * Creates a human-readable summary of the integer solution including
         * selected reactions, objective value, MIP gap, and solver statistics.
         *
         * @return std::string Summary text describing the solution
         */
        std::string summary() override;

    private:
        /**
         * @brief Calculate the biomass multiplier for integer formulation
         *
         * Computes an appropriate biomass objective multiplier based on the
         * continuous solution. This multiplier is used in formulation 1 to
         * scale the biomass objective relative to reaction addition penalties.
         */
        void calc_biomass_mult();
        /**
         * @brief Formulate the MIP problem
         *
         * Constructs the mixed integer programming formulation including:
         * - Binary variables for reaction selection
         * - Continuous variables for metabolic fluxes
         * - Stoichiometric and thermodynamic constraints
         * - Objective function minimizing reaction additions
         */
        void formulate ();
        /**
         * @brief Apply rounding heuristics to improve integer solution
         *
         * Performs iterative rounding to fix fractional binary variables and
         * improve solution quality. Tracks the number of rounding iterations
         * performed in rounding_iters_.
         */
        void fix_rounding();

        /** @brief Formulation type selector (1: with biomass multiplier, 2: without) */
        int which_formulation_;
        /** @brief Shared pointer to the underlying LP/MIP solver implementation */
        LPSolverImpPtr imp_;

        /** @brief Continuous solution for each experiment used as warm start */
        MultiSolPtr cts_sol_;

        // Solution attributes
        /** @brief Biomass objective multiplier for scaling in formulation 1 */
        double biomass_obj_mult_;
        /** @brief Number of rounding iterations performed during solution */
        int rounding_iters_;
        /** @brief Final MIP optimality gap (relative difference between best bound and incumbent) */
        double mip_gap_;
    };
}
