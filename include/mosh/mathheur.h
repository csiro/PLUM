/**
 * @file mathheur.h
 * @brief Mathematical heuristic solver for metabolic gap-filling problems
 *
 * This file defines the MathHeur class, which implements a mathematical programming
 * heuristic approach for identifying minimal sets of reactions to add to a metabolic
 * network to satisfy growth requirements. The solver uses Gurobi's mixed-integer
 * programming capabilities to iteratively improve solutions.
 */

#pragma once

// Compile with
// g++ fn.cpp -lgurobi81 -lgurobi_g++5.2

#include <vector>
#include <memory>

#include "gurobi_c++.h"

#include "mosh/gapsolver.h"
#include "mosh/grblpsolver.h"

#include "lime/numutil.h"
#include "lime/rand.h"

/**
 * @brief Namespace for metabolic optimization and solver heuristics
 */
namespace mosh
{
    /**
     * @class MathHeur
     * @brief Mathematical heuristic solver for metabolic gap-filling
     *
     * MathHeur extends GapSolver to provide a mathematical programming-based heuristic
     * for solving metabolic gap-filling problems. It formulates the problem as a mixed-integer
     * linear program (MILP) using Gurobi, where binary variables indicate reaction usage and
     * continuous variables represent metabolic fluxes. The solver iteratively refines solutions
     * to find minimal reaction sets that enable biomass production.
     *
     * The approach minimizes the number of reactions added while satisfying:
     * - Stoichiometric mass balance constraints
     * - Flux bounds for each reaction
     * - Minimum biomass production threshold
     */
    class MathHeur : public GapSolver 
    {
    public:
        /**
         * @brief Constructs a MathHeur solver instance
         *
         * Initializes the mathematical heuristic solver with the given scenario and parameters.
         * Sets up the Gurobi model with appropriate optimization settings including minimization
         * objective, thread count, random seed, MIP focus on feasibility, and tolerance parameters.
         *
         * @param scenario Pointer to the metabolic scenario defining the network and gap-filling problem
         * @param params Pointer to solver parameters controlling optimization behavior
         * @param env Shared pointer to Gurobi environment for model creation
         * @param seed Random seed for reproducibility of stochastic components
         */
        MathHeur (
            Scenario* scenario, Params* params,
            GrbEnvPtr env, int seed
        ) :
            GapSolver(scenario, params, seed),
            iter_(0),
            max_iters_(params->max_mh_iters),
            model_(*env),
            env_(env),
            use_(num_reactions()),
            flux_(num_reactions()),
            biomass_()
        {
            model_.set(GRB_StringAttr_ModelName, "plum_mh");
            model_.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
            model_.set(GRB_IntParam_Threads, params->num_threads);
            model_.set(GRB_IntParam_Seed, grb_seed());
            model_.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_FEASIBILITY);
            model_.set(GRB_DoubleParam_MIPGap, params->max_mip_gap);
            model_.set(GRB_DoubleParam_FeasibilityTol, params->feasibility_tol);
            model_.set(GRB_DoubleParam_IntFeasTol, params->int_feas_tol);
            if (params->time_limit > 0) {
                double dbl_limit = (double) params->time_limit;
                model_.set(GRB_DoubleParam_TimeLimit, dbl_limit);
            }
            //limeSetEpsilon(1e-8);
        }

        /**
         * @brief Solves the metabolic gap-filling problem using mathematical heuristics
         *
         * Executes the main solving algorithm to identify a minimal set of reactions that
         * enables biomass production in the metabolic network. The method iteratively improves
         * the solution by solving MILP formulations.
         *
         * @return Shared pointer to Solution object containing the gap-filled reaction set
         */
        SolutionPtr solve() override;
        
        /**
         * @brief Accesses the flux variable for a specific reaction
         *
         * @param k Index of the reaction in the metabolic network
         * @return Reference to the Gurobi variable representing flux through reaction k
         */
        GRBVar& flux(size_t k) {return flux_[k];}
        
        /**
         * @brief Generates a summary string of the solver's execution
         *
         * @return String containing solver statistics and solution information
         */
        std::string summary() override;
        
        /**
         * @brief Writes the Gurobi model to a file for inspection or debugging
         *
         * @param filename Path to the output file (e.g., .lp, .mps, or .rlp format)
         */
        void write_model (std::string filename)
        {
            model_.write (filename);
        }

    private:
        /**
         * @brief Initializes the MILP model structure with variables and constraints
         *
         * Creates decision variables for reaction usage and fluxes, and formulates
         * stoichiometric balance constraints, flux bounds, and biomass requirements.
         */
        void set_up_model();
        /**
         * @brief Iteratively improves the solution for a specific metabolite
         *
         * Refines the current solution by focusing on constraints related to the given
         * metabolite, attempting to reduce the number of reactions needed.
         *
         * @param met Pointer to the metabolite to focus improvement efforts on
         * @param sol Reference to shared pointer of current solution (updated in-place)
         */
        void improve(const Metabolite* met, SolutionPtr& sol);
        
        int iter_; /**< Current iteration count of the heuristic algorithm */
        int max_iters_; /**< Maximum number of iterations allowed before termination */
        int num_threads_; /**< Number of threads for parallel Gurobi optimization */
        int time_limit_; /**< Time limit in seconds for optimization (0 for no limit) */
        lime::Rand rand_; /**< Random number generator for stochastic heuristic components */"
        
        GRBModel model_; /**< Gurobi MILP model for gap-filling optimization */
        GrbEnvPtr env_; /**< Shared pointer to Gurobi environment */"
        
        std::vector<GRBVar> use_; /**< Binary decision variables indicating whether each reaction is used */"
        std::vector<GRBVar> flux_; /**< Continuous variables representing flux through each reaction */"

        GRBVar biomass_; /**< Variable representing biomass production flux */"

    };
}
