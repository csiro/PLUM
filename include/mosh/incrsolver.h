/**
 * @file incrsolver.h
 * @brief Incremental solver for metabolic gap-filling problems
 *
 * This file contains the IncrSolver class which implements an incremental
 * approach to gap-filling by adding reactions in order of increasing cost
 * until biomass production is achieved without dummy reactions.
 */
#pragma once

#include <set>

#include "mosh/gapsolver.h"
#include "mosh/lpsolver.h"
#include "mosh/reaction.h"

/**
 * @namespace mosh
 * @brief Main namespace for metabolic optimization and simulation heuristics
 */
namespace mosh
{
    /**
     * @class IncrSolver
     * @brief Incremental solver for metabolic gap-filling problems
     *
     * IncrSolver extends GapSolver to implement an incremental approach where
     * reactions are added to the metabolic network in order of increasing cost
     * (objective coefficient). The solver iteratively adds reactions until the
     * target biomass can be produced without requiring dummy reactions.
     *
     * The algorithm proceeds in cost iterations, where each cost level is processed
     * through multiple levels of metabolite propagation. This approach allows for
     * efficient exploration of the reaction space while minimizing the total cost
     * of added reactions.
     */
    class IncrSolver : public GapSolver 
    {
    public:
        /**
         * @brief Create incremental solver
         *
         * Constructs an IncrSolver instance with the specified parameters for
         * solving metabolic gap-filling problems incrementally.
         *
         * @param scenario Pointer to the metabolic network scenario
         * @param params Pointer to solver parameters
         * @param seed Random seed for stochastic operations
         * @param imp Linear programming solver implementation for LP relaxation
         * @param imp2 Integer programming solver implementation (non-null only if integer solution required)
         * @param min_cost Minimum cost threshold for reactions to consider
         * @param selected_react_file Output stream for logging selected reactions
         */
        IncrSolver(
            Scenario* scenario, Params* params, int seed,
            LPSolverImpPtr imp, LPSolverImpPtr imp2,
            double min_cost, std::ofstream& selected_react_file
        ) :
            GapSolver(scenario, params, seed),
            lp_solver_ (scenario, params, rand_seed(), 1, imp),
            imp2_(imp2),
            return_int_sol_(imp2 != nullptr),
            time_limit_(params->time_limit),
            min_cost_(min_cost),
            max_used_cost_(min_cost),
            costs_(),
            selected_react_file_(selected_react_file),
            available_(scenario->num_metabolites(), false),
            summary_(),
            cost_iters_(0),
            level_iters_(0),
            cycle_mets_(0)
        {
            init();
        }

        /**
         * @brief Solve the metabolic gap-filling problem
         *
         * Executes the incremental solving algorithm to find a minimal set of
         * reactions that enable biomass production.
         *
         * @return Pointer to the solution containing selected reactions and objective value
         */
        SolutionPtr solve() override;
        
        /**
         * @brief Get summary statistics of the solving process
         *
         * Returns a string containing statistics about the incremental solving process,
         * including number of cost iterations, level iterations, and cycle metabolites processed.
         *
         * @return String containing summary statistics
         */
        std::string summary() override;
        
    private:
        /**
         * @brief Initialize solver data structures
         *
         * Initializes internal data structures including cost sets, metabolite
         * availability tracking, and prepares the solver for execution.
         */
        void init();
        /**
         * @brief Execute the main incremental solving algorithm
         *
         * Implements the core incremental solving logic, iterating through cost
         * levels and adding reactions until a feasible solution is found.
         *
         * @return Pointer to the solution or nullptr if no solution found
         */
        SolutionPtr do_solve();
        /**
         * @brief Process reactions at a given cost level
         *
         * Iterates through reactions at the specified cost, attempting to enable
         * biomass production by adding reactions and propagating metabolite availability.
         *
         * @param cost The cost level to process
         * @return Pointer to solution if found at this cost level, nullptr otherwise
         */
        SolutionPtr cost_iter (double cost);
        /**
         * @brief Process one level of metabolite propagation
         *
         * Performs one iteration of metabolite availability propagation at the given
         * cost level, attempting to add reactions whose substrates are now available.
         *
         * @param cost The current cost level being processed
         * @param level The propagation level (distance from initially available metabolites)
         * @return True if at least one reaction was added, false otherwise
         */
        bool level_iter (double cost, int level);
        /**
         * @brief Attempt to add a metabolite involved in cycles
         *
         * Handles special case of metabolites that participate in metabolic cycles,
         * which may require special treatment to break circular dependencies.
         *
         * @param cost The current cost level
         * @return True if a cycle metabolite was successfully processed, false otherwise
         */
        bool add_cycle_met (double cost);

        /**
         * @brief Check if a reaction can be enabled
         *
         * Determines whether a reaction can be activated based on the availability
         * of its substrate metabolites.
         *
         * @param react Pointer to the reaction to check
         * @return True if all substrates are available, false otherwise
         */
        bool is_enabled (const Reaction* react);
        /**
         * @brief Check if a reaction produces carbon-containing compounds
         *
         * Determines whether a reaction generates carbon-containing metabolites,
         * which is important for biomass production in metabolic networks.
         *
         * @param react Pointer to the reaction to check
         * @return True if the reaction produces carbon compounds, false otherwise
         */
        bool makes_carbon (const Reaction* react);
        /**
         * @brief Propagate metabolite availability from a reaction
         *
         * Updates the availability status of product metabolites when a reaction
         * is added to the network, enabling subsequent reactions.
         *
         * @param react Pointer to the reaction whose products should be marked available
         */
        void propagate (const Reaction* react);
        /**
         * @brief Check if a metabolite is available
         *
         * Queries whether a metabolite has been produced by a reaction at a
         * previous level and is therefore available for use.
         *
         * @param met Pointer to the metabolite to check
         * @return True if the metabolite is available, false otherwise
         */
        bool is_available (const Metabolite* met) const
        {
            return available_[met->index()];
        }
            

        LPSolver lp_solver_; /**< Linear programming solver for flux balance analysis */"
        LPSolverImpPtr imp2_; /**< Integer programming solver implementation (used only for integer solutions) */"
        bool return_int_sol_; /**< Flag indicating whether to return integer solution */"
        int time_limit_; /**< Maximum time allowed for solving (seconds) */"
        double min_cost_; /**< Minimum cost threshold for considering reactions */"
        double max_used_cost_; /**< Maximum cost of reactions actually used in the solution */"
        std::set<double> costs_; /**< Set of unique reaction costs for incremental processing */"
        lime::Rand rand_; /**< Random number generator for stochastic operations */"
        std::ofstream& selected_react_file_; /**< Output stream for logging selected reactions */"

        // Array over metabolites
        // Has the metabolite been produced by a reaction at
        // a previous level
        std::vector<bool> available_; /**< Tracks whether each metabolite has been produced at a previous level */"

        // Summary stats
        std::string summary_; /**< Summary statistics string for the solving process */"
        size_t cost_iters_; /**< Number of cost iterations performed */"
        size_t level_iters_; /**< Number of level iterations performed */"
        size_t cycle_mets_; /**< Number of cycle metabolites processed */"
    };
}
