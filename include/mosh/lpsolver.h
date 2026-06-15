/**
 * @file lpsolver.h
 * @brief Linear programming solver for metabolic gap-filling problems
 *
 * This file defines the LPSolver class, which implements a linear programming
 * approach to solving metabolic network gap-filling problems. It extends the
 * GapSolver interface and uses an LP solver implementation to find optimal
 * reaction sets that enable target metabolite production or biomass flux.
 */
#pragma once

#include "mosh/gapsolver.h"
#include "mosh/lpsolverimp.h"
#include "mosh/dualvals.h"

#include "lime/debug.h"

/**
 * @brief Metabolic optimization and solver heuristics namespace
 */
namespace mosh
{
    /**
     * @class LPSolver
     * @brief Linear programming solver for metabolic gap-filling
     *
     * LPSolver implements a linear programming approach to solve metabolic network
     * gap-filling problems. It formulates the problem as an LP optimization that
     * finds minimal sets of reactions to add to enable target metabolic functions.
     * The solver supports multiple formulations and can optimize biomass production
     * or target metabolite flux.
     */
    class LPSolver : public GapSolver
    {
    public:
        /**
         * @brief Construct an LPSolver instance
         * @param scenario The metabolic scenario containing reactions and metabolites
         * @param params Solver parameters controlling optimization behavior
         * @param seed Random seed for stochastic elements
         * @param which_formulation The LP formulation type to use
         * @param imp Implementation pointer to the underlying LP solver
         */
        LPSolver(
            Scenario* scenario, Params* params, int seed,
            int which_formulation, 
            LPSolverImpPtr imp
        ) :
            GapSolver(scenario, params, seed),
            which_experiment_(0), 
            which_formulation_(which_formulation),
            imp_(imp),
            target_flux_(0.0f),
            biomass_search_iters_(0),
            biomass_obj_mult_(params->init_biomass_obj_mult)
        {
        }

        /**
         * @brief Get the current experiment index
         * @return The index of the experiment being solved
         */
        size_t which_experiment() const {return which_experiment_;}
        /**
         * @brief Set the current experiment to solve
         * @param which_experiment Index of the experiment to use
         */
        void set_which_experiment (size_t which_experiment) {
            which_experiment_ = which_experiment;
        }
        /**
         * @brief Get the current LP formulation type
         * @return The formulation index being used
         */
        size_t which_formulation() const {return which_formulation_;}
        /**
         * @brief Set the LP formulation type to use
         * @param which_formulation Index specifying the formulation
         */
        void set_which_formulation (size_t which_formulation);

        /**
         * @brief Solve the gap-filling problem
         * @return Pointer to the solution containing added reactions and flux values
         */
        SolutionPtr solve() override;

        /**
         * @brief Formulate the LP problem for a given experiment
         * @param exp Pointer to the experiment defining constraints and objectives
         */
        // ALternatively, solve via formulate,
        void formulate (const Experiment* exp);
        /**
         * @brief Solve the previously formulated LP problem
         * @return Pointer to the solution with optimal reaction set and fluxes
         */
        // then re-solve after adding reactions
        SolutionPtr do_solve();
        
        /**
         * @brief Enable a reaction in the LP formulation
         * @param react Pointer to the reaction to enable
         */
        // Enable a reaction
        void enable_reaction (const Reaction* react);
        /**
         * @brief Disable a reaction in the LP formulation
         * @param react Pointer to the reaction to disable
         */
        // Disable a reaction
        void disable_reaction (const Reaction* react);

        /**
         * @brief Set the cost coefficient for all non-biomass reactions
         * @param exp Experiment index
         * @param cost Cost value to assign to each reaction
         */
        void set_all_react_cost (size_t exp, double cost)
        {
            DEBUG (
                'l', "Set all react cost for exp " << exp << " to " << cost
            );
            for (auto& react_ptr : scenario_->reactions()) {
                Reaction* react = react_ptr.get();
                if (!react->is_biomass())
                    imp_->set_react_cost (exp, react, cost);
            }
        }
        /**
         * @brief Set the cost coefficient for a specific reaction
         * @param exp Experiment index
         * @param react Pointer to the reaction
         * @param cost Cost value to assign
         */
        void set_react_cost (size_t exp, Reaction* react, double cost)
        {
            DEBUG (
                'l', "Set react cost for exp " << exp <<
                " react " << react->name() << " to " << cost
            );
            assert (!react->is_biomass());
            imp_->set_react_cost (exp, react, cost);
        }

        /**
         * @brief Get the current biomass objective multiplier
         * @return The multiplier applied to the biomass objective term
         */
        double biomass_obj_mult() const {return biomass_obj_mult_;}
        /**
         * @brief Set the biomass objective multiplier
         * @param exp Experiment index
         * @param mult Multiplier value (sign will be negated internally)
         * @note The sign is negated when passed to the implementation
         */
        // Set obj multiplier. NOTE : Negates the sign
        void set_biomass_obj_mult (size_t exp, double mult)
        {
            DEBUG (
                'l', "Set biomass obj mult for exp " << exp << " to " << mult <<
                " was " << biomass_obj_mult_
            );
            biomass_obj_mult_ = mult;
            imp_->set_biomass_mult (exp, -mult);
        }
        /**
         * @brief Get the target flux value
         * @return The target flux constraint value
         */
        double target_flux() const {return target_flux_;}
        /**
         * @brief Set the target flux constraint value
         * @param target_flux The flux value to target
         */
        void set_target_flux (double target_flux) {
            target_flux_ = target_flux;
        }

        /**
         * @brief Write the LP model to a file
         * @param fn Filename for the output model file
         */
        void write_model (std::string fn)
        {
            imp_->write_model (fn);
        }

        /**
         * @brief Generate a summary string of the solver state and results
         * @return String containing solver summary information
         */
        std::string summary() override;

        /**
         * @brief Set quiet mode to suppress solver output
         * @param quiet True to suppress output, false for verbose mode
         */
        void set_quiet (bool quiet) override
        {
            GapSolver::set_quiet (quiet);
            imp_->local_set_quiet (quiet);
        }
        

    protected:
        /**
         * @brief Perform biomass objective search to find optimal biomass multiplier
         * @param sol Initial solution to refine
         * @return Pointer to solution with optimized biomass objective
         */
        SolutionPtr biomass_obj_search (SolutionPtr sol);

        size_t which_experiment_; /**< Index of the current experiment being solved */
        size_t which_formulation_; /**< Type of LP formulation being used */
        LPSolverImpPtr imp_; /**< Pointer to the underlying LP solver implementation */

        // Target flux for forumation 2
        double target_flux_; /**< Target flux constraint value for formulation 2 */
        
        // Solution attributes
        int biomass_search_iters_; /**< Number of iterations in biomass objective search */
        double biomass_obj_mult_; /**< Multiplier applied to biomass objective term */
        static int num_model_writes_; /**< Static counter for the number of model writes performed */
    };
    /**
     * @typedef LPSolverPtr
     * @brief Shared pointer type for LPSolver instances
     */
    using LPSolverPtr = std::shared_ptr<LPSolver>;
}
