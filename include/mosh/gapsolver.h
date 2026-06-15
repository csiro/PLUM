/**
 * @file gapsolver.h
 * @brief Abstract base class for metabolic gap-filling solvers
 *
 * This file defines the GapSolver interface for solving metabolic gap-filling
 * problems using flux balance analysis (FBA). Gap-filling algorithms identify
 * missing reactions that enable a metabolic network to produce target metabolites.
 */
#pragma once

#include <vector>
#include <memory>

#include "lime/dig.h"
#include "lime/timekeeper.h"
#include "lime/constants.h"
#include "lime/rand.h"
#include "lime/debug.h"
#include "lime/strutil.h"

#include "mosh/scenario.h"
#include "mosh/params.h"
#include "mosh/solution.h"
#include "mosh/multisol.h"
#include "mosh/reaction.h"
#include "mosh/metabolite.h"

/**
 * @namespace mosh
 * @brief Metabolic Optimization and Simulation Hub namespace
 */
namespace mosh
{
    /**
     * @class GapSolver
     * @brief Abstract base class for metabolic gap-filling solvers
     *
     * GapSolver provides a common interface for different gap-filling algorithms
     * that identify minimal sets of reactions needed to complete metabolic pathways.
     * It manages the optimization scenario, parameters, timing, and random number
     * generation for derived solver implementations.
     *
     * Derived classes must implement the solve() method to provide specific
     * gap-filling algorithms (e.g., MILP-based, heuristic-based approaches).
     */
    class GapSolver 
    {
    public:
        /**
         * @brief Constructs a GapSolver with scenario, parameters, and random seed
         * @param scenario Pointer to the metabolic scenario containing network and experiments
         * @param params Pointer to solver parameters (time limits, tolerances, etc.)
         * @param seed Random seed for reproducible stochastic behavior
         */
        GapSolver(
            Scenario* scenario, Params* params, int seed
        ) :
            scenario_(scenario),
            params_(params),
            quiet_(false),
            rand_(seed),
            status_(UNSOLVED),
            best_sol_fn_(""),
            timer_()
        {
            if (params_->time_limit > 0)
                timer_.setTimeLimit (params_->time_limit);
        }

        /**
         * @brief Returns the number of metabolites in the scenario
         * @return Number of metabolites in the metabolic network
         */
        size_t num_metabolites() const {return scenario_->num_metabolites();}
        /**
         * @brief Retrieves a metabolite by index
         * @param k Zero-based index of the metabolite
         * @return Pointer to the metabolite at index k
         */
        const Metabolite* metabolite(size_t k)
        {
            return scenario_->metabolite(k);
        }
        
        /**
         * @brief Returns the number of reactions in the scenario
         * @return Number of reactions in the metabolic network
         */
        size_t num_reactions() const {return scenario_->num_reactions();}
        /**
         * @brief Retrieves a reaction by index
         * @param k Zero-based index of the reaction
         * @return Pointer to the reaction at index k
         */
        Reaction* reaction(size_t k) const {return scenario_->reaction(k);}
        
        /**
         * @brief Returns the number of experimental conditions in the scenario
         * @return Number of experiments (growth conditions, media, etc.)
         */
        size_t num_experiments() const {return scenario_->num_experiments();}
        /**
         * @brief Retrieves an experiment by index
         * @param k Zero-based index of the experiment
         * @return Pointer to the experiment at index k
         */
        Experiment* experiment(size_t k) const
        {
            return scenario_->experiment(k);
        }
        
        /**
         * @brief Gets the filename of the best solution found
         * @return Filename where the best solution is/will be saved
         */
        std::string best_sol_fn() const {return best_sol_fn_;}
        /**
         * @brief Sets the filename for saving the best solution
         * @param best_sol_fn Filename path for the best solution output
         */
        void set_best_sol_fn (std::string best_sol_fn) {
            best_sol_fn_ = best_sol_fn;
        }

        /**
         * @brief Generates a random seed for Gurobi optimizer
         * @return Random integer seed in range [0, 1000000)
         */
        int grb_seed() {return rand_.uniform0n_1 (1000000);}
        /**
         * @brief Generates a general-purpose random seed
         * @return Random integer seed for stochastic algorithms
         */
        int rand_seed() {return rand_.generateSeed();}

        /**
         * @brief Pure virtual method to solve the gap-filling problem
         *
         * Derived classes implement specific gap-filling algorithms to find
         * minimal sets of reactions that satisfy the experimental constraints.
         *
         * @return Shared pointer to the solution containing selected reactions and fluxes
         */
        virtual SolutionPtr solve() = 0;

        /**
         * @brief Returns a human-readable summary of the solver state
         * @return String summarizing solver configuration and results
         */
        virtual std::string summary() {return "";}

        /**
         * @brief Gets elapsed solving time in seconds
         * @return Number of seconds since solver started
         */
        double elapsed_time_secs () const {return timer_.elapsedTimeSecs();}
        /**
         * @brief Checks if time limit has not been exceeded
         * @return True if solver still has time remaining, false otherwise
         */
        bool has_time_left() const {return timer_.hasTimeLeft();}
        /**
         * @brief Gets remaining time before time limit is reached
         * @return Number of seconds remaining (negative if limit exceeded)
         */
        double time_left_seconds() const {return timer_.timeLeftSecs();}
        /**
         * @brief Gets the current solution status
         * @return StatusEnum indicating solution state (optimal, infeasible, etc.)
         */
        StatusEnum status() const {return status_;}
        /**
         * @brief Converts solution status to human-readable string
         * @return String representation of the current status
         */
        std::string status_str() const
        {
            std::vector<std::string> stat_str =
                {
                    "unsolved",
                    "cts_optimal",
                    "int_optimal",
                    "suboptimal",
                    "infeasible",
                    "random"
                };
            return stat_str[status_];
        }

        /**
         * @brief Checks if solver is in quiet mode
         * @return True if console output is suppressed, false otherwise
         */
        bool quiet() const {return quiet_;}
        /**
         * @brief Sets quiet mode to suppress or enable console output
         * @param quiet True to suppress output, false to enable
         */
        virtual void set_quiet (bool quiet) {
            quiet_ = quiet;
        }
        /**
         * @brief Displays a message to console (respects quiet mode)
         * @param message Text to display with colored formatting
         */
        void on_screen (std::string message) {
            DEBUG ('A', message);
            if (!quiet_) {
                std::cout <<
                    escSeqColour (lime::CYAN) <<
                    message <<
                    escSeqColour (lime::RESET) <<
                    std::endl;
            }
        }

    protected:
        Scenario* scenario_; /**< Pointer to metabolic scenario (network, experiments) */"
        Params* params_; /**< Pointer to solver parameters (tolerances, limits) */"
        bool quiet_; /**< Flag to suppress console output */"

        lime::Rand rand_; /**< Random number generator for stochastic operations */"

        StatusEnum status_; /**< Current solution status (optimal, infeasible, etc.) */"

        std::string best_sol_fn_; /**< Filename for best solution output */"
        
    private:
        lime::TimeKeeper timer_; /**< Timer for tracking elapsed time and enforcing limits */"
    };

    /**
     * @typedef GapSolverPtr
     * @brief Shared pointer type for GapSolver instances
     */
    using GapSolverPtr = std::shared_ptr<GapSolver>;
}
