/**
 * @file lpsolverimp.h
 * @brief Abstract base class for LP solver implementations in metabolic gap-filling
 *
 * This file defines the LPSolverImp interface, which provides an abstraction layer
 * for different linear programming solver backends (e.g., CPLEX, Gurobi, GLPK).
 * The interface supports both continuous (LP) and integer (MILP) formulations for
 * metabolic network gap-filling and flux balance analysis problems.
 */

#pragma once

/* An abstraction of the implementation of an LP solver
 */

#include "mosh/gapsolver.h"
#include "mosh/dualvals.h"

/**
 * @namespace mosh
 * @brief Main namespace for the MOSH (Metabolic Optimization and gap-filling Solver Heuristics) library
 */
namespace mosh
{
    /**
     * @class LPSolverImp
     * @brief Abstract base class providing interface for LP/MILP solver implementations
     *
     * This class defines the interface that must be implemented by concrete LP solver
     * backends. It handles formulation of metabolic network optimization problems including:
     * - Flux variables for reactions across experiments
     * - Binary use variables for reaction selection
     * - Metabolite mass balance constraints
     * - Reaction-use variable linking constraints
     *
     * The solver supports multi-experiment scenarios where the same network is optimized
     * under different experimental conditions simultaneously.
     */
    class LPSolverImp
    {
    public:
        /**
         * @brief Constructor for LP solver implementation
         * @param scenario Pointer to the metabolic scenario containing network and experiments
         * @param params Pointer to solver parameters and configuration
         */
        LPSolverImp (Scenario* scenario, Params* params) :
            scenario_(scenario),
            params_(params)
        {
        }

        /**
         * @brief Initialize the continuous (LP) relaxation of the problem
         *
         * Sets up the LP formulation without integer constraints. Used for the initial
         * continuous relaxation phase before solving the MILP.
         */
        virtual void init_cts() = 0;
        /**
         * @brief Initialize the integer (MILP) formulation
         * @param cts_sol Solution from the continuous relaxation, used as warmstart or for branching priorities
         */
        virtual void init_int(SolutionPtr cts_sol) = 0;

        /**
         * @brief Create a flux variable for a reaction in a specific experiment
         * @param react Pointer to the reaction for which to create a flux variable
         * @param exp Index of the experiment (0-based)
         * @param name Variable name for the solver
         * @param obj_coeff Objective function coefficient for this variable
         * @param lb Lower bound on flux value
         * @param ub Upper bound on flux value
         */
        virtual void make_flux_var (
            const Reaction* react, size_t exp, std::string name,
            double obj_coeff, double lb, double ub
        ) = 0;
        /**
         * @brief Create flux variables for a reaction across all experiments
         * @param react Pointer to the reaction for which to create flux variables
         * @param name Base variable name (will be suffixed with experiment index)
         * @param obj_coeff Objective function coefficient for these variables
         * @param lb Lower bound on flux values
         * @param ub Upper bound on flux values
         */
        void make_flux_vars (
            const Reaction* react, std::string name,
            double obj_coeff, double lb, double ub
        ) {
            for (size_t exp = 0; exp < num_experiments(); exp++)
                make_flux_var (react, exp, name, obj_coeff, lb, ub);
        }
        /**
         * @brief Create a binary use variable for a reaction
         * @param react Pointer to the reaction for which to create a use variable
         * @param name Variable name for the solver
         * @param obj_coeff Objective function coefficient (typically a penalty cost)
         * @param lb Lower bound (typically 0)
         * @param ub Upper bound (typically 1 for binary variables)
         *
         * Use variables indicate whether a reaction is present/active in the gap-filled network.
         * They are shared across all experiments (a reaction is either added or not).
         */
        virtual void make_use_var (
            const Reaction* react, std::string name, double obj_coeff,
            double lb, double ub
        ) = 0;
        /**
         * @brief Add a metabolite mass balance constraint for steady-state flux
         * @param met Pointer to the metabolite for which to enforce mass balance
         * @param experiment Index of the experiment
         * @param lb Lower bound for constraint (typically 0 for steady state)
         * @param ub Upper bound for constraint (typically 0 for steady state)
         * @param reacts Vector of reactions involving this metabolite
         * @param coeffs Stoichiometric coefficients corresponding to each reaction
         *
         * Enforces: lb <= sum(coeffs[i] * flux[reacts[i]]) <= ub
         */
        virtual void add_met_constraint (
            const Metabolite* met, size_t experiment, double lb, double ub,
            std::vector<Reaction*>& reacts, std::vector<double>& coeffs
        ) = 0;
        /**
         * @brief Add constraint linking reaction flux to its use variable for a specific experiment
         * @param react Pointer to the reaction
         * @param exp Index of the experiment
         *
         * Enforces that flux can only be non-zero if the use variable is 1, typically:
         * lb * use <= flux <= ub * use
         */
        virtual void add_react_link_constraint (
            const Reaction* react, size_t exp
        ) = 0;
        /**
         * @brief Add flux-use linking constraints for a reaction across all experiments
         * @param react Pointer to the reaction
         */
        void add_react_link_constraints (const Reaction* react)
        {
            for (size_t exp = 0; exp < num_experiments(); exp++) {
                add_react_link_constraint (react, exp);
            }
        }
        /**
         * @brief Update the bounds on a reaction flux variable
         * @param react Pointer to the reaction
         * @param exp Index of the experiment
         * @param lb New lower bound on flux
         * @param ub New upper bound on flux
         */
        virtual void set_react_bounds (
            const Reaction* react, size_t exp, double lb, double ub
        ) = 0;
        /**
         * @brief Set the objective function coefficient for a reaction flux
         * @param exp Index of the experiment
         * @param react Pointer to the reaction
         * @param cost New objective coefficient value
         */
        virtual void set_react_cost (
            size_t exp, const Reaction* react, double cost
        ) = 0;
        /**
         * @brief Set the biomass production multiplier for an experiment
         * @param exp Index of the experiment
         * @param mult Multiplier value (e.g., for enforcing minimum growth rates)
         */
        virtual void set_biomass_mult (size_t exp, double mult) = 0;
        /**
         * @brief Finalize the problem formulation before optimization
         *
         * Optional hook for solver-specific setup after all variables and constraints
         * have been added. Default implementation does nothing.
         */
        virtual void finalise_formulation() {};
        
        /**
         * @brief Solve the optimization problem
         * @return Status code indicating success, infeasibility, timeout, etc.
         */
        virtual StatusEnum optimize() = 0;

        /**
         * @brief Extract solution for all experiments from the solved model
         * @return Shared pointer to Solution object (MultiSol if multiple experiments, otherwise single-experiment solution)
         */
        virtual SolutionPtr make_sol()
        {
            if (num_experiments() > 1) {
                MultiSolPtr sol =
                    std::make_shared<MultiSol> (scenario_, params_);
                for (size_t k = 0; k < num_experiments(); k++)
                    sol->set_sol (k, make_sol(k));
                return sol;
            }
            return make_sol (0);
        }
        
        /**
         * @brief Extract solution for a specific experiment from the solved model
         * @param which_experiment Index of the experiment to extract
         * @return Shared pointer to Solution object for the specified experiment
         */
        virtual SolutionPtr make_sol(size_t which_experiment) = 0;
        /**
         * @brief Get the objective function value of the optimal solution
         * @return Objective value
         */
        virtual double get_objective () = 0;
        /**
         * @brief Get the MIP optimality gap
         * @return Gap value (typically as fraction: (upper_bound - lower_bound) / upper_bound)
         *
         * Returns 0 for LP problems or when optimal solution is proven.
         */
        virtual double get_mip_gap () = 0;
        
        /**
         * @brief Write the optimization model to a file for debugging
         * @param filename Output file path (format depends on solver: .lp, .mps, etc.)
         */
        virtual void write_model (std::string filename) = 0;

        /**
         * @brief Get the total number of variables in the formulation
         * @return Number of variables
         */
        virtual size_t num_vars() = 0;
        /**
         * @brief Get the total number of constraints in the formulation
         * @return Number of constraints
         */
        virtual size_t num_constraints() = 0;

        /**
         * @brief Get the number of metabolites in the scenario
         * @return Number of metabolites
         */
        size_t num_metabolites() const {return scenario_->num_metabolites();}
        /**
         * @brief Access a metabolite by index
         * @param k Metabolite index (0-based)
         * @return Pointer to the Metabolite object
         */
        const Metabolite* metabolite(size_t k)
        {
            return scenario_->metabolite(k);
        }
        
        /**
         * @brief Get the number of reactions in the scenario
         * @return Number of reactions
         */
        size_t num_reactions() const {return scenario_->num_reactions();}
        /**
         * @brief Access a reaction by index
         * @param k Reaction index (0-based)
         * @return Pointer to the Reaction object
         */
        Reaction* reaction(size_t k) const {return scenario_->reaction(k);}
        
        /**
         * @brief Get the number of experiments in the scenario
         * @return Number of experiments
         */
        size_t num_experiments() const {return scenario_->num_experiments();}
        /**
         * @brief Access an experiment by index
         * @param k Experiment index (0-based)
         * @return Pointer to the Experiment object
         */
        const Experiment* experiment(size_t k) const
        {
            return scenario_->experiment(k);
        }
        
        /**
         * @brief Control solver output verbosity
         * @param quiet If true, suppress solver output; if false, enable verbose output
         *
         * Default implementation does nothing. Concrete solvers should override to
         * configure their logging/output behavior.
         */
        virtual void local_set_quiet (bool quiet)
        {
        }
        
    protected:
        /** @brief Pointer to the metabolic scenario (network, experiments, etc.) */
        Scenario* scenario_;
        /** @brief Pointer to solver parameters and configuration */
        Params* params_;        
    };
    /** @brief Shared pointer type for LPSolverImp objects */
    using LPSolverImpPtr = std::shared_ptr<LPSolverImp>;
}
