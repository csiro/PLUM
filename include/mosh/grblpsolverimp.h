/**
 * @file grblpsolverimp.h
 * @brief Gurobi-based LP solver implementation for metabolic network optimization
 *
 * This file provides the Gurobi solver implementation for flux balance analysis
 * and metabolic gap-filling problems. It implements the LPSolverImp interface
 * using the Gurobi optimization library.
 *
 * Compile with: g++ fn.cpp -lgurobi81 -lgurobi_g++5.2
 */

#pragma once

#include <vector>
#include <memory>

#include "gurobi_c++.h"

#include "mosh/lpsolverimp.h"
#include "mosh/dualvals.h"

#include "lime/numutil.h"

/**
 * @brief Metabolic optimization and solver handling namespace
 */
namespace mosh
{
    /** @brief Shared pointer to Gurobi environment for resource management */
    using GrbEnvPtr = std::shared_ptr<GRBEnv>;

    /**
     * @class GrbLPSolverImp
     * @brief Gurobi-based linear programming solver implementation
     *
     * Implements the LPSolverImp interface using the Gurobi optimization library
     * to solve flux balance analysis and metabolic gap-filling problems. This class
     * manages decision variables for reaction fluxes and usage indicators, along with
     * metabolic constraints and optimization objectives.
     *
     * The solver supports both continuous (FBA) and integer (gap-filling) formulations,
     * with multi-experiment scenarios for simultaneous optimization across different
     * experimental conditions.
     */
    class GrbLPSolverImp : public LPSolverImp
    {
    public:
        /**
         * @brief Constructs a Gurobi LP solver instance
         *
         * Initializes the Gurobi model with solver parameters including thread count,
         * random seed, time limits, MIP gap tolerance, and feasibility tolerances.
         * Allocates variable storage for reaction fluxes across all experiments.
         *
         * @param scenario Pointer to the metabolic scenario containing reactions and metabolites
         * @param params Pointer to solver parameters (threads, tolerances, time limits)
         * @param seed Random seed for reproducible optimization
         * @param env Shared pointer to Gurobi environment
         */
        GrbLPSolverImp(
            Scenario* scenario, Params* params, int seed, GrbEnvPtr env
        ) :
            LPSolverImp(scenario, params),
            model_(*env),
            cts_sol_(nullptr),
            flux_(
                scenario_->num_reactions(),
                std::vector<GRBVar>(scenario->num_experiments())
            ),
            use_(scenario_->num_reactions()),
            constr_(),
            biomass_(scenario->num_experiments())
        {
            model_.set(GRB_StringAttr_ModelName, "plum_lp");
            model_.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
            model_.set(GRB_IntParam_Threads, params->num_threads);
            model_.set(GRB_IntParam_Seed, seed);
            if (params->time_limit > 0) {
                double dbl_limit = (double) params->time_limit;
                model_.set(GRB_DoubleParam_TimeLimit, dbl_limit);
            }
            model_.set(GRB_DoubleParam_MIPGap, params_->max_mip_gap);
            model_.set(GRB_DoubleParam_IntFeasTol, params_->int_feas_tol);
            model_.set(GRB_DoubleParam_FeasibilityTol, params_->feasibility_tol);
            //limeSetEpsilon(1e-8);
        }

        /**
         * @brief Initializes the continuous (LP relaxation) formulation
         *
         * Sets up the solver for continuous flux balance analysis without
         * integer constraints. Currently has no implementation body.
         */
        void init_cts() override
        {
        }

        /**
         * @brief Initializes the integer (gap-filling) formulation
         *
         * Sets up the solver for mixed-integer programming with reaction usage indicators.
         * Stores the continuous solution for reference during integer optimization.
         *
         * @param cts_sol Shared pointer to the continuous solution from LP relaxation
         */
        void init_int(SolutionPtr cts_sol) override
        {
            cts_sol_ = cts_sol;
            // Switch off presolve... (for debugging)
            //model_.set(GRB_IntParam_Presolve, 0);
        }

        /**
         * @brief Creates a decision variable for reaction flux
         *
         * Adds a continuous variable representing the flux (flow rate) through a
         * metabolic reaction in a specific experiment.
         *
         * @param react Pointer to the reaction
         * @param exp Experiment index
         * @param name Variable name for debugging/output
         * @param obj_coeff Objective function coefficient (cost/weight)
         * @param lb Lower bound on flux value
         * @param ub Upper bound on flux value
         */
        void make_flux_var (
            const Reaction* react, size_t exp, std::string name,
            double obj_coeff, double lb, double ub
        ) override;
        /**
         * @brief Creates a binary decision variable for reaction usage
         *
         * Adds an integer variable (typically binary) indicating whether a reaction
         * is used in the gap-filled network. Used in MIP formulations for minimizing
         * the number of added reactions.
         *
         * @param react Pointer to the reaction
         * @param name Variable name for debugging/output
         * @param obj_coeff Objective function coefficient (penalty for adding reaction)
         * @param lb Lower bound (typically 0)
         * @param ub Upper bound (typically 1 for binary)
         */
        void make_use_var (
            const Reaction* react, std::string name, double obj_coeff, double lb, double ub
        ) override;
        /**
         * @brief Adds a metabolite mass balance constraint
         *
         * Enforces steady-state mass balance for a metabolite in a specific experiment.
         * The constraint ensures that the sum of production minus consumption equals zero
         * (or falls within specified bounds).
         *
         * @param met Pointer to the metabolite
         * @param experiment Experiment index
         * @param lb Lower bound on net production/consumption
         * @param ub Upper bound on net production/consumption
         * @param reacts Vector of reactions involving this metabolite
         * @param coeffs Stoichiometric coefficients for each reaction
         */
        void add_met_constraint (
            const Metabolite* met, size_t experiment, double lb, double ub,
            std::vector<Reaction*>& reacts, std::vector<double>& coeffs
        ) override;
        /**
         * @brief Adds a constraint linking flux and usage variables
         *
         * Creates a constraint that couples the continuous flux variable with the
         * binary usage indicator, typically enforcing that flux can only be nonzero
         * if the usage variable is 1.
         *
         * @param react Pointer to the reaction
         * @param exp Experiment index
         */
        void add_react_link_constraint (
            const Reaction* react, size_t exp
        ) override;
        /**
         * @brief Updates the bounds on a reaction flux variable
         *
         * Modifies the lower and upper bounds of an existing flux variable,
         * useful for incorporating experimental constraints or thermodynamic
         * irreversibility.
         *
         * @param react Pointer to the reaction
         * @param exp Experiment index
         * @param lb New lower bound on flux
         * @param ub New upper bound on flux
         */
        void set_react_bounds (
            const Reaction* react, size_t exp, double lb, double ub
        ) override
        {
            DEBUG (
                'z', "      Set bounds for react " << *react <<
                " exp " << exp <<
                " to [" << lb << "," << ub << "]"
            );
            flux_[react->index()][exp].set (GRB_DoubleAttr_LB, lb);
            flux_[react->index()][exp].set (GRB_DoubleAttr_UB, ub);
        }
        /**
         * @brief Updates the objective coefficient for a reaction flux
         *
         * Changes the cost/weight of a flux variable in the objective function,
         * allowing dynamic adjustment of optimization priorities.
         *
         * @param exp Experiment index
         * @param react Pointer to the reaction
         * @param cost New objective coefficient
         */
        void set_react_cost (
            size_t exp, const Reaction* react, double cost
        ) override
        {
            DEBUG (
                'z', "      Set cost for react " << *react <<
                " to " << cost
            );
            flux_[react->index()][exp].set (GRB_DoubleAttr_Obj, cost);
        }
        /**
         * @brief Sets the objective multiplier for biomass production
         *
         * Adjusts the coefficient of the biomass flux variable in the objective,
         * typically used to maximize or ensure minimum biomass production.
         *
         * @param exp Experiment index
         * @param mult Multiplier value (negative to maximize biomass)
         */
        void set_biomass_mult (size_t exp, double mult) override
        {
            biomass_[exp].set (GRB_DoubleAttr_Obj, mult);
        }

        /**
         * @brief Finalizes the optimization model before solving
         *
         * Updates the Gurobi model to incorporate all variables and constraints
         * that have been added, preparing it for optimization.
         */
        void finalise_formulation() override
        {
            model_.update();
        }

        /**
         * @brief Executes the optimization solver
         *
         * Runs the Gurobi optimizer on the formulated model and returns the
         * solution status (optimal, infeasible, timeout, etc.).
         *
         * @return StatusEnum indicating the optimization result
         */
        StatusEnum optimize() override;

        /**
         * @brief Extracts solution for a specific experiment
         *
         * Retrieves the optimized flux values from the Gurobi solution and
         * packages them into a Solution object.
         *
         * @param which_experiment Index of the experiment to extract
         * @return SolutionPtr shared pointer to the extracted solution
         */
        SolutionPtr  make_sol (size_t which_experiment) override;
        /**
         * @brief Retrieves the optimal objective value
         *
         * @return The objective function value at the optimal solution
         */
        double get_objective () override
        {
            return model_.get(GRB_DoubleAttr_ObjVal);
        }
        /**
         * @brief Retrieves the MIP optimality gap
         *
         * Returns the relative gap between the best integer solution and the
         * best bound, indicating solution quality for integer programs.
         *
         * @return MIP gap as a fraction (0.0 = proven optimal)
         */
        double get_mip_gap () override
        {
            return model_.get(GRB_DoubleAttr_MIPGap);
        }

        /**
         * @brief Writes the optimization model to a file
         *
         * Exports the model in Gurobi's native format for inspection or debugging.
         *
         * @param filename Output file path (format determined by extension: .lp, .mps, etc.)
         */
        void write_model (std::string filename) override
        {
            model_.write (filename);
        }
        /**
         * @brief Returns the number of decision variables in the model
         *
         * @return Total count of variables (continuous and integer)
         */
        size_t num_vars() override
        {
            return (size_t)model_.get (GRB_IntAttr_NumVars);
        }
        /**
         * @brief Returns the number of constraints in the model
         *
         * @return Total count of constraints (equalities and inequalities)
         */
        size_t num_constraints() override
        {
            return (size_t)model_.get (GRB_IntAttr_NumConstrs);
        }

        /**
         * @brief Accesses flux variable for single-experiment scenario
         *
         * Convenience accessor for the flux variable of a reaction when only
         * one experiment is present (experiment 0).
         *
         * @param k Reaction index
         * @return Reference to the Gurobi variable for the reaction flux
         */
        GRBVar& flux(size_t k) {return flux_[k][0];}
        /**
         * @brief Accesses flux variable for multi-experiment scenario
         *
         * @param k Reaction index
         * @param exp Experiment index
         * @return Reference to the Gurobi variable for the reaction flux
         */
        GRBVar& flux(size_t k, size_t exp) {return flux_[k][exp];}

        /**
         * @brief Controls solver output verbosity
         *
         * Enables or disables Gurobi console output during optimization.
         *
         * @param quiet If true, suppresses solver output; if false, enables output
         */
        void local_set_quiet (bool quiet) override
        {
            model_.set(GRB_IntParam_OutputFlag, quiet ? 0 : 1);
        }
        
        /**
         * @brief Selects a subset of high-quality reactions for gap-filling
         *
         * Uses heuristics or optimization to identify the most promising reactions
         * to add to the metabolic network.
         *
         * @param num_to_select Number of reactions to select
         * @return Count of successfully selected reactions
         */
        int select_good_reactions (int num_to_select);
        /**
         * @brief Retrieves dual values (shadow prices) from the solution
         *
         * Extracts dual values for constraints, which represent the marginal value
         * of relaxing each constraint. Useful for sensitivity analysis and identifying
         * metabolic bottlenecks.
         *
         * @return DualValsPtr shared pointer to dual values structure
         */
        DualValsPtr get_dual_vals();
        /**
         * @brief Computes the reduced cost for a reaction
         *
         * Returns the reduced cost of a reaction's flux variable, indicating how much
         * the objective would worsen if the reaction were forced to have nonzero flux.
         *
         * @param react Pointer to the reaction
         * @return Reduced cost value
         */
        double reduced_cost (const Reaction* react);

    protected:
        /** @brief Gurobi optimization model containing all variables and constraints */
        GRBModel model_;

        /** @brief Continuous solution from LP relaxation, used as reference for MIP */
        SolutionPtr cts_sol_;

        /** @brief Decision variables for reaction fluxes [reaction_index][experiment_index] */
        std::vector<std::vector<GRBVar>> flux_;
        /** @brief Binary decision variables indicating reaction usage in gap-filling */
        std::vector<GRBVar> use_;
        /** @brief Gurobi constraint objects for programmatic constraint manipulation */
        std::vector<GRBConstr> constr_;

        /** @brief Decision variables for biomass production flux in each experiment */
        std::vector<GRBVar> biomass_;
    };
    /** @brief Shared pointer to GrbLPSolverImp for resource management */
    using GrbLPSolverImpPtr = std::shared_ptr<GrbLPSolverImp>;
}
