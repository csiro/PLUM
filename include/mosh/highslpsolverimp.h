/**
 * @file highslpsolverimp.h
 * @brief HiGHS-based implementation of the linear programming solver interface
 *
 * This file defines the HighsLPSolverImp class, which provides a concrete implementation
 * of the LPSolverImp interface using the HiGHS optimization solver. It is used for
 * solving mixed-integer linear programs in metabolic network gap-filling and flux
 * balance analysis problems.
 */
#pragma once

#include "lime/numutil.h"

#include "mosh/lpsolverimp.h"
#include "mosh/highsutil.h"

#include "highs/Highs.h"

/**
 * @namespace mosh
 * @brief Metabolic optimization solver harness namespace
 */
namespace mosh
{
    /**
     * @class HighsLPSolverImp
     * @brief HiGHS-based linear programming solver implementation
     *
     * This class implements the LPSolverImp interface using the HiGHS optimization library.
     * It manages the construction and solution of mixed-integer linear programs for metabolic
     * network analysis, including flux balance analysis and metabolic gap-filling problems.
     * The class handles variable creation, constraint addition, optimization execution, and
     * solution retrieval for multiple experimental conditions.
     */
    class HighsLPSolverImp : public LPSolverImp
    {
    public:
        /**
         * @brief Constructs a HiGHS LP solver instance with specified parameters
         *
         * Initializes the HiGHS model with solver parameters including time limits,
         * thread count, MIP gap tolerance, feasibility tolerances, and random seed.
         * Allocates internal data structures for managing reaction variables across
         * multiple experimental conditions.
         *
         * @param scenario Pointer to the metabolic network scenario containing reactions and experiments
         * @param params Pointer to solver parameters (time limits, tolerances, threads, etc.)
         * @param seed Random seed for reproducible solver behavior
         */
        HighsLPSolverImp(
            Scenario* scenario, Params* params, int seed
        ) :
            LPSolverImp(scenario, params),
            model_(),
            flux_var_id_(
                scenario_->num_reactions(),
                std::vector<HighsInt> (scenario_->num_experiments())
            ),
            use_var_id_(scenario_->num_reactions()),
            highs_coeffs_(new double[scenario_->num_reactions()]),
            col_id_(new HighsInt[scenario_->num_reactions()]),
            cts_sol_(nullptr),
            biomass_id_(scenario_->num_experiments(), -1)
        {
            if (params->time_limit > 0) {
                check_highs_status (
                    model_.setOptionValue("time_limit", params->time_limit),
                    "setOptionValue-time_limit"
                );
            }
            check_highs_status (
                model_.setOptionValue("threads", params->num_threads),
                "setOptionValue-threads"
            );
            check_highs_status (
                model_.setOptionValue("mip_rel_gap", params->max_mip_gap),
                "setOptionValue-mip_rel_gap"
            );
            check_highs_status (
                model_.setOptionValue(
                    "primal_feasibility_tolerance", params->feasibility_tol
                ),
                "setOptionValue-primal_feasibility_tolerance"
            );
            check_highs_status (
                model_.setOptionValue(
                    "mip_feasibility_tolerance", params->int_feas_tol
                ),
                "setOptionValue-mip_feasibility_tolerance"
            );
            check_highs_status (
                model_.setOptionValue(
                    "random_seed", seed
                ),
                "setOptionValue-random_seed"
            );
        }
        /**
         * @brief Destroys the HiGHS LP solver instance and frees allocated resources
         *
         * Deallocates working arrays used for constraint construction.
         */
        virtual ~HighsLPSolverImp()
        {
            delete [] highs_coeffs_;
            delete [] col_id_;
        }

        /**
         * @brief Initializes the solver for continuous optimization problems
         *
         * This method is called before solving continuous (non-integer) linear programs.
         * Currently a no-op for the HiGHS implementation.
         */
        void init_cts() override
        {
        }
        
        /**
         * @brief Initializes the solver for mixed-integer optimization problems
         *
         * Stores a reference to the continuous solution for use in constructing
         * the mixed-integer program, such as for variable fixing or warm-starting.
         *
         * @param cts_sol Shared pointer to the continuous solution from a prior solve
         */
        void init_int(SolutionPtr cts_sol) override
        {
            cts_sol_ = cts_sol;
        }

        /**
         * @brief Creates a flux variable for a reaction in a specific experiment
         *
         * Adds a decision variable representing the flux (reaction rate) through a
         * metabolic reaction under a particular experimental condition. The variable
         * is added to the optimization model with the specified bounds and objective coefficient.
         *
         * @param react Pointer to the reaction for which to create a flux variable
         * @param exp Experiment index for this flux variable
         * @param name Variable name for debugging and model export
         * @param obj_coeff Objective function coefficient for this variable
         * @param lb Lower bound on the flux value
         * @param ub Upper bound on the flux value
         */
        void make_flux_var (
            const Reaction* react, size_t exp, std::string name,
            double obj_coeff, double lb, double ub
        ) override;
        /**
         * @brief Creates a binary usage variable for a reaction
         *
         * Adds a binary decision variable indicating whether a reaction is used
         * (active) in the metabolic network. This is typically used in gap-filling
         * formulations to penalize the addition of reactions.
         *
         * @param react Pointer to the reaction for which to create a usage variable
         * @param name Variable name for debugging and model export
         * @param obj_coeff Objective function coefficient (typically a penalty cost)
         * @param lb Lower bound (typically 0 for binary variables)
         * @param ub Upper bound (typically 1 for binary variables)
         */
        void make_use_var (
            const Reaction* react, std::string name, double obj_coeff,
            double lb, double ub
        ) override;
        /**
         * @brief Adds a metabolite steady-state mass balance constraint
         *
         * Creates a constraint enforcing metabolite mass balance (steady-state assumption)
         * for a specific metabolite in a given experiment. The constraint is a linear
         * combination of reaction flux variables weighted by stoichiometric coefficients.
         *
         * @param met Pointer to the metabolite for which to add the constraint
         * @param experiment Experiment index for this constraint
         * @param lb Lower bound for the metabolite balance (typically 0 for steady-state)
         * @param ub Upper bound for the metabolite balance (typically 0 for steady-state)
         * @param reacts Vector of reaction pointers participating in the constraint
         * @param coeffs Vector of stoichiometric coefficients corresponding to reactions
         */
        void add_met_constraint (
            const Metabolite* met, size_t experiment, double lb, double ub,
            std::vector<Reaction*>& reacts, std::vector<double>& coeffs
        ) override;
        /**
         * @brief Adds a constraint linking reaction flux to its usage indicator
         *
         * Creates constraints that activate the binary usage variable when the
         * corresponding flux variable is non-zero. This is used in MILP formulations
         * for metabolic gap-filling to track which reactions are active.
         *
         * @param react Pointer to the reaction to link
         * @param exp Experiment index for the flux variable to link
         */
        void add_react_link_constraint (
            const Reaction* react, size_t exp
        ) override;
        /**
         * @brief Updates the flux bounds for a reaction in a specific experiment
         *
         * Modifies the lower and upper bounds on a reaction's flux variable,
         * typically used to enforce experimental observations or thermodynamic constraints.
         *
         * @param react Pointer to the reaction whose bounds to update
         * @param exp Experiment index for the flux variable
         * @param lb New lower bound on the flux
         * @param ub New upper bound on the flux
         */
        void set_react_bounds (
            const Reaction* react, size_t exp, double lb, double ub
        ) override;
        /**
         * @brief Sets the objective coefficient for a reaction's flux variable
         *
         * Updates the cost (objective function coefficient) associated with a
         * reaction's flux in a specific experiment, used to penalize or favor
         * certain reaction activities.
         *
         * @param exp Experiment index for the flux variable
         * @param react Pointer to the reaction whose cost to set
         * @param cost New objective coefficient value
         */
        void set_react_cost (
            size_t exp, const Reaction* react, double cost
        ) override;
        /**
         * @brief Sets the biomass production multiplier for an experiment
         *
         * Adjusts the objective coefficient for the biomass flux variable in a
         * specific experiment, typically used to weight biomass production in
         * multi-objective or multi-experiment formulations.
         *
         * @param exp Experiment index
         * @param mult Multiplier to apply to the biomass objective coefficient
         */
        void set_biomass_mult (size_t exp, double mult) override;
        
        /**
         * @brief Executes the optimization and returns the solution status
         *
         * Solves the linear or mixed-integer program constructed through prior
         * variable and constraint additions. Returns a status code indicating
         * optimality, infeasibility, or solver errors.
         *
         * @return StatusEnum indicating the optimization result (optimal, infeasible, error, etc.)
         */
        StatusEnum optimize() override;

        /**
         * @brief Extracts the solution for a specific experiment
         *
         * Retrieves the optimal flux values for all reactions in the specified
         * experiment from the solved model and packages them into a Solution object.
         *
         * @param which_experiment Experiment index for which to extract the solution
         * @return SolutionPtr shared pointer to the solution object containing flux values
         */
        SolutionPtr make_sol(size_t which_experiment) override;
        /**
         * @brief Returns the optimal objective function value
         *
         * Retrieves the objective value from the most recent optimization solve.
         *
         * @return Optimal objective function value
         */
        double get_objective () override
        {
            return model_.getObjectiveValue();
        }
        /**
         * @brief Returns the mixed-integer programming optimality gap
         *
         * For MIP problems, returns the relative gap between the best solution found
         * and the best possible objective value. Returns 0 for proven optimal solutions.
         *
         * @return Relative MIP gap (0 means proven optimal)
         */
        double get_mip_gap () override
        {
            auto info = model_.getInfo();
            return info.mip_gap;
        }
        
        /**
         * @brief Exports the optimization model to a file
         *
         * Writes the current LP/MIP model to disk in a standard format (e.g., MPS, LP)
         * for debugging, inspection, or solving with external tools.
         *
         * @param filename Path to the output file (format determined by extension)
         */
        void write_model (std::string filename) override
        {
            model_.writeModel (filename);
        }

        /**
         * @brief Returns the number of variables in the optimization model
         *
         * @return Total count of decision variables (columns) in the model
         */
        size_t num_vars() override
        {
            return model_.getNumCol();
        }
        /**
         * @brief Returns the number of constraints in the optimization model
         *
         * @return Total count of constraints (rows) in the model
         */
        size_t num_constraints() override
        {
            return model_.getNumRow();
        }

        /**
         * @brief Controls solver output verbosity
         *
         * Disables or enables solver output messages and progress logs.
         *
         * @param quiet If true, suppresses solver output; if false, enables output
         */
        void local_set_quiet (bool quiet) override
        {
            check_highs_status (
                model_.setOptionValue("output_flag", false),
                "setOptionValue-output_flag"
            );
            
        }
        
    protected:
        Highs model_; /**< HiGHS solver instance managing the optimization model */

        // The flux and use far indices for each reaction in each experiment
        std::vector<std::vector<HighsInt>> flux_var_id_; /**< Variable IDs for reaction fluxes indexed by [reaction][experiment] */
        std::vector<HighsInt> use_var_id_; /**< Binary usage indicator variable IDs indexed by reaction */

        /* Working area for calls to "addRow" */
        double *highs_coeffs_; /**< Working buffer for constraint coefficients in addRow calls */
        HighsInt *col_id_; /**< Working buffer for column indices in addRow calls */

        SolutionPtr cts_sol_; /**< Pointer to the continuous relaxation solution for MIP warm-starting */

        std::vector<HighsInt> biomass_id_; /**< Variable IDs for biomass objective flux in each experiment */
    };
    /**
     * @typedef HighsLPSolverImpPtr
     * @brief Shared pointer type for HighsLPSolverImp instances
     */
    using HighsLPSolverImpPtr = std::shared_ptr<HighsLPSolverImp>;
}
