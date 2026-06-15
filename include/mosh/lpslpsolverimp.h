/**
 * @file lpslpsolverimp.h
 * @brief LPSolve-based implementation of the LP solver interface for metabolic network optimization
 *
 * This file defines the LpsLPSolverImp class, which provides a concrete implementation
 * of the LPSolverImp interface using the lp_solve library for solving linear and
 * mixed-integer linear programs in metabolic gap-filling and flux balance analysis.
 */

#pragma once


#include "lp_lib.h"

#include "mosh/lpsolverimp.h"

#include "lime/numutil.h"

/**
 * @namespace mosh
 * @brief Main namespace for metabolic optimization and solver handling
 */
namespace mosh
{
    /**
     * @class LpsLPSolverImp
     * @brief LPSolve-based implementation of the LP solver for metabolic network problems
     *
     * This class wraps the lp_solve library to provide linear programming and mixed-integer
     * programming capabilities for flux balance analysis and metabolic gap-filling. It handles
     * the creation of flux variables, metabolite balance constraints, and reaction bounds
     * required for metabolic network optimization.
     */
    class LpsLPSolverImp : public LPSolverImp
    {
    public:
        /**
         * @brief Constructor for LPSolve solver implementation
         * @param scenario Pointer to the metabolic scenario containing reactions and metabolites
         * @param params Pointer to solver parameters including time limits and options
         */
        LpsLPSolverImp (Scenario* scenario, Params* params) :
            LPSolverImp (scenario, params),
            model_(NULL),
            add_row_mode_set_(false),
            expect_num_rows_(0),
            lps_coeffs_(new REAL[scenario_->num_reactions()]),
            col_id_(new int[scenario_->num_reactions()]),
            biomass_id_(-1)
        {
        }
        
        /**
         * @brief Destructor that cleans up lp_solve model and allocated arrays
         */
        virtual ~LpsLPSolverImp()
        {
            DEBUG ('A', "Cleanup");
            if (model_ != NULL)
                delete_lp (model_);
            delete [] lps_coeffs_;
            delete [] col_id_;
        }

        /**
         * @brief Initialize continuous linear program for flux balance analysis
         *
         * Sets up the lp_solve model for a continuous optimization problem with
         * variables representing reaction fluxes. Creates constraints for metabolite
         * balance equations.
         */
        void init_cts() override
        {
            // Set up for continuous prob
            model_ = make_lp(0, scenario_->num_reactions());
            expect_num_rows_ = 2 * scenario_->num_metabolites();
            const char* name = "plum";
            set_lp_name (model_, (char*)name);
            set_verbose (model_, IMPORTANT);
            if (params_->time_limit > 0) {
                long long_limit = (long) params_->time_limit;
                set_timeout (model_, long_limit);
            }
        }

        /**
         * @brief Initialize mixed-integer linear program for gap-filling
         * @param cts_sol Solution from the continuous relaxation (may be used for initialization)
         *
         * Sets up the lp_solve model for a mixed-integer problem with both flux variables
         * and binary indicator variables for reaction usage. Includes linking constraints
         * between flux and indicator variables.
         */
        void init_int(SolutionPtr cts_sol) override
        {
            if (model_ != NULL)
                delete_lp (model_);
            model_ = make_lp(0, 2 * scenario_->num_reactions());
            const char* name = "plum-int";
            set_lp_name (model_, (char*)name);
            set_verbose (model_, IMPORTANT);
            if (params_->time_limit > 0) {
                long long_limit = (long) params_->time_limit;
                set_timeout (model_, long_limit);
            }
            // The number of constraints is
            // - 1 per reaction (linking)
            // - 2 per metabolite per experiment (balance: LB + UB)
            expect_num_rows_ =
                scenario_->num_reactions() +
                2 * (scenario_->num_metabolites() + scenario_->num_experiments());
        }

        /**
         * @brief Create a flux variable for a reaction in a specific experiment
         * @param react Pointer to the reaction
         * @param exp Experiment index
         * @param name Variable name for debugging and output
         * @param obj_coeff Objective function coefficient for this variable
         * @param lb Lower bound on flux value
         * @param ub Upper bound on flux value
         */
        void make_flux_var (
            const Reaction* react, size_t exp, std::string name,
            double obj_coeff, double lb, double ub
        ) override;
        /**
         * @brief Create a binary indicator variable for reaction usage
         * @param react Pointer to the reaction
         * @param name Variable name for debugging and output
         * @param obj_coeff Objective function coefficient (typically a penalty cost)
         * @param lb Lower bound (typically 0 for binary variables)
         * @param ub Upper bound (typically 1 for binary variables)
         */
        void make_use_var (
            const Reaction* react, std::string name, double obj_coeff,
            double lb, double ub
        ) override;
        /**
         * @brief Add metabolite mass balance constraint for an experiment
         * @param met Pointer to the metabolite
         * @param experiment Experiment index
         * @param lb Lower bound on metabolite production/consumption rate
         * @param ub Upper bound on metabolite production/consumption rate
         * @param reacts Vector of reactions involving this metabolite
         * @param coeffs Stoichiometric coefficients for each reaction
         */
        void add_met_constraint (
            const Metabolite* met, size_t experiment, double lb, double ub,
            std::vector<Reaction*>& reacts, std::vector<double>& coeffs
        ) override;
        /**
         * @brief Add linking constraint between flux and indicator variables
         * @param react Pointer to the reaction
         * @param exp Experiment index
         *
         * Creates a constraint that enforces flux = 0 when indicator = 0,
         * allowing non-zero flux only when the reaction is marked as used.
         */
        void add_react_link_constraint (
            const Reaction* react, size_t exp
        ) override;
        /**
         * @brief Set or update flux bounds for a reaction in an experiment
         * @param react Pointer to the reaction
         * @param exp Experiment index
         * @param lb New lower bound on flux
         * @param ub New upper bound on flux
         */
        void set_react_bounds (
            const Reaction* react, size_t exp, double lb, double ub
        ) override;
        /**
         * @brief Set or update objective coefficient for a reaction flux
         * @param exp Experiment index
         * @param react Pointer to the reaction
         * @param cost Objective coefficient (e.g., gap-filling penalty or flux weight)
         */
        void set_react_cost (
            size_t exp, const Reaction* react, double cost
        ) override;
        /**
         * @brief Set multiplier for biomass objective function
         * @param exp Experiment index
         * @param mult Multiplier value (typically 1 for maximization, -1 for constraints)
         */
        void set_biomass_mult (size_t exp, double mult) override;
        
        /**
         * @brief Solve the linear or mixed-integer program
         * @return Status code indicating optimal, infeasible, time limit, or other outcomes
         */
        StatusEnum optimize() override;

        /**
         * @brief Extract solution flux values for a specific experiment
         * @param which_experiment Experiment index to extract solution for
         * @return Shared pointer to Solution object containing flux values
         */
        SolutionPtr make_sol(size_t which_experiment) override;
        /**
         * @brief Get the objective function value from the solution
         * @return Objective value (e.g., biomass production or gap-filling cost)
         */
        double get_objective () override
        {
            return ::get_objective(model_);
        }
        /**
         * @brief Get the MIP optimality gap
         * @return Relative gap between best solution and best bound, or -1 if not available
         *
         * Note: lp_solve does not provide MIP gap information, so this always returns -1.
         */
        double get_mip_gap () override
        {
            return -1.0f;
        }
        
        /**
         * @brief Write the LP model to a file for debugging or external analysis
         * @param filename Path to output file (typically .lp format)
         */
        void write_model (std::string filename) override
        {
            write_lp (model_, (char*)filename.c_str());
        }

        /**
         * @brief Get the number of variables in the model
         * @return Total number of decision variables (fluxes and indicators)
         */
        size_t num_vars() override
        {
            return get_Ncolumns(model_);
        }
        /**
         * @brief Get the number of constraints in the model
         * @return Total number of constraints (metabolite balances and linking constraints)
         */
        size_t num_constraints() override
        {
            return get_Nrows(model_);
        }

    protected:
        lprec* model_; /**< Pointer to lp_solve model structure */

        bool add_row_mode_set_; /**< Flag indicating whether batch row-add mode is active */
        int expect_num_rows_; /**< Expected number of constraints for memory pre-allocation */
        REAL *lps_coeffs_; /**< Temporary array for constraint coefficients */
        int *col_id_; /**< Temporary array for column indices in sparse constraints */

        int biomass_id_; /**< Column index of the biomass reaction flux variable */"
    };
    /**
     * @typedef LpsLPSolverImpPtr
     * @brief Shared pointer type for LpsLPSolverImp instances
     */
    using LpsLPSolverImpPtr = std::shared_ptr<LpsLPSolverImp>;
}
