/**
 * @file lpslpsolverimp.cpp
 * @brief Implementation of LP solver wrapper using lp_solve library for metabolic network optimization
 *
 * This file contains the implementation of LpsLPSolverImp class methods, which provides
 * an interface to the lp_solve library for solving flux balance analysis (FBA) and
 * metabolic gap-filling problems. It handles both continuous and mixed-integer programming
 * formulations for optimizing metabolic network flux distributions.
 */
#ifdef PLUM_LPSOLVER

#include <sstream>
#include <list>
#include <iomanip>

#include "lime/debug.h"
#include "lime/numutil.h"
#include "lime/error.h"

#include "mosh/constants.h"
#include "mosh/lpsutil.h"
#include "mosh/lpslpsolverimp.h"
#include "mosh/reaction.h"

using namespace std;
using namespace lime;
using namespace mosh;

/**
 * @brief Creates a binary decision variable for reaction usage in mixed-integer formulations
 *
 * Creates a binary variable that indicates whether a reaction is active (used) in the
 * optimal solution. This is used in gap-filling problems to minimize the number of
 * reactions added to make the network feasible.
 *
 * @param react Pointer to the reaction for which the use variable is created
 * @param name Name identifier for the variable in the LP model
 * @param obj_coeff Objective function coefficient (typically penalty cost for using this reaction)
 * @param lb Lower bound for the variable (typically 0 for binary)
 * @param ub Upper bound for the variable (typically 1 for binary)
 */
void
LpsLPSolverImp::make_use_var (
    const Reaction* react, string name, double obj_coeff, double lb, double ub
)
{
    int var_id = react->lps_id() + scenario_->num_reactions();
        
    REAL lps_obj_coeff = (REAL)obj_coeff;
    REAL lps_lb = (REAL)lb;
    REAL lps_ub = (REAL)ub;
    
    check_lps_return (
        set_column (model_, var_id, &lps_obj_coeff),
        "set_column"
    );
    check_lps_return (
        set_col_name (model_, var_id, (char*)name.c_str()),
        "set_col_name"
    );
    check_lps_return (
        set_bounds (model_, var_id, lps_lb, lps_ub),
        "set_bounds"
    );
    check_lps_return (
        set_binary (model_, var_id, TRUE),
        "set_binary"
    );
}

/**
 * @brief Creates a continuous flux variable for a reaction in the metabolic network
 *
 * Creates a continuous variable representing the flux (reaction rate) through a given
 * reaction. The flux variable is bounded by the reaction's minimum and maximum flux
 * constraints and may have an objective coefficient for optimization.
 *
 * @param react Pointer to the reaction for which the flux variable is created
 * @param exp Experiment index for multi-experiment scenarios
 * @param name Name identifier for the variable in the LP model
 * @param obj_coeff Objective function coefficient (e.g., for maximizing biomass production)
 * @param lb Lower bound for the flux (typically negative for reversible reactions)
 * @param ub Upper bound for the flux
 */
void
LpsLPSolverImp::make_flux_var (
    const Reaction* react, size_t exp, string name,
    double obj_coeff, double lb, double ub
)
{
    int var_id = react->lps_id();
        
    REAL lps_obj_coeff = (REAL)obj_coeff;
    REAL lps_lb = (REAL)lb;
    REAL lps_ub = (REAL)ub;
    
    check_lps_return (
        set_column (model_, var_id, &lps_obj_coeff),
        "set_column"
    );
    check_lps_return (
        set_col_name (model_, var_id, (char*)name.c_str()),
        "set_col_name"
    );
    check_lps_return (
        set_bounds (model_, var_id, lps_lb, lps_ub),
        "set_bounds"
    );
    if (react->is_biomass())
        biomass_id_ = var_id;
}

/**
 * @brief Adds mass balance constraints for a metabolite in the metabolic network
 *
 * Creates two constraints (lower bound and upper bound) that enforce mass balance
 * for a given metabolite. The constraints ensure that the sum of fluxes producing
 * and consuming the metabolite satisfies steady-state assumptions in flux balance analysis.
 *
 * @param met Pointer to the metabolite for which constraints are added
 * @param exp Experiment index for multi-experiment scenarios
 * @param lb Lower bound for the metabolite balance (typically 0 for steady-state)
 * @param ub Upper bound for the metabolite balance (typically 0 for steady-state)
 * @param reacts Vector of reactions that participate in this metabolite's balance
 * @param coeffs Vector of stoichiometric coefficients corresponding to each reaction
 */
void
LpsLPSolverImp::add_met_constraint (
    const Metabolite* met, size_t exp, double lb, double ub,
    vector<Reaction*>& reacts, vector<double>& coeffs
)
{
    if (!add_row_mode_set_) {
        // Start of constraints
        set_add_rowmode(model_, TRUE);
        add_row_mode_set_ = true;
        check_lps_return (
            resize_lp (
                model_,  expect_num_rows_,
                get_Ncolumns(model_)
            ),
            "resize_lp"
        );
    }
    REAL lps_lb = (REAL)lb;
    REAL lps_ub = (REAL)ub;

    int num_coeffs = (int)coeffs.size();
    for (size_t k = 0; k < coeffs.size(); k++) {
        lps_coeffs_[k] = (REAL)coeffs[k];
        col_id_[k] = reacts[k]->lps_id();
    }

    int constraint_id = get_Nrows(model_);
    check_lps_return (
        add_constraintex (model_, num_coeffs, lps_coeffs_, col_id_, GE, lps_lb),
        "add_constraintex"
    );
    string name = met->name() + "_LB";
    check_lps_return (
        set_row_name (model_, constraint_id, (char*)name.c_str()),
        "set_row_name"
    );
    constraint_id++;
    
    check_lps_return (
        add_constraintex (model_, num_coeffs, lps_coeffs_, col_id_, LE, lps_ub),
        "add_constraintex"
    );
    name = met->name() + "_UB";
    check_lps_return (
        set_row_name (model_, constraint_id, (char*)name.c_str()),
        "set_row_name"
    );
}

/**
 * @brief Adds linking constraint between flux and binary use variables for a reaction
 *
 * Creates a constraint that links the continuous flux variable with the binary use variable
 * for a reaction in mixed-integer formulations. The constraint enforces that if the flux
 * is non-zero, the use variable must be 1, implementing the big-M formulation:
 * flux - flux_ub * use <= 0
 *
 * @param react Pointer to the reaction for which the linking constraint is created
 * @param exp Experiment index for multi-experiment scenarios
 */
void
LpsLPSolverImp::add_react_link_constraint (const Reaction* react, size_t exp)
{
    size_t idx = react->index();
    int constraint_id = get_Nrows(model_);

    int flux_id = react->lps_id();
    int use_id = flux_id + scenario_->num_reactions();
        
    // 1.0 * flux 
    col_id_[0] = flux_id;
    lps_coeffs_[0] = 1.0;
    // - flux-ub * use
    col_id_[1] = use_id;
    lps_coeffs_[1] = -react->flux_ub();
    
    // flux - bigM use <= 0
    check_lps_return (
        add_constraintex (model_, 2, lps_coeffs_, col_id_, LE, 0.0),
        "add_constraintex"
    );
    string name = react->name() + "-link";
    check_lps_return (
        set_row_name (model_, constraint_id, (char*)name.c_str()),
        "set_row_name"
    );
}

/**
 * @brief Solves the linear or mixed-integer programming problem
 *
 * Invokes the lp_solve solver to optimize the constructed LP/MIP model. The method
 * handles solver configuration, execution, and status interpretation. It returns
 * the optimization status indicating whether an optimal, suboptimal, or infeasible
 * solution was found.
 *
 * @return StatusEnum indicating the optimization result (CTS_OPTIMAL, SUB_OPTIMAL, or INFEASIBLE_)
 */
StatusEnum
LpsLPSolverImp::optimize()
{
    if (add_row_mode_set_) {
        // End of constraints
        set_add_rowmode(model_, FALSE);
        add_row_mode_set_ = false;
    }

    set_verbose (model_, DETAILED);
    
    int lps_status = ::solve (model_);

    DEBUG ('A', "Finished solve with status " << lps_status);

    set_verbose (model_, NORMAL);
    
    if (
        lps_status != OPTIMAL &&
        lps_status != SUBOPTIMAL &&
        lps_status != TIMEOUT
    ) {
        DEBUG ('A', "The model cannot be solved - status " << lps_status);
        return INFEASIBLE_;
    }

    return (lps_status == OPTIMAL ? CTS_OPTIMAL : SUB_OPTIMAL);
}
    

/**
 * @brief Updates the lower and upper bounds for a reaction flux variable
 *
 * Modifies the bounds on a reaction's flux variable after the variable has been created.
 * This is used for dynamic constraint modification during iterative optimization procedures.
 *
 * @param react Pointer to the reaction whose bounds are being updated
 * @param exp Experiment index for multi-experiment scenarios
 * @param lb New lower bound for the flux variable
 * @param ub New upper bound for the flux variable
 */
void
LpsLPSolverImp::set_react_bounds (
    const Reaction* react, size_t exp, double lb, double ub
)
{
    int var_id = react->lps_id();
    REAL lps_lb = (REAL)lb;
    REAL lps_ub = (REAL)ub;
    set_bounds (model_, var_id, lps_lb, lps_ub);
}

/**
 * @brief Sets the objective function coefficient for a reaction flux variable
 *
 * Updates the objective coefficient for a reaction's flux variable in the LP model.
 * This is used to change optimization objectives dynamically, such as switching from
 * biomass maximization to flux minimization in multi-phase optimization.
 *
 * @param exp Experiment index for multi-experiment scenarios
 * @param react Pointer to the reaction whose objective coefficient is being set
 * @param cost New objective function coefficient for the reaction
 */
void
LpsLPSolverImp::set_react_cost (
    size_t exp, const Reaction* react, double cost
)
{
    int var_id = react->lps_id();
    REAL coeff = cost;
    check_lps_return (
        set_mat (model_, 0, var_id, coeff),
        "set_mat"
    );
}

/**
 * @brief Sets the objective function coefficient for the biomass reaction
 *
 * Updates the objective coefficient specifically for the biomass flux variable.
 * This is commonly used to maximize or constrain biomass production in flux balance
 * analysis, or to set a minimum biomass threshold in gap-filling problems.
 *
 * @param exp Experiment index for multi-experiment scenarios
 * @param mult Multiplier (coefficient) for the biomass objective term
 */
void
LpsLPSolverImp::set_biomass_mult (size_t exp, double mult) 
{
    REAL coeff = mult;
    assert(biomass_id_ >= 0);
    check_lps_return (
        set_mat (model_, 0, biomass_id_, coeff),
        "set_mat"
    );
}

/**
 * @brief Extracts the solution from the optimized LP model
 *
 * Retrieves the optimal flux values from the solved LP/MIP model and constructs
 * a Solution object containing the flux distribution for all reactions in the
 * metabolic network. Only selected reactions have their flux values extracted.
 *
 * @param exp Experiment index for multi-experiment scenarios
 * @return SolutionPtr Shared pointer to a Solution object containing the optimal flux distribution
 */
SolutionPtr
LpsLPSolverImp::make_sol(size_t exp)
{
    SolutionPtr sol = make_shared<Solution> (scenario_, params_);

    REAL* val = new REAL [num_vars()];
    DEBUG ('A', "  Val is " << val);

    check_lps_return (
        get_variables (model_, val),
        "get_variables"
    );
    
    for (size_t k = 0; k < scenario_->num_reactions(); k++) {
        if (scenario_->reaction(k)->is_selected()) {
            sol->set_flux (k, val[k]);
            DEBUG (
                'f', "  Val for flux " << k <<
                " " << scenario_->reaction(k)->name() <<
                " is " << sol->flux(k)
            );
            if (sol->flux(k) > 0.0f) {
                DEBUG (
                    'F', "  Val for flux " << k <<
                    " " << scenario_->reaction(k)->name() <<
                    " is " << sol->flux(k) <<
                    " dummy " << scenario_->reaction(k)->is_dummy()
                 );
            }
        }
    }
    
    DEBUG ('A', "  Val is " << val << " delete");
    delete [] val;
    
    DEBUG ('A', "  Done");
    return sol;
}

#endif // PLUM_LPSOLVER
