/**
 * @file highslpsolverimp.cpp
 * @brief Implementation of HiGHS-based LP solver for metabolic network gap-filling
 *
 * This file contains the implementation of the HighsLPSolverImp class, which provides
 * a concrete implementation of the LPSolverImp interface using the HiGHS optimization library.
 * It handles the formulation and solution of mixed-integer linear programming problems for
 * metabolic network reconstruction, including flux balance analysis (FBA) constraints,
 * reaction usage variables, and multi-experiment scenarios.
 */

#include <sstream>
#include <list>
#include <iomanip>

#include "lime/debug.h"
#include "lime/numutil.h"
#include "lime/error.h"

#include "mosh/constants.h"
#include "mosh/highslpsolverimp.h"
#include "mosh/reaction.h"

using namespace std;
using namespace lime;
using namespace mosh;

/**
 * @brief Create a binary usage variable for a reaction in the optimization model
 *
 * Creates an integer variable that indicates whether a reaction is used (active) in the
 * metabolic network. This variable is shared across all experiments and is used in
 * linking constraints to enforce that flux can only occur if the reaction is selected.
 *
 * @param react Pointer to the reaction for which to create the usage variable
 * @param name Name identifier for the variable in the model
 * @param obj_coeff Objective function coefficient (penalty for using this reaction)
 * @param lb Lower bound for the usage variable (typically 0)
 * @param ub Upper bound for the usage variable (typically 1)
 */
void
HighsLPSolverImp::make_use_var (
    const Reaction* react, string name, double obj_coeff, double lb, double ub
)
{
    HighsInt id = num_vars();
    use_var_id_[react->index()] = id;
    check_highs_status (
        model_.addVar (lb, ub),
        "addVar"
    );
    check_highs_status (
        model_.changeColCost (id, obj_coeff),
        "changeColCost"
    );
    check_highs_status (
        model_.passColName (id, name),
        "passColName"
    );
    check_highs_status (
        model_.changeColIntegrality(id, HighsVarType::kInteger),
        "changeColIntegrality"
    );
}

/**
 * @brief Create a flux variable for a reaction in a specific experiment
 *
 * Creates a continuous variable representing the flux (reaction rate) through a reaction
 * in a particular experimental condition. For biomass reactions, the variable ID is also
 * stored separately for quick access during objective function modifications.
 *
 * @param react Pointer to the reaction for which to create the flux variable
 * @param exp Index of the experimental condition
 * @param name Base name identifier for the variable (experiment index appended if multiple experiments)
 * @param obj_coeff Objective function coefficient for this flux variable
 * @param lb Lower bound for the flux (typically negative for reversible reactions)
 * @param ub Upper bound for the flux (maximum reaction rate)
 */
void
HighsLPSolverImp::make_flux_var (
    const Reaction* react, size_t exp, string name,
    double obj_coeff, double lb, double ub
)
{
    HighsInt id = num_vars();
    flux_var_id_[react->index()][exp] = id;
    if (react->is_biomass()) {
        biomass_id_[exp] = id;
    }
    
    check_highs_status (
        model_.addVar (lb, ub),
        "addVar"
    );
    check_highs_status (
        model_.changeColCost (id, obj_coeff),
        "changeColCost"
    );
    if (num_experiments() > 1)
        name += "_" + to_string(exp);
    check_highs_status (
        model_.passColName (id, name),
        "passColName"
    );
}

/**
 * @brief Add a metabolite mass balance constraint to the model
 *
 * Creates a constraint enforcing steady-state mass balance for a metabolite in a specific
 * experiment. The constraint ensures that the sum of fluxes producing the metabolite equals
 * the sum of fluxes consuming it (weighted by stoichiometric coefficients), implementing
 * the fundamental constraint of flux balance analysis.
 *
 * @param met Pointer to the metabolite for which to create the constraint
 * @param exp Index of the experimental condition
 * @param lb Lower bound for the constraint (typically 0 for strict steady-state)
 * @param ub Upper bound for the constraint (typically 0 for strict steady-state)
 * @param reacts Vector of reactions participating in the metabolite's mass balance
 * @param coeffs Vector of stoichiometric coefficients corresponding to each reaction
 */
void
HighsLPSolverImp::add_met_constraint (
    const Metabolite* met, size_t exp, double lb, double ub,
    vector<Reaction*>& reacts, vector<double>& coeffs
)
{
    HighsInt num_coeffs = (HighsInt)coeffs.size();
    for (size_t k = 0; k < coeffs.size(); k++) {
        highs_coeffs_[k] = coeffs[k];
        col_id_[k] = flux_var_id_[reacts[k]->index()][exp];
    }

    HighsInt constraint_id = (HighsInt)num_constraints();
    check_highs_status (
        model_.addRow (lb, ub, num_coeffs, col_id_, highs_coeffs_),
        "addRow"
    );

    string name = met->name() + "-bounds";
    check_highs_status (
        model_.passRowName (constraint_id, name),
        "passRowName"
    );
}

/**
 * @brief Add a constraint linking reaction flux to its usage indicator variable
 *
 * Creates a big-M constraint that enforces the relationship between a reaction's flux
 * variable and its binary usage variable. The constraint ensures that flux can only be
 * non-zero if the usage variable is 1, implementing the logical constraint:
 * flux <= flux_ub * use, or equivalently: flux - flux_ub * use <= 0
 *
 * @param react Pointer to the reaction for which to create the linking constraint
 * @param exp Index of the experimental condition
 */
void
HighsLPSolverImp::add_react_link_constraint (const Reaction* react, size_t exp)
{
    HighsInt use_id = use_var_id_[react->index()];
    HighsInt flux_id = flux_var_id_[react->index()][exp];
    
    // 1.0 * flux 
    col_id_[0] = flux_id;
    highs_coeffs_[0] = 1.0;
    // - flux-ub * use
    col_id_[1] = use_id;
    highs_coeffs_[1] = -react->flux_ub();
        
    // -inf <= flux - bigM use <= 0
    check_highs_status (
        model_.addRow (-model_.getInfinity(), 0, 2, col_id_, highs_coeffs_),
        "addRow"
    );
}

/**
 * @brief Update the flux bounds for a reaction in a specific experiment
 *
 * Modifies the lower and upper bounds on a reaction's flux variable for a particular
 * experimental condition. This is used to enforce experimental constraints, such as
 * gene knockouts (bounds set to 0) or measured flux ranges.
 *
 * @param react Pointer to the reaction whose bounds are being modified
 * @param exp Index of the experimental condition
 * @param lb New lower bound for the reaction flux
 * @param ub New upper bound for the reaction flux
 */
void
HighsLPSolverImp::set_react_bounds (
    const Reaction* react, size_t exp, double lb, double ub
)
{
    DEBUG (
        'H', "               Set bounds for " << react->name() <<
        " to " << lb << ", " << ub << " in exp " << exp
    );
    HighsInt id = flux_var_id_[react->index()][exp];
    check_highs_status (
        model_.changeColBounds (id, lb, ub),
        "changeColBounds"
    );
}


/**
 * @brief Update the objective function coefficient for a reaction flux variable
 *
 * Modifies the coefficient of a reaction's flux variable in the objective function
 * for a specific experiment. This is used to weight different objectives, such as
 * maximizing biomass production or minimizing total flux.
 *
 * @param exp Index of the experimental condition
 * @param react Pointer to the reaction whose objective coefficient is being modified
 * @param cost New objective function coefficient for the reaction flux
 */
void
HighsLPSolverImp::set_react_cost (
    size_t exp, const Reaction* react, double cost
)
{
    DEBUG (
        'H', "               Set obj coeff for " << react->name() <<
        " to " << cost << " in exp " << exp
    );
    HighsInt id = flux_var_id_[react->index()][exp];
    check_highs_status (
        model_.changeColCost(id, cost),
        "changeColCost"
    );
}

/**
 * @brief Update the objective function coefficient for the biomass reaction
 *
 * Modifies the coefficient of the biomass flux variable in the objective function
 * for a specific experiment. This is commonly used to maximize biomass production
 * during flux balance analysis.
 *
 * @param exp Index of the experimental condition
 * @param mult New objective coefficient multiplier for the biomass flux
 */
void
HighsLPSolverImp::set_biomass_mult (size_t exp, double mult) 
{
    assert(biomass_id_[exp] >= 0);
    check_highs_status (
        model_.changeColCost(biomass_id_[exp], mult),
        "changeBiomassCost"
    );
}

/**
 * @brief Solve the optimization problem using the HiGHS solver
 *
 * Invokes the HiGHS solver to optimize the formulated mixed-integer linear program.
 * The function maps HiGHS status codes to the internal StatusEnum type and handles
 * various solver outcomes including optimal solutions, time limits, and infeasibility.
 *
 * @return StatusEnum indicating the solution status:
 *         - CTS_OPTIMAL: Optimal solution found
 *         - SUB_OPTIMAL: Feasible solution found but not proven optimal (e.g., time limit)
 *         - INFEASIBLE_: Problem is infeasible or solver failed
 */
StatusEnum
HighsLPSolverImp::optimize()
{
    check_highs_status (
        model_.run(),
        "run"
    );
    
    HighsModelStatus status = model_.getModelStatus();
    
    if (
        status != HighsModelStatus::kOptimal &&
        status != HighsModelStatus::kTimeLimit &&
        status != HighsModelStatus::kUnknown
    ) {
        DEBUG (
            'A', "The model cannot be solved - status " << (int)status << " " <<
            model_.modelStatusToString(status)
        );
        cout << "HIGHS model returned status " << (int) status << ": " <<
            model_.modelStatusToString(status) << endl;
        return INFEASIBLE_;
    }

    return  (status == HighsModelStatus::kOptimal ? CTS_OPTIMAL : SUB_OPTIMAL);
}

/**
 * @brief Extract the solution for a specific experiment from the solver
 *
 * Creates a Solution object containing the flux values for all selected reactions
 * in a particular experimental condition. If the solver did not find a feasible solution,
 * an empty Solution object is returned. Small flux values below the integrality tolerance
 * threshold are rounded to zero.
 *
 * @param exp Index of the experimental condition for which to extract the solution
 * @return SolutionPtr Shared pointer to a Solution object containing flux values
 */
SolutionPtr
HighsLPSolverImp::make_sol(size_t exp)
{
    SolutionPtr sol = make_shared<Solution> (scenario_, params_);

    if (
        model_.getModelStatus() != HighsModelStatus::kOptimal &&
        model_.getModelStatus() != HighsModelStatus::kTimeLimit &&
        model_.getModelStatus() != HighsModelStatus::kUnknown
    ) {
        // Solve failed - return empty sol
        return sol;
    }

    const HighsSolution& hsol = model_.getSolution();

    for (size_t k = 0; k < num_reactions(); k++) {
        if (reaction(k)->is_selected()) {
            double this_flux = hsol.col_value[flux_var_id_[k][exp]];
            if (this_flux <= params_->int_feas_tol * 10.0f)
                this_flux = 0.0f;
                
            DEBUG (
                'f', "  Val for flux " << k <<
                " " << reaction(k)->name() <<
                "  experiment  " << exp <<
                " is " << hsol.col_value[k] << 
                " = " << this_flux
            );
            sol->set_flux (k, this_flux);
        }
    }
    return sol;
}

