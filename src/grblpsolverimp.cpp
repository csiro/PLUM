/**
 * @file grblpsolverimp.cpp
 * @brief Implementation of Gurobi-based LP solver for metabolic network optimization
 *
 * This file provides the concrete implementation of the GrbLPSolverImp class,
 * which uses the Gurobi optimizer to solve linear and mixed-integer programming
 * problems for metabolic gap-filling and flux balance analysis. It handles
 * variable creation for reaction fluxes, metabolite mass-balance constraints,
 * and optimization of metabolic network models.
 */

#include <sstream>
#include <list>
#include <iomanip>

#include "lime/debug.h"
#include "lime/numutil.h"
#include "lime/error.h"

#include "mosh/constants.h"
#include "mosh/grblpsolverimp.h"
#include "mosh/reaction.h"

using namespace std;
using namespace lime;
using namespace mosh;


/**
 * @brief Create a continuous decision variable representing flux through a reaction
 *
 * Creates a Gurobi continuous variable for the flux value of a specific reaction
 * in a given experiment. The variable is stored in the flux_ matrix and indexed
 * by reaction and experiment number. If this is a biomass reaction, the variable
 * is also stored in the biomass_ vector for quick access.
 *
 * @param react Pointer to the reaction for which to create the flux variable
 * @param exp Index of the experimental condition
 * @param name Base name for the variable (will be suffixed with experiment number if multiple experiments)
 * @param obj_coeff Objective function coefficient for this variable
 * @param lb Lower bound on the flux value
 * @param ub Upper bound on the flux value
 */
void
GrbLPSolverImp::make_flux_var (
    const Reaction* react,  size_t exp, string name,
    double obj_coeff, double lb, double ub
)
{
    if (num_experiments() > 1)
        name += "_" + to_string(exp);
    flux_[react->index()][exp] =
        model_.addVar (lb, ub, obj_coeff, GRB_CONTINUOUS, name);
    if (react->is_biomass())
        biomass_[exp] = flux_[react->index()][exp];
}

/**
 * @brief Create a binary decision variable indicating whether a reaction is used
 *
 * Creates a Gurobi binary variable that indicates whether a reaction is active
 * in the solution (used in mixed-integer formulations for gap-filling). If a
 * continuous solution is provided, the variable's starting value is initialized
 * based on whether the reaction was used in that solution to warm-start the MIP.
 *
 * @param react Pointer to the reaction for which to create the use variable
 * @param name Name for the variable
 * @param obj_coeff Objective function coefficient (typically penalty for adding reactions)
 * @param lb Lower bound (typically 0 for binary variables)
 * @param ub Upper bound (typically 1 for binary variables)
 */
void
GrbLPSolverImp::make_use_var (
    const Reaction* react,  string name,
    double obj_coeff, double lb, double ub
)
{
    use_[react->index()] =
        model_.addVar (lb, ub, obj_coeff, GRB_BINARY, name);
    if (cts_sol_ != nullptr) {
        double start_guess = 0.0f;
        if (cts_sol_->uses_react(react))
            start_guess = 1.0f;
        use_[react->index()].set (GRB_DoubleAttr_Start, start_guess);
    }
}

/**
 * @brief Add mass-balance constraints for a metabolite in a specific experiment
 *
 * Creates linear constraints enforcing the steady-state mass-balance condition
 * for a metabolite: the sum of fluxes producing the metabolite must equal the
 * sum of fluxes consuming it. This is implemented as two inequality constraints
 * (lower and upper bounds) on the net production rate.
 *
 * @param met Pointer to the metabolite for which to add constraints
 * @param exp Index of the experimental condition
 * @param lb Lower bound on the net production rate (typically 0 for steady-state)
 * @param ub Upper bound on the net production rate (typically 0 for steady-state)
 * @param reacts Vector of reactions that involve this metabolite
 * @param coeffs Vector of stoichiometric coefficients (positive for production, negative for consumption)
 */
void
GrbLPSolverImp::add_met_constraint (
    const Metabolite* met, size_t exp, double lb, double ub,
    vector<Reaction*>& reacts, vector<double>& coeffs
)
{
    GRBLinExpr sum_expr = 0;

    for (size_t k = 0; k < reacts.size(); k++) {
        sum_expr += coeffs[k] * flux_[reacts[k]->index()][exp];
    }

    string name = met->name();
    if (num_experiments() > 1)
        name += "_" + to_string(exp);

    constr_.push_back (
        model_.addConstr (lb <= sum_expr, name + "_LB")
    );
    constr_.push_back (
        model_.addConstr (sum_expr <= ub, name + "_UB")
    );
}

/**
 * @brief Add indicator constraint linking reaction usage and flux variables
 *
 * Creates a logical constraint that forces the flux through a reaction to be
 * zero when the reaction's use variable is 0. This is implemented as an indicator
 * constraint: if use[react] = 0, then flux[react][exp] = 0. This constraint type
 * is used in mixed-integer formulations for gap-filling.
 *
 * @param react Pointer to the reaction for which to add the linking constraint
 * @param exp Index of the experimental condition
 */
void
GrbLPSolverImp::add_react_link_constraint (const Reaction* react, size_t exp)
{
    string name = react->name();
    if (num_experiments() > 1)
        name += "_" + to_string(exp);

    model_.addGenConstrIndicator (
        use_[react->index()],
        0, flux_[react->index()][exp], GRB_EQUAL, 0.0f,
        name + "-link"
    );
}

/**
 * @brief Solve the optimization problem using Gurobi
 *
 * Updates the model and invokes the Gurobi optimizer to solve the LP or MIP.
 * If the problem is infeasible, computes an Irreducible Inconsistent Subsystem
 * (IIS) and writes it to a file for debugging. Returns a status code indicating
 * whether an optimal or sub-optimal solution was found, or if the problem is
 * infeasible.
 *
 * @return StatusEnum indicating the solution status (CTS_OPTIMAL, SUB_OPTIMAL, or INFEASIBLE_)
 */
StatusEnum
GrbLPSolverImp::optimize()
{
    model_.update();

    DEBUG (
        'g', "  Solving problem with " << num_vars() <<
        " vars and " << num_constraints() << " constraints"
    )
    model_.optimize();

    int grb_status = model_.get(GRB_IntAttr_Status);
    if (grb_status == GRB_INFEASIBLE) {
        model_.computeIIS();
        model_.write("infease.ilp");
    }
    if (grb_status != GRB_OPTIMAL && grb_status != GRB_TIME_LIMIT) {
        DEBUG ('A', "The model cannot be solved - status " << grb_status);
        return INFEASIBLE_;
    }
    return (grb_status == GRB_OPTIMAL ? CTS_OPTIMAL : SUB_OPTIMAL);
}
    

/**
 * @brief Extract solution flux values for a specific experiment from Gurobi model
 *
 * Creates a Solution object and populates it with the optimal flux values
 * obtained from the Gurobi solver for a given experiment. Flux values below
 * a small threshold (10 * int_feas_tol) are treated as zero to handle
 * numerical precision issues.
 *
 * @param exp Index of the experimental condition to extract solution for
 * @return SolutionPtr Shared pointer to the Solution object containing flux values
 */
SolutionPtr
GrbLPSolverImp::make_sol(size_t exp)
{
    SolutionPtr sol =
        make_shared<Solution> (scenario_, params_);

    for (size_t k = 0; k < num_reactions(); k++) {
        if (reaction(k)->is_selected()) {
            double this_flux = flux_[k][exp].get(GRB_DoubleAttr_X);
            if (this_flux <= params_->int_feas_tol * 10.0f)
                this_flux = 0.0f;
            DEBUG (
                'f', "  Val for flux " << k <<
                " " << reaction(k)->name() <<
                "  experiment  " << exp <<
                " is " << flux_[k][exp].get(GRB_DoubleAttr_X) <<
                " = " << this_flux
            );
            sol->set_flux (k, this_flux);
        }
    }
    return sol;
}


// Return dual vals for each constraint
// Constraint 2m   = LB for metabolite m
// Constraint 2m+1 = UB for metabolite m
/**
 * @brief Retrieve dual values (shadow prices) for metabolite mass-balance constraints
 *
 * Extracts the dual values from the Gurobi LP solution for all metabolite
 * constraints. Each metabolite has two constraints (lower and upper bound on
 * net production), so dual values are retrieved for both. These dual values
 * represent the shadow prices and can be used for economic analysis or to
 * identify limiting metabolites.
 *
 * @return DualValsPtr Shared pointer to DualVals object containing dual values for all metabolites
 */
DualValsPtr
GrbLPSolverImp::get_dual_vals()
{
    DualValsPtr vals = make_shared<DualVals>(num_metabolites());

    DEBUG ('g', "Get dual vals");
    for (size_t m = 0; m < num_metabolites(); m++) {
        size_t m2 = 2 * m;
        double lb_dual = constr_[m2].get(GRB_DoubleAttr_Pi);
        double ub_dual = constr_[m2+1].get(GRB_DoubleAttr_Pi);
        DEBUG (
            'g', "  Dual vals for " << metabolite(m)->name() <<
            " are LB-dual " << lb_dual << " UB-dual " << ub_dual
        );
        vals->set_duals (m, lb_dual, ub_dual);
    }
    return vals;
}

/**
 * @brief Select reactions with favorable reduced costs for column generation
 *
 * Implements a column generation heuristic by computing reduced costs for all
 * unselected reactions and selecting up to num_to_select reactions with the
 * best (most negative) reduced costs. These reactions are candidates for adding
 * to the model in the next iteration. Reactions with reduced cost at or below
 * the threshold are marked as selected.
 *
 * @param num_to_select Maximum number of reactions to select based on reduced cost
 * @return int Number of reactions newly selected in this call
 */
int
GrbLPSolverImp::select_good_reactions (int num_to_select)
{
    int num_selected = 0;

    DEBUG ('C', "Select good cols");
    vector<double> duals;
    for (size_t k = 0; k < num_reactions(); k++) {
        auto react = reaction(k);
        double cost = reduced_cost (react);
        react->set_reduced_cost(cost);
        if (!react->is_selected()) {
            DEBUG (
                'C', "  Reaction " << *react <<
                " reduced cost " << cost
            );
            duals.push_back (cost);
        }
    }
    // Select a maximum of num_to_select columns
    std::sort(duals.begin(), duals.end());
    double min_val = duals[num_to_select];
    DEBUG ('C', "Min reduced cost for selection is " << min_val);

    for (size_t k = 0; k < num_reactions(); k++) {
        if (!reaction(k)->is_selected()) {
            if (reaction(k)->reduced_cost() <= min_val) {
                // Yay! - column is useful
                DEBUG ('C', "    Select " << *reaction(k));
                reaction(k)->set_selected (true);
                num_selected++;
            }
        }
    }
    DEBUG ('C', "Num newly selected is " << num_selected);
    return num_selected;
}

/**
 * @brief Compute the reduced cost of a reaction in the current LP solution
 *
 * Retrieves the reduced cost from Gurobi for the flux variable of the given
 * reaction. Reduced cost indicates how much the objective function would change
 * per unit increase in the variable. Negative reduced costs indicate reactions
 * that could improve the objective if added to the basis.
 *
 * @param react Pointer to the reaction for which to compute reduced cost
 * @return double The reduced cost value from the LP dual solution
 */
double
GrbLPSolverImp::reduced_cost (const Reaction* react)
{
    return flux_[react->index()][0].get (GRB_DoubleAttr_RC);
}





