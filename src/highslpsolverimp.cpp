
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

void
HighsLPSolverImp::set_biomass_mult (size_t exp, double mult) 
{
    assert(biomass_id_[exp] >= 0);
    check_highs_status (
        model_.changeColCost(biomass_id_[exp], mult),
        "changeBiomassCost"
    );
}

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

