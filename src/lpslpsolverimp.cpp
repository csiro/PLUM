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
