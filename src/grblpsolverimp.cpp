
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

double
GrbLPSolverImp::reduced_cost (const Reaction* react)
{
    return flux_[react->index()][0].get (GRB_DoubleAttr_RC);
}





