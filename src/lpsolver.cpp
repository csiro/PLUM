#include <sstream>
#include <list>
#include <iomanip>

#include "lime/debug.h"
#include "lime/numutil.h"
#include "lime/fileutil.h"
#include "lime/error.h"

#include "mosh/constants.h"
#include "mosh/lpsolver.h"
#include "mosh/reaction.h"

using namespace std;
using namespace lime;
using namespace mosh;

int LPSolver::num_model_writes_ = 0;

void
LPSolver::formulate (const Experiment* exp)
{
    DEBUG ('L', "LP formulation for experiment " << exp->name());
    on_screen ("LP Formulation");

    imp_->init_cts();
    
    biomass_obj_mult_ = params_->init_biomass_obj_mult;
    
    // Add vars for flux for each selected reaction
    for (size_t k = 0; k < num_reactions(); k++) {
        auto react = scenario_->reaction(k);
        DEBUG (
            'L', "  React " << k << " " << react->name() <<
            " obj " << react->obj_coeff() <<
            " selected " << react->is_selected()
        );
        // Add vars for flux for every reaction
        double lb = 0.0f;
        double ub = 0.0f;
        // Only allow flux if it is selected
        if (reaction(k)->is_selected()) {
            ub = react->flux_ub();
            
            if (react->has_known_flux() && !react->is_biomass()) {
                double error = 0.0f;
                if (params_->use_error)
                    error = react->known_flux_error();
                else
                    error =  params_->epsilon;
                lb = react->known_flux() - error;
                ub = react->known_flux() + error;
                if (lb < 0.0f)
                    lb = 0.0f;
            }
        }
        double obj_coeff = react->obj_coeff();
        if (react->is_biomass()) {
            switch (which_formulation_) {
            case 1: 
                obj_coeff *= biomass_obj_mult_;
                break;
            case 2:
                lb = target_flux_;
                break;
            case 3:
                obj_coeff = 0.0f;
                lb = target_flux_;
                ub = target_flux_;
                break;
            }
        }
        DEBUG (
            'L', "      Adding var for reaction " << k <<
            " " << react->name() <<
            " obj-coeff " << obj_coeff << 
            " lb " << lb << " ub " << ub
        );
        
        imp_->make_flux_var (
            react, exp->index(), react->name(), obj_coeff, lb, ub
        );
    }
    
    // Add constraints for each metabolite
    DEBUG ('L', "  Add metabolite constraints");
    auto biomass = scenario_->biomass_react();
    
    for (size_t m = 0; m < num_metabolites(); m++) {
        auto met = metabolite(m);
        DEBUG (
            'L', "    Add constraints for " << met->name() <<
            " exp " << exp->name()
        );
        
        vector<Reaction*> reacts;
        vector<double> coeffs;
        
        for (size_t k = 0; k < num_reactions(); k++) {
            auto react = reaction(k);
            if (react->met_has_coeff (m)) {
                coeffs.push_back (react->met_coeff (m));
                reacts.push_back (react);
                DEBUG (
                    'L', "      React " << *react <<
                    " coeff " << react->met_coeff (m)
                );
            }
        }
        double lb = exp->lb(met);
        double ub = exp->ub(met);
        
        if (biomass->uses(met)) {
            // Add wriggle-room
            double coeff = fabs(biomass->met_coeff (met));
            double wiggle = coeff * params_->biomass_wiggle_pc / 100.0f;
            ub += wiggle;
            //lb -= wiggle;
            DEBUG ('L', "        Allowing wiggle " << wiggle);
        }
        DEBUG (
            'L', "      " << met->name() << " lb " << lb << " ub " << ub 
        );
        imp_->add_met_constraint (
            met, exp->index(), lb, ub, reacts, coeffs
        );
    }
    
    DEBUG ('L', "Cts formulation completed");
    imp_->finalise_formulation();
}

void
LPSolver::set_which_formulation (size_t which_formulation)
{
    if (which_formulation_ == which_formulation)
        return; // Nothing to do

    which_formulation_ = which_formulation;
    if (which_formulation == 1) {
        limeCrash ("Can't (yet) switch back to formulation 1");
    }
    else if (which_formulation == 2) {
        // Switch to formulation 2
        set_biomass_obj_mult (which_experiment_, 0.0f);
        size_t idx = scenario_->biomass_react()->index();
        imp_->set_react_bounds (
            scenario_->biomass_react(), which_experiment_,
            target_flux_, scenario_->biomass_react()->flux_ub()
        );
    }
    else {
        // Switch to formulation 2
        set_biomass_obj_mult (which_experiment_, 0.0f);
        size_t idx = scenario_->biomass_react()->index();
        imp_->set_react_bounds (
            scenario_->biomass_react(), which_experiment_,
            target_flux_, target_flux_
        );
    }
}

SolutionPtr
LPSolver::solve ()
{
    formulate (experiment(which_experiment_));
        
    return do_solve();
}

void
LPSolver::enable_reaction (const Reaction* react)
{
    imp_->set_react_bounds (react, which_experiment_, 0.0f, react->flux_ub());
}

void
LPSolver::disable_reaction (const Reaction* react)
{
    imp_->set_react_bounds (react, which_experiment_, 0.0f, 0.0f);
}

SolutionPtr
LPSolver::do_solve()
{
    DEBUG (
        'l', "    LP solve exp " << which_experiment_ << " " <<
        scenario_->experiment(which_experiment_)->name()
    );
    on_screen (
        "  Solving LP with " +
        to_string (imp_->num_vars()) + " vars and " +
        to_string (imp_->num_constraints()) + " constraints" +
        " biomass obj mult " + to_string (biomass_obj_mult_)
    );
    
    if (Debug::doDebug('w')) {
        num_model_writes_++;
        string fn = "lpsolver-" + to_string(num_model_writes_) + ".lp";
        DEBUG ('w', "Writing " + fn);
        imp_->write_model (fn);
    }      

    DEBUG ('l', "      Optimize");
    status_ = imp_->optimize();
    
    if (status_ == INFEASIBLE_) {
        DEBUG ('A', "The model cannot be solved");
        return nullptr;
    }

    double obj = imp_->get_objective();
    DEBUG ('l', "        Solved. Obj value " << obj);
    on_screen ("    Finished. Objective is " + to_string (obj));

    auto sol = imp_->make_sol(which_experiment_);
 
    if (params_->max_biomass_search_iters > 0) {
        sol = biomass_obj_search (sol);
    }

    return sol;
}

SolutionPtr 
LPSolver::biomass_obj_search (SolutionPtr sol)
{
    DEBUG (
        'l', "Biomass obj search - num iters: " <<
        params_->max_biomass_search_iters
    );
    biomass_search_iters_ = 0;

    double curr_obj_mult = params_->init_biomass_obj_mult;

    double lb = curr_obj_mult;
    double ub = -1.0f; //  No UB found yet
    double sol_mult = curr_obj_mult;
    double prev_biomass_flux = 0.0f;
    double best_dummy_flux = LIME_BIG_DOUBLE;

    imp_->local_set_quiet (true);
    
    DEBUG ('g', "Positive biomass search");
    on_screen ("Check biomass_flux");
    on_screen ("  biomass flux is " + to_string (sol->biomass_flux()));
    while (biomass_search_iters_ < params_->max_biomass_search_iters) {
        DEBUG (
            'g', "Biomass search iter " << biomass_search_iters_ <<
            " lb " << lb << " ub " << ub << " mult " << curr_obj_mult <<
            " biomass " << sol->biomass_flux() << " dummy " << sol->dummy_flux()
        );
        if (limeIsZero (sol->biomass_flux())) {
            lb = curr_obj_mult;
            prev_biomass_flux = 0.0f;
            if (ub < 0.0f) // No UB seen
                curr_obj_mult *= params_->biomass_search_mult1;
            else 
                curr_obj_mult = lb + (ub - lb) / 2.0;
            DEBUG (
                'g', "    Biomax flux is still zero - increase mult to " <<
                curr_obj_mult
            );
        }
        else if (params_->use_dummy && !limeIsZero(sol->dummy_flux())) {
            // We've gone too far
            if (sol->dummy_flux() < best_dummy_flux) {
                best_dummy_flux = sol->dummy_flux();
                sol_mult = curr_obj_mult;
            }
            ub = curr_obj_mult;
            curr_obj_mult = lb + (ub - lb) / 2.0;
            prev_biomass_flux = 0.0f;
            DEBUG (
                'g', "    Dummy flux is " << sol->dummy_flux() <<
                " - decrease mult to " << curr_obj_mult <<
                " best dummy flux is " << best_dummy_flux
            );
        }
        else {
            // Positive biomass flux. See if we have converged
            DEBUG ('g', "    Biomax flux is " << sol->biomass_flux());
            if (limeDblEqual (prev_biomass_flux, sol->biomass_flux())) {
                DEBUG ('g', "  Flux converged - stop");
                on_screen ("  Flux converged - stop");
                break;
            }
            // Otherwise, try a little higher
            sol_mult = curr_obj_mult;
            prev_biomass_flux = sol->biomass_flux();
            best_dummy_flux = 0.0f;
            lb = curr_obj_mult;
            if (ub < 0.0f)  // No UB seen
                curr_obj_mult *= params_->biomass_search_mult2;
            else 
                curr_obj_mult = lb + (ub - lb) / 2.0;
            DEBUG (
                'g', "    Biomax flux not converged - increase mult to " <<
                curr_obj_mult
            );
        }
        biomass_search_iters_++;
        
        set_biomass_obj_mult (which_experiment_, curr_obj_mult);
        
        on_screen (
            limeFormat (
                "New biomass mult %g lb %d ub %d", 
                curr_obj_mult, (int)lb, (int) ub
            )
        );

        StatusEnum status = imp_->optimize();
        if (status == INFEASIBLE_)
            limeCrash ("Got infeasible model??");
        
        if (Debug::doDebug('w')) {
            string model_fn =
                unique_filename ("lpsolver-biomass-%03d.lp", 1000);
            if (model_fn.length() > 0) {
                DEBUG ('w', "  Write model in " << model_fn);
                imp_->write_model (model_fn);
            };
        }
        sol = imp_->make_sol(which_experiment_);
        const char* is_runaway =
            params_->is_runaway_flux(sol->biomass_flux())
            ? "yes"
            : "no";
        on_screen (
            "  iter " + to_string (biomass_search_iters_) +
            " biomass flux " + to_string (sol->biomass_flux()) +
            " dummy flux " + to_string (sol->dummy_flux()) +
            " is runaway " + is_runaway
        );
        DEBUG (
            'g', "  iter " + to_string (biomass_search_iters_) +
            " biomass flux " + to_string (sol->biomass_flux()) +
            " dummy flux " + to_string (sol->dummy_flux()) +
            " is runaway " + is_runaway
        );
    }
    if (
        !limeIsZero(sol->dummy_flux()) ||
        (
            limeIsZero(sol->biomass_flux()) &&
            best_dummy_flux < LIME_BIG_DOUBLE
        )
    ) {
        // Return to last-known non-dummy solution
        biomass_search_iters_++;
        curr_obj_mult = sol_mult;
        set_biomass_obj_mult (which_experiment_, curr_obj_mult);
        on_screen (
            "Return to feasibility - biomass mult " +
            to_string (curr_obj_mult)
        );
        StatusEnum status = imp_->optimize();
        if (status == INFEASIBLE_)
            limeCrash ("Got infeasible model??");
        
        sol = imp_->make_sol(which_experiment_);
        for (auto& react : scenario_->reactions()) {
            if (react->is_dummy() && sol->flux(react->index()) > 0.0f) {
                on_screen (
                    "Still have dummy flux of " + to_string(sol->flux(react->index())) + 
                    " in " + react->name()
                );
            }
        }
    }
    biomass_obj_mult_ = curr_obj_mult;
    on_screen (
        "  Biomass obj mult " + to_string (biomass_obj_mult_)
    );
    DEBUG ('l', "Biomass obj mult is " << biomass_obj_mult_);
    DEBUG ('l', "Biomass flux is " << sol->biomass_flux());
    if (!quiet())
        imp_->local_set_quiet (false);
    return sol;
}

string
LPSolver::summary() 
{
    stringstream str;
    str << " biomass_obj_mult " <<
        fixed << setprecision(1) << biomass_obj_mult_ <<
        " biomass_search_iters " << biomass_search_iters_;
    return str.str();
}

