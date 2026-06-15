/**
 * @file intsolver.cpp
 * @brief Integer linear programming solver implementation for metabolic gap-filling
 *
 * This file implements the IntSolver class which formulates and solves integer linear
 * programming (ILP) problems for metabolic network gap-filling. It supports multiple
 * formulations for optimizing reaction selection while maintaining flux balance constraints
 * and biomass production requirements.
 */

#include <list>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cmath>

#include "lime/debug.h"
#include "lime/numutil.h"
#include "lime/constants.h"
#include "lime/sorter.h"
#include "lime/rand.h"
#include "lime/error.h"

#include "mosh/constants.h"
#include "mosh/intsolver.h"

using namespace std;
using namespace lime;
using namespace mosh;

/**
 * @brief Formulates the integer linear programming problem for gap-filling
 *
 * Constructs the ILP formulation including:
 * - Flux variables for biomass and metabolic reactions
 * - Binary use variables for reaction selection
 * - Metabolite balance constraints
 * - Reaction linking constraints
 *
 * The specific formulation depends on which_formulation_:
 * - Formulation 1: Optimizes biomass with multiplier
 * - Formulation 2: Fixes lower bound for biomass flux from continuous solution
 * - Formulation 3: Minimizes reaction cost then flux with cost multiplier
 */
void
IntSolver::formulate ()
{
    DEBUG ('L', "Formulate integer prob - formulation #" << which_formulation_);
    imp_->init_int(cts_sol_);

    for (size_t k = 0; k < num_experiments(); k++) {
        experiment(k)->set_base_flux (
            cts_sol_->sol(k)->biomass_flux()
        );
    }
    double rank0_target_flux =
        scenario_->biolog_rank(0)->base_flux() *
        params_->biomass_opt_mult;
    if (!limeIsZero (params_->target_flux))
        rank0_target_flux = params_->target_flux;
    vector<double> target_flux;
    scenario_->calc_target_flux (target_flux, rank0_target_flux);

    // Used for int formulation 3.
    // Formulation 3 miimises sum react cost; THEN sum flux
    // So, multiply react cost by the sum flux in cts sol to make
    // sure sum react cost dominates
    double form_3_mult = 2.0f * cts_sol_->sum_flux();
    if (which_formulation_ == 3)
        DEBUG (
            'L', "  React cost multiplier for formulation 3 is " << form_3_mult
        );


    // Add flux vars for biomass reaction, for each experiment
    if (scenario_->biomass_react() == nullptr)
        limeCrash ("No biomass reaction");

    auto biomass = scenario_->biomass_react();
    for (size_t exp = 0; exp < num_experiments(); exp++) {
        double lb = 0.0f;
        double ub = biomass->flux_ub();
        double obj_coeff = 0;

        switch (which_formulation_) {
        case 1:
            obj_coeff = biomass->obj_coeff() * biomass_obj_mult_;
            break;
        case 2:
        case 3:
            // Fix a lower bound for flux from the cts solutions
            lb = target_flux[exp];
        }
        DEBUG (
            'L', "      Adding cts var for biomass reaction " <<
            biomass->name() <<
            " obj-coeff " << obj_coeff << 
            " lb " << lb << " ub " << ub
        )
        imp_->make_flux_var (
            biomass, exp, biomass->name(), obj_coeff, lb, ub
        );
    }
        
    // Add flux vars for each (non-biomass) reaction
    for (size_t k = 0; k < num_reactions(); k++) {
        auto react = reaction(k);
        if (react->is_biomass()) // Handled above
            continue;
        double lb = 0.0f;
        double ub = 0.0f;
        if (react->is_selected()) {
            DEBUG ('G', "      " << react->name() << " is selected");
            ub = react->flux_ub();
            if (react->has_known_flux()) {
                double error = 0.0f;
                if (params_->use_error)
                    error = react->known_flux_error();
                else
                    error =  params_->epsilon;
                lb = react->known_flux() - error;
                ub = react->known_flux() + error;
                if (lb < 0)
                    lb = 0;
            }
        }
        else {
            DEBUG ('g', "      " << react->name() << " is NOT selected");
        }
        double obj_coeff = 0.0f;
        if (which_formulation_ == 3)
            obj_coeff = react->obj_coeff();
        DEBUG (
            'L', "      Adding cts vars for reaction " << k <<
            " " << react->name() <<
            " obj-coeff " << obj_coeff << 
            " lb " << lb << " ub " << ub
        );
        imp_->make_flux_vars (react, react->name(), obj_coeff, lb, ub);
    }

    // Add use vars for each reaction
    for (size_t k = 0; k < num_reactions(); k++) {
        auto react = reaction(k);
        double lb = 0.0f;
        double ub = 1.0f;
        if (!react->is_selected()) 
            ub = 0.0f;
            
        double obj_coeff = react->obj_coeff();
        if (react->is_biomass()) {
            obj_coeff = 0.0f; // Already counted above
        }
        else if (which_formulation_ == 3)
            obj_coeff *= form_3_mult;
            
        string name = react->name() + "-use";
        DEBUG (
            'L', "      Adding binary var for reaction " << k <<
            " " << name << " obj-coeff " << obj_coeff
        );
        
        imp_->make_use_var (react, name, obj_coeff, lb, ub);
    }

    for (size_t exp = 0; exp < num_experiments(); exp++) {
        for (size_t m = 0; m < num_metabolites(); m++) {
            auto met = metabolite(m);
            
            DEBUG (
                'L', "    Add constraints for " << met->name() <<
                " exp " << experiment(exp)->name()
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
            DEBUG (
                'L', "      Adding constraint for met " << m <<
                " " << met->name() <<
                " lb " << experiment(exp)->lb(met) <<
                " ub " << experiment(exp)->ub(met)
            );
            
            imp_->add_met_constraint (
                met, exp, experiment(exp)->lb(met), experiment(exp)->ub(met),
                reacts, coeffs
            );
        }
    }
    
    for (size_t k = 0; k < num_reactions(); k++) {
        auto react = scenario_->reaction(k);
        imp_->add_react_link_constraints (react);
    }
    
    DEBUG ('L', "Int formulation completed");
    imp_->finalise_formulation();
}

/**
 * @brief Solves the integer linear programming problem
 *
 * Performs the complete solving workflow:
 * 1. Calculates biomass multiplier for formulation 1
 * 2. Formulates the ILP problem
 * 3. Writes model to file if debug mode enabled
 * 4. Optimizes using the LP solver implementation
 * 5. Fixes rounding issues in the solution
 * 6. Creates and returns solution object
 *
 * @return SolutionPtr Pointer to solution object, or nullptr if infeasible or no continuous solution
 */
SolutionPtr
IntSolver::solve ()
{
    on_screen ("Solve with ILP solver");
    DEBUG ('g', "    Solving ILP with " << num_reactions() << " reactions");

    if (cts_sol_ == nullptr)
        return nullptr;
    
    if (which_formulation_ == 1)
        calc_biomass_mult();

    formulate();
        
    on_screen (
        "  Solving LP with " +
        to_string (imp_->num_vars()) + " vars and " +
        to_string (imp_->num_constraints()) + " constraints" +
        " biomass obj mult " + to_string (biomass_obj_mult_)
    );

    if (Debug::doDebug('w')) {
        string fn = "intsolver" + to_string(which_formulation_) + ".lp";
        DEBUG ('w', "Writing " << fn);
        imp_->write_model (fn);
        on_screen ("  Wrote " + fn);
    }

    status_ = imp_->optimize();

    if (status_ == INFEASIBLE_) {
        DEBUG ('w', "Bad _status: " << status_);
        return nullptr;
    }

    mip_gap_ = imp_->get_mip_gap();

    on_screen (
        "    Finished. Objective is " + to_string (imp_->get_objective()) + 
        " status " + to_string (status_)
    );
    DEBUG (
        'g', "Solved. Objective " << imp_->get_objective() <<
        " status " << status_
    );

    fix_rounding ();

    return imp_->make_sol();
}

/**
 * @brief Calculates the biomass objective coefficient multiplier for formulation 1
 *
 * Computes a multiplier to scale the biomass reaction objective coefficient to ensure
 * it dominates the objective function. The multiplier is based on the absolute objective
 * value from the continuous solution or parameter settings, and is further scaled by
 * a heuristic integer factor to avoid slow phase transition zones.
 *
 * Uses the maximum multiplier required across all experiments.
 */
void
IntSolver::calc_biomass_mult()
{
    if (cts_sol_ == nullptr && limeIsZero (params_->abs_obj_ub)) {
        DEBUG ('g', "No sol - can't do anything");
        return; 
    }
    
    biomass_obj_mult_ = params_->init_biomass_obj_mult;
    // Use the max mult required
    for (size_t exp = 0; exp < num_experiments(); exp++) {
        double abs_val = params_->abs_obj_ub;
        if (cts_sol_ != nullptr)
            abs_val = cts_sol_->sol(exp)->abs_obj_value();
        double base_mult = abs(scenario_->biomass_react()->obj_coeff());
        double this_mult = abs_val / base_mult;
        DEBUG (
            'g', "Exp " << exp << " abs obj val is " << abs_val <<
            " base mult is " << base_mult <<
            " this mult is " << this_mult
        );
        if (this_mult > biomass_obj_mult_)
            biomass_obj_mult_ = this_mult;
    }
    DEBUG (
        'g', "Max multiplier is " << biomass_obj_mult_ <<
        " after int factor applied " <<
        biomass_obj_mult_ * params_->biomass_mult_int
    );
    // HEURISTIC! Multiplying by a bit seems to take it out of the slow
    // "phase transition" zone
    biomass_obj_mult_ *= params_->biomass_mult_int;
}


/**
 * @brief Fixes rounding errors in the integer solution
 *
 * Post-processes the integer solution to handle fractional binary variables that
 * should be 0 or 1. This method is currently commented out but would iteratively
 * identify reactions with non-zero flux but fractional use variables, then test
 * forcing them up or down to find the best integer-feasible solution.
 *
 * The commented implementation uses a heuristic approach that:
 * - Identifies reactions with fractional use variables
 * - Tests forcing each variable up (to 1) or down (to 0)
 * - Selects the forcing direction that yields better objective and flux
 */
void
IntSolver::fix_rounding ()
{
/*
    Rand rand(0);
    int uniq_id = rand.uniform0n_1(LIME_BIG_INT);
    
    double obj_val = model_.get(_DoubleAttr_ObjVal);
    
    rounding_iters_ = 0;
    if (Debug::doDebug('w')) {
        string model_fn = limeFormat ("intsolver-%d.lp", rounding_iters_);
        DEBUG ('w', "Writing " << model_fn);
        model_.write (model_fn);
        string tmp_sol_fn =
            limeFormat (
                "intsolver-%s-%d.sol", "init", rounding_iters_
            );
        model_.write (tmp_sol_fn);
        DEBUG ('w', "Save sol in " << tmp_sol_fn);
    }
    
    while (true) {
        DEBUG ('g', "Check rounding");

        // Calc dummy flow
        // Choose biggesst non-zero flux, sorting dummies higher
        Sorter<size_t, double> sorter;
        for (size_t k = 0; k < num_reactions(); k++) {
            if (reaction(k)->is_selected()) {
                double u = use_[k].get(_DoubleAttr_X);
                DEBUG ('G', "  Reaction " << k << " u is " << u);
                if (u < 0.5) {
                    // Not used
                    double f = flux_[k].get(_DoubleAttr_X);
                    if (!limeIsZero (f)) {
                        double ub = flux_[k].get(_DoubleAttr_UB);
                        DEBUG (
                            'g', "  Reaction " << k << " " << *reaction(k) <<
                            " is 'not used' but has non-zero flux " << f <<
                            " ub is " << ub
                        );
                        //
                        if (ub < 0.5) {
                           // Already set UB to 0. Ignore for now.
                            continue;
                        }
                        if (reaction(k)->is_dummy())
                            f += 1.0; // All fluxes should be < 1e-6
                        sorter.add (k, -f);
                    }
                }
            }
        }
        if (sorter.size() == 0) {
            DEBUG ('g', "All good");
            break;
        }
        on_screen ("Found " + to_string (sorter.size()) + " fractional vars");
        sorter.doSort();
        size_t idx = sorter[0];
        auto react = reaction(idx);
        
        double u = use_[idx].get(_DoubleAttr_X);
        double force_up_obj_est =
            obj_val + (1.0f - u) * react->obj_coeff();
        on_screen ("Resolve LP with updated UB for " + react->name());
        DEBUG ('g', "  Test reaction idx " << idx << " " << react->name());

        rounding_iters_++;
        // Force up
        on_screen ("    Force up");
        DEBUG ('g', "    Test force up ");
        use_[idx].set(_DoubleAttr_LB, 1.0f);
        
        model_.update();

        if (Debug::doDebug('w')) {
            string model_fn =
                limeFormat ("intsolver-%s-%d.lp", "up", rounding_iters_);
            DEBUG ('w', "  Write model in " << model_fn);
            model_.write (model_fn);
        }
        
        model_.optimize();

        // Save sol for later use in re-solve"
        string sol_fn =
            limeFormat (
                "/tmp/intsolver-%d.sol", uniq_id
            );
        model_.write (sol_fn);
        DEBUG ('g', "Save sol in " << sol_fn);
        
        if (Debug::doDebug('w')) {
            string tmp_sol_fn =
                limeFormat (
                    "intsolver-%s-%d.sol", "up", rounding_iters_
                );
            model_.write (tmp_sol_fn);
            DEBUG ('w', "  Write sol in " << tmp_sol_fn);
        }
        
        DEBUG ('g', "      Force-up est is " << force_up_obj_est);
        double force_up_obj = model_.get(_DoubleAttr_ObjVal);
        on_screen ("      Force-up obj is " + to_string (force_up_obj));
        DEBUG ('g', "      Force-up obj is " << force_up_obj);
        double force_up_dummy = dummy_flux();
        
        // Force down
        on_screen ("    Force down");
        use_[idx].set(_DoubleAttr_LB, 0.0f);
        use_[idx].set(_DoubleAttr_UB, 0.0f);
        flux_[idx].set(_DoubleAttr_UB, 0.0f);
        model_.update();
        
        if (Debug::doDebug('w')) {
            string model_fn =
                limeFormat ("intsolver-%s-%d.lp", "down", rounding_iters_);
            DEBUG ('w', "  Write model in " << model_fn);
            model_.write (model_fn);
        }
        
        model_.optimize();

        if (Debug::doDebug('w')) {
            string tmp_sol_fn =
                limeFormat (
                    "intsolver-%s-%d.sol", "down", rounding_iters_
                );
            DEBUG ('w', "  Write sol in " << tmp_sol_fn);
            model_.write (tmp_sol_fn);
        }

        enum {NO_RESOLVE, REREAD_AND_RESOLVE, RESOLVE_ONLY} resolve =
                                                                NO_RESOLVE;
        
        int _status = model_.get(_IntAttr_Status);
        if (_status == _INFEASIBLE) {
            // No sol - must be included
            DEBUG ('g', "    Model is infeasible - lock in force up");
            on_screen ("    Model is infeasible - lock in force up");
            use_[idx].set(_DoubleAttr_LB, 1.0f);
            use_[idx].set(_DoubleAttr_UB, 1.0f);
            flux_[idx].set(_DoubleAttr_UB, react->flux_ub());
            resolve = REREAD_AND_RESOLVE;
        }
        else {
            // Recalc dummy flux
            double force_down_obj = model_.get(_DoubleAttr_ObjVal);
            double force_down_dummy = dummy_flux();
            DEBUG ('g', "    Force down obj is " << force_down_obj);
            on_screen ("      Force down obj is " + to_string (force_down_obj));
            
            // Favour force-up, since it will not lead to infeasibilities
            DEBUG ('g', "      Up dummy flux: " << force_up_dummy);
            DEBUG ('g', "    Down dummy flux: " << force_down_dummy);
            
            if (
                force_up_obj < force_down_obj + 0.1 ||
                force_down_dummy > force_up_dummy
            ) {
                // Lock in force-up
                on_screen ("      Lock in force up");
                DEBUG ('g', "      Lock in force up");
                use_[idx].set(_DoubleAttr_LB, 1.0f);
                use_[idx].set(_DoubleAttr_UB, 1.0f);
                flux_[idx].set(_DoubleAttr_UB, react->flux_ub());
                resolve = REREAD_AND_RESOLVE;
            }
            else {
                on_screen ("      Lock in force down");
                DEBUG ('g', "      Lock in force down");
                // UBs already set
            }
        }
        if (resolve == REREAD_AND_RESOLVE || resolve == RESOLVE_ONLY) {
            DEBUG ('g', "    Resolve to restore solution");
            on_screen ("    Resolve to restore solution");
            model_.update();
            if (resolve == REREAD_AND_RESOLVE)
                model_.read (sol_fn);
            model_.optimize();
            
            on_screen (
                "      obj is " +
                to_string (model_.get(_DoubleAttr_ObjVal))
            );
            DEBUG ('g', "      obj is " << model_.get(_DoubleAttr_ObjVal));
            if (Debug::doDebug('w')) {
                string model_fn =
                    limeFormat ("intsolver-resolve-%d.lp", rounding_iters_);
                DEBUG ('w', "  Writing " << model_fn);
                model_.write (model_fn);
            }
        }
        remove (sol_fn.c_str());
    }
    
    if (Debug::doDebug('w')) {
        string fn = "intsolver-end.lp";
        DEBUG ('w', "Writing " << fn);
        model_.write (fn);
        on_screen ("  Wrote " + fn);
    }
*/
}

/**
 * @brief Generates a summary string of the solver results
 *
 * @return std::string Summary containing biomass objective multiplier, MIP gap,
 *         and number of rounding iterations
 */
string
IntSolver::summary()
{
    stringstream str;
    str <<
        " biomass_obj_mult " <<
        fixed << setprecision(1) << biomass_obj_mult_ <<
        " mip_gap " << mip_gap_ << 
        " rounding_iters " << rounding_iters_;
    return str.str();
}
