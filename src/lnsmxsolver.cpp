/* LNS style multi-experiment solver

   - Reactions with lowest cost are assumed to be "in the model"
   - Find the set of reactions to match a set of experiments.

 */

 
#include <list>
#include <set>
#include <chrono>
#include <thread>
#include <math.h>

#include "lime/debug.h"
#include "lime/numutil.h"
#include "lime/strutil.h"
#include "lime/fileutil.h"
#include "lime/constants.h"
#include "lime/sortpairt.h"
#include "lime/error.h"

#include "mosh/lnsmxsolver.h"

#ifdef PLUM_LPSOLVER
#include "mosh/lpslpsolverimp.h"
#endif
#ifdef PLUM_HIGHS
#include "mosh/highslpsolverimp.h"
#endif

using namespace std;
using namespace lime;
using namespace mosh;

LnsMxSolver::LnsMxSolver (
    Scenario* scenario, Params* params, Flavour flavour, WhichObj which_obj,
    int seed, int max_iters, bool make_all_unit, string progname
) :
    GapSolver (scenario, params, seed),
    flavour_(flavour),
    which_obj_(which_obj),
    obj_elt_({COST, TAU, ERROR}),
    grbenv_(nullptr),
    imp_(scenario->num_experiments()),
    lp_solver_(scenario->num_experiments()),
    simanneal_(
        params->sa_targ_accept, params->sa_targ_prob, 0.05f, rand_.generateSeed()
    ),
    pareto_front_(),
    best_sol_(nullptr),
    incumb_sol_(nullptr),
    progname_(progname),
    save_best_fn_(""),
    save_best_sol_fn_(""),
    progress_(nullptr),
    max_iters_(max_iters),
    make_all_unit_(make_all_unit), 
    quiet_(1),
    biomass_obj_mult_(0.0f),
    base_cost_(0.0f),
    iter_(0),
    inner_iters_(0),
    num_improves_(0),
    num_simanneal_(0),
    iter_found_best_(0),
    target_flux_(scenario->num_experiments()),
    react_chooser_(rand_.generateSeed()),
    exp_chooser_(rand_.generateSeed()),
    adapt_choice_(
        (int)NUM_METHODS,
        AdaptCfg (
            params->adapt_sigma1, params->adapt_sigma2, params->adapt_sigma3,
            params->adapt_learn_rate, 100, true
        ),
        rand_.generateSeed()
    ), 
    tabu_list_(scenario->num_reactions(), never),
    target_delta_(0.0f),
    interval_len_((max_iters_ / params->num_intervals) + 1),
    unit_cost_(params->gene_ind_cost),
    non_unit_react_(scenario_->num_reactions(), false),
    bad_reacts_(),
    unit_fails_(),
    why_failed_(NO_FAIL),
    rank0_flux_(0.0f)
{
    DEBUG ('x', "LnsMxSolver");

    if (flavour == GUROBI) {
        try {
            grbenv_ = make_shared<GRBEnv>();
        }
        catch (GRBException e) {
            limeCrash (
                "Caught Gurobi exception " << e.getErrorCode() << ": " <<
                e.getMessage()
            );
        }
    }
    // Initialise solvers, one for each experiment
    for (size_t k = 0; k < num_exp(); k++) {
        switch (flavour) {
        case GUROBI:
            imp_[k] =
                make_shared<GrbLPSolverImp> (
                    scenario_, params, rand_.generateSeed(), grbenv_
                );
            break;
        case LP_SOLVE:
#ifdef PLUM_LPSOLVER
            imp_[k] =
                make_shared<LpsLPSolverImp> (scenario_, params);
#else
            limeCrash ("LPSOLVER not compiled in");
#endif
            break;
        case HIGHS:
#ifdef PLUM_HIGHS
            imp_[k] =
                make_shared<HighsLPSolverImp> (
                    scenario_, params, rand_.generateSeed()
                );
#else
            limeCrash ("HIGHS solver not compiled in");
#endif
        }
            
        lp_solver_[k] =
            make_shared<LPSolver> (
                scenario_, params, rand_.generateSeed(), 1, imp_[k]
            );
        lp_solver_[k]->set_which_experiment(k);
    }

    double max_biolog = scenario_->biolog_rank(0)->biolog_score();
#ifndef NDEBUG
    DEBUG ('x', "Ranking");
    for (size_t k = 0; k < num_exp(); k++) {
        auto exp = scenario_->biolog_rank(k);
        DEBUG (
            'x', "  Exp " << exp->index() << " " << exp->name() <<
            " " << exp->biolog_score() <<
            " is_growth " << exp->is_growth()
        );
    }
#endif

    // Set up select method chooser, with chance of selection in [0,1] 
    adapt_choice_.setWeight (     BIAS_SEL, 20.0f / 250.0f);
    adapt_choice_.setWeight (     COST_SEL, 20.0f / 250.0f);
    adapt_choice_.setWeight (     FLUX_SEL, 20.0f / 250.0f);
    adapt_choice_.setWeight (COST_FLUX_SEL, 20.0f / 250.0f);
    adapt_choice_.setWeight (     RAND_SEL, 20.0f / 250.0f);
    // Don't choose bad-sel randomly
    adapt_choice_.setWeight (      BAD_SEL,  0.0f / 250.0f); 
    // Make add choices sum to same as select choices
    adapt_choice_.setWeight (     ADD_COST, 50.0f / 250.0f);
    adapt_choice_.setWeight (     ADD_RAND, 50.0f / 250.0f);
    // And add a bit for making unit costs
    adapt_choice_.setWeight (    MAKE_UNIT, 50.0f / 250.0f);
}

// Are two doubles close enough, relatively speaking
bool rel_equal (double a, double b, double close)
{
    double base = limeMax (fabs(a), fabs(b));
    return (fabs (a - b) <= close * base);
}

SolutionPtr 
LnsMxSolver::solve ()
{
    if (quiet_ == 2)
        cout << "Initial solve" << endl;
    
    init_solve();

    if (quiet_ == 2) {
        cout << "LNS solve" << endl;
        cout <<
            limeFormat (
                "Initial sol    obj %6.2lf"
                " cost %6d tau %5.3lf error %6.2lf"
                " mis-growth%2d runaways%2d nonunit %4d",
                best_sol_->obj(), (int) best_sol_->cost(),
                best_sol_->tau(), best_sol_->error(),
                best_sol_->growth_mismatch(),
                best_sol_->num_runaways(),
                best_sol_->non_unit_used()
            ) << endl;
    }

    if (make_all_unit_)
        make_all_unit_cost();

    return improve();
}

SolutionPtr 
LnsMxSolver::improve ()
{
    int num_unselected = 0;
    int iters_no_improve = 0;

    if (quiet_ == 2)
        cout << "Start improvement" << endl;
    
    while (
        iter_ < max_iters_ &&
        (params_->time_limit == 0 || has_time_left()) 
    ) {
        iter_++;
        DEBUG ('x', "Iteration " << iter_);

        if (quiet_ == 2) {
            // 'quiet', but occasionally show prog
            if (iter_ % 100 == 0)
                cout <<
                    limeFormat (
                        "Iter %4d best obj %6.2lf"
                        " cost %6d tau %5.3lf error %6.2lf"
                        " mis-growth%2d runaways%2d nonunit %4d",
                        iter_, best_sol_->obj(), (int) best_sol_->cost(),
                        best_sol_->tau(), best_sol_->error(),
                        best_sol_->growth_mismatch(),
                        best_sol_->num_runaways(),
                        best_sol_->non_unit_used()
                    ) << endl;
        }
        else {
            screen_pos (0, 0);
            on_screen (progname_);
            on_screen (limeFormat ("Iteration %d/%d", iter_, max_iters_));
        }

        check_tabu_list (iter_);

        if ((iter_ % interval_len_) == 0) {
            // New interval. Recalc target flux
            // Calc target flux for each exp
            size_t rank0 = scenario_->biolog_rank0();
            double new_rank0_flux = target_flux_[rank0] - target_delta_;
            scenario_->calc_target_flux (target_flux_, new_rank0_flux);
        }

        int num_to_remove = 1;
        int num_to_add = 1;
        int num_to_change = 1;

        // List of reactions to remove, add, or make unit
        // Only ONE list will have elements
        ReactList remove_list;
        ReactList add_list;
        ReactList unit_list;
 
        Method method = NUM_METHODS;
        int count = 0;
        bool ok = false;
        do {
            iters_no_improve++;
            //remove_list.clear();
            //add_list.clear();
            //unit_list.clear();
            if (bad_reacts_.size() > 0) {
                // Try to get rid of the bad reactions
                method = BAD_SEL;
            }
            else if (incumb_sol_->non_unit_used() > 0 && count == 0) {
                // Try to get rid of the non-unit reactions
                // (first time round the loop, anyway)
                method = MAKE_UNIT;
            }
            else {
                method = (Method) adapt_choice_.suggest();
            }
            DEBUG (
                'x', "  Method is " << method_name_[method]
            );

            switch (method) {
            case BIAS_SEL:
            case COST_SEL:
            case FLUX_SEL:
            case COST_FLUX_SEL:
            case RAND_SEL:
            case BAD_SEL:
                if (
                    select_reactions_to_remove (
                        num_to_remove, method, remove_list
                    )
                ) {
                    ok = true;
                    disable_reactions (remove_list);
                }
                break;
            case ADD_COST:
            case ADD_RAND:
                if (
                    select_reactions_to_add (num_to_add, method, add_list)
                ) {
                    ok = true;
                    enable_reactions (add_list);
                }
                break;
            case MAKE_UNIT:
                if (
                    select_reactions_to_make_unit (
                        num_to_change, method, unit_list
                    )
                ) {
                    ok = true;
                    set_unit_react_cost (unit_list);
                }
                break;
            }
            if (count++ > 10) {
                // Looping - try removing a tabu
                if (!remove_tabu())  {
                    on_screen ("Looping, but no more tabu to remove");
                    break;
                }
            }
        }
        while (!ok);
        
        if (
            remove_list.size() == 0 &&
            add_list.size() == 0 &&
            unit_list.size() == 0
        ) {
            on_screen ("Nothing selected - give up");
            break;
        }
        adapt_choice_.use (method);
        string str = "   Add reaction";
        for (auto react : add_list)
            str += " " + react->name();
        on_screen (str);
        str = "Remove reaction";
        for (auto react : remove_list)
            str += " " + react->name();
        on_screen (str);
        str = "      Make unit";
        for (auto react : unit_list)
            str += " " + react->name();
        on_screen (str);
        
        enum Outcome {
            FAIL, NEW_BEST, NEW_BEST_BY_SECOND,
            NEW_INCUMB, NEW_INCUMB_BY_ADD, NEW_INCUMB_BY_SECOND, 
            IGNORE, RESTORE_BEST
        };
        vector<const char*> outcomeName = {
            "Fail", "New_best", "New_best_by_2nd",
            "New_incumb", "New_incumb_by_add", "New_incumb_by_2nd", 
            "Ignore", "Restore_best"
        };
        
        bool make_incumbant = false;
        Outcome outcome = FAIL;
        bool is_pareto = false;

        simanneal_.iter (iter_, best_sol_->obj());

        MultiSolPtr sol = do_solve();
        SolAttributesPtr curr_sol = nullptr;
        if (sol == nullptr) {
            DEBUG ('x', "do_solve failed");
        }
        else {
            // We found a new, feasible solution
            outcome = IGNORE;
            curr_sol =
                make_shared<SolAttributes> (this, sol, non_unit_react_);

            // "Equal" here is within 1%
            bool bestObj1Equal =
                rel_equal (curr_sol->obj(), best_sol_->obj(), 0.01f);
            bool incumbObj1Equal =
                rel_equal (curr_sol->obj(), incumb_sol_->obj(), 0.01f);
            
            DEBUG ('x', "New sol " << curr_sol->to_string());
            DEBUG (
                'x', "     best obj1 equal: " << bestObj1Equal <<
                " " << fabs (curr_sol->obj() - best_sol_->obj()) <<
                " " << curr_sol->obj() << " / " << best_sol_->obj()
            );
            DEBUG (
                'x', "   incumb obj1 equal: " << incumbObj1Equal <<
                " " << fabs (curr_sol->obj() - incumb_sol_->obj()) <<
                " " << curr_sol->obj() << " / " << incumb_sol_->obj()
            );

            is_pareto = pareto_front_.add (curr_sol->multi_obj(), curr_sol);
            DEBUG ('x', "Sol is on the pareto front");
            
            if (curr_sol->obj() < best_sol_->obj()) {
                DEBUG ('x', "Sol is new best :-)");
                save_best (curr_sol);
                make_incumbant = true;
                outcome = NEW_BEST;
                iters_no_improve = 0;
            }
            else if (bestObj1Equal && curr_sol->obj2() < best_sol_->obj2()) {
                DEBUG ('x', "New best by second");
                save_best (curr_sol);
                make_incumbant = true;
                outcome = NEW_BEST_BY_SECOND;
            }
            else if (bestObj1Equal && remove_list.size() == 0) {
                // Might be by add or by making unit
                DEBUG ('x', "New best by add (called by-2nd)");
                save_best (curr_sol);
                make_incumbant = true;
                outcome = NEW_BEST_BY_SECOND;
            }
            else if (incumbObj1Equal && remove_list.size() == 0) {
                DEBUG ('x', "Sol is same obj with no removes - keep it");
                make_incumbant = true;
                outcome = NEW_INCUMB_BY_ADD;
            }
            else if (
                incumbObj1Equal &&
                curr_sol->obj2() < incumb_sol_->obj2()
            ) {
                DEBUG ('x', "Sol has same obj but better obj2");
                make_incumbant = true;
                outcome = NEW_INCUMB_BY_SECOND;
            }
            else if (simanneal_.accept (curr_sol->obj(), incumb_sol_->obj())) {
                DEBUG ('x', "SimAnneal accept");
                make_incumbant = true;
                num_simanneal_++;
                outcome = NEW_INCUMB;
            }
        }
        if (
            params_->max_no_improve > 0 &&
            iters_no_improve > params_->max_no_improve
        ) {
            outcome = RESTORE_BEST;
        }

        if (sol == nullptr) {
            on_screen ("");
            on_screen (
                limeFormat (
                    "Iter %4d Incumb %8.2lf %s Obj %8s %s",
                    iter_, incumb_sol_->obj(), obj_name_[which_obj_],
                    "---", outcomeName[outcome]
                )
            );
            on_screen (
                limeFormat (
                    "          error %8s tau %5s cost %8s "
                    "mis-growth%4s runaways%3s nonunit %4s ", 
                    "---", "---", "---", "---", "--", "---"
                )
            );
            on_screen (
                limeFormat (
                    "          num unsel %4d improve %4d "
                    "simanneal %4d noimprove %4d ",
                    num_unselected, num_improves_, num_simanneal_,
                    iters_no_improve
                )
            );
        }
        else {
            on_screen ("");
            on_screen (
                limeFormat (
                    "Iter %4d Incumb %8.2lf %s Obj %8.2lf %s",
                    iter_, incumb_sol_->obj(),
                    obj_name_[which_obj_],
                    curr_sol->obj(), outcomeName[outcome]
                )
            );
            on_screen (
                limeFormat (
                    "          error %8.2lf tau %5.3lf cost %8d "
                    "mis-growth%4d runaways%3d nonunit %4d ", 
                    curr_sol->error(), curr_sol->tau(),
                    (int) curr_sol->cost(), curr_sol->growth_mismatch(),
                    curr_sol->num_runaways(), curr_sol->non_unit_used()
                )
            );
            on_screen (
                limeFormat (
                    "          num unsel %4d improve %4d "
                    "simanneal %4d noimprove %4d",
                    num_unselected, num_improves_, num_simanneal_,
                    iters_no_improve
                )
            );
            on_screen ("");
            on_screen (
                limeFormat (
                    "Iter %4d Incumb %8.2lf %s Obj %8.2lf %s",
                    iter_found_best_, best_sol_->obj(),
                    obj_name_[which_obj_],
                    best_sol_->obj(), "Best"
                )
            );
            on_screen (
                limeFormat (
                    "          error %8.2lf tau %5.3lf cost %8d "
                    "mis-growth%4d runaways%3d nonunit %4d ", 
                    best_sol_->error(), best_sol_->tau(),
                    (int) best_sol_->cost(), best_sol_->growth_mismatch(),
                    best_sol_->num_runaways(), best_sol_->non_unit_used()
                )
            );
            clear_screen (false);
        }
        
        if (progress_ != nullptr || Debug::doDebug('I')) {
            string name = "";
            double cost = 0.0f;
            double maxflux = 0.0f;
            ReactList* the_list = &add_list;
            size_t the_size = add_list.size();
            if (remove_list.size() > 0) {
                the_list = &remove_list;
                the_size = remove_list.size();
            }
            else if (unit_list.size() > 0) {
                the_list = &unit_list;
                the_size = unit_list.size();
            }
            for (auto react : *the_list) {
                name += react->name() + " ";
                cost += react->obj_coeff();
                if (incumb_sol_->sol() != nullptr)
                    maxflux += incumb_sol_->sol()->max_flux (react);
            }
            string sol_string =
                sol == nullptr ? string("FAILED") : curr_sol->to_string();
            stringstream prog;
            prog <<
                "Iter " << iter_ <<
                " method " << method_name_[method] << 
                " outcome " << outcomeName[outcome] << 
                " added " << add_list.size() << 
                " removed " << remove_list.size() <<
                " makeunit " << unit_list.size() <<
                " name " << name <<
                " cost " << cost <<
                " maxflux " << maxflux <<
                " rank0flux " << rank0_flux_ <<
                " ok " << (sol != nullptr) << 
                " whyfailed " << why_failed_ << 
                " " << sol_string << 
                " make_incumb " << make_incumbant <<
                " num_improve " << num_improves_ <<
                " simanneal " << num_simanneal_ <<
                " iters_no_improve " << iters_no_improve <<
                " target_flux " << target_flux_[scenario_->biolog_rank0()] <<
                " adaptchoice " << adapt_choice_ <<
                " " << simanneal_;
            
            DEBUG ('I', prog.str());
            if (progress_ != nullptr) 
                (*progress_) << prog.str() << endl;
        }
            
        if (outcome == RESTORE_BEST) {
            DEBUG ('x', "Restore best sol");
            incumb_sol_ = best_sol_;

            // First enable all reactions that _aren't_ in best_unselected_
            for (auto& react_ptr : scenario_->reactions()) {
                Reaction* react = react_ptr.get();
                if (!react->is_selected()) 
                    if (!list_contains (best_sol_->unselected(), react))
                        enable_reaction(react);
                // And clear any tabu on the reaction
                clear_tabu (react);
            }
            // Now disable reactions in best_unselected_
            for (auto react : best_sol_->unselected()) {
                if (react->is_selected()) 
                    disable_reaction(react);
            }
            // And reset costs on non-unit reacts
            reset_react_cost (best_sol_->non_unit_react());
            iters_no_improve = 0;
            bad_reacts_.clear();
        }
        else if (make_incumbant) {
            incumb_sol_ = curr_sol;
            
            // Don't consider re-adding/re-removing for a while
            DEBUG ('x', "Keep new sol as incumbant");
            DEBUG (
                'x', "  Mark new incumbant reactions as tabu for " <<
                params_->incumb_tabu_tenure << " iters"
            );
            if (remove_list.size() > 0) {
                make_tabu (
                    remove_list,
                    iter_ + params_->incumb_tabu_tenure / remove_list.size()
                );
            }
            make_tabu (add_list, iter_ + params_->incumb_tabu_tenure);
            
            if (bad_reacts_.size() > 0) {
                DEBUG ('x', "    Successfully removed bad reaction");
            }
        }
        else {
            DEBUG ('x', "Revert to incumbant");
            enable_reactions (remove_list);
            disable_reactions (add_list);
            reset_react_cost (unit_list);
            int tenure = 0;
            size_t divisor = 1;
            if (remove_list.size() > 1) {
                divisor = remove_list.size();
                DEBUG ('x', "    Divisor is " << divisor);
            }                
            if (sol == nullptr) {
                tenure = params_->fail_tabu_tenure;
                DEBUG (
                    'x', "  Mark failed reactions as tabu for " <<
                    tenure << " iters"
                );
            }
            else {
                tenure = params_->ignore_tabu_tenure / divisor;
                DEBUG (
                    'x', "  Mark unsuccessful reactions as tabu for " <<
                    tenure << " iters"
                );
            }
            make_tabu (remove_list, iter_ + tenure);
            make_tabu (add_list, iter_ + tenure);
            if (bad_reacts_.size() > 0) {
                DEBUG (
                    'x', "    Failed to remove bad reactions - revert costs"
                );
                reset_react_cost (remove_list);
            }
            if (unit_list.size() > 0) {
                // If I failed to make a reaction unit cost,
                // try to get rid of it
                DEBUG (
                    'x', "  Failed to make react unit cost - mark as bad"
                );
                bad_reacts_.insert (
                    bad_reacts_.begin(), unit_list.begin(), unit_list.end()
                );
            }
        }
        if (bad_reacts_.size() > 0) {
            // Remove reacts in 'remove_list' from bad_reacts
            for (auto react : remove_list) {
                auto it =
                    find (bad_reacts_.begin(), bad_reacts_.end(), react);
                if (it != bad_reacts_.end())
                    bad_reacts_.erase(it);
            }
        }
        num_unselected = incumb_sol_->unselected().size();

        // Update adaptation scores
        switch (outcome) {
        case FAIL:
        case RESTORE_BEST:
            // No update required
            break;
        case NEW_BEST:
        case NEW_BEST_BY_SECOND:
            adapt_choice_.score1(method);
            break;
        case NEW_INCUMB:
        case NEW_INCUMB_BY_ADD:
        case NEW_INCUMB_BY_SECOND:
            adapt_choice_.score2(method);
            break;
        case IGNORE:
            if (is_pareto)
                adapt_choice_.score3(method);
            break;
        }
    }

    // Reset target flux to that which gave the best result
    scenario_->calc_target_flux (target_flux_, best_sol_->rank0_target());

    on_screen ("Best sol: " + best_sol_->to_string());
    for (size_t ik = 0; ik < lp_solver_.size(); ik++) {
        size_t k = scenario_->biolog_rank(ik)->index();
        on_screen (
            limeFormat (
                "  Experiment %2d %10s target %6.4lf achieved %6.4lf "
                " reqd growth %d is growth %d is runaway %d", 
                k, experiment(k)->name().c_str(), 
                target_flux_[k], best_sol_->sol(k)->biomass_flux(),
                experiment(k)->is_growth(),
                params_->is_growth_flux(best_sol_->sol()->biomass_flux(k)),
                best_sol_->sol(k)->is_runaway()
            )
        );
    }
    
    return best_sol_->sol();
}

SolutionPtr 
LnsMxSolver::solve_by_add ()
{
    if (quiet_ == 2)
        cout << "Initial solve" << endl;
    
    init_solve();

    int num_unselected = 0;
    int iters_no_improve = 0;

    unselect_non_gene_ind(num_unselected);

    // We have just invalidated the best and incumb solutions. Reset them
    best_sol_ = nullptr;
    incumb_sol_ = nullptr;
    
    if (quiet_ == 2)
        cout << "Start iterations" << endl;

    // Form a list of reactions to add
    SortPairT<Reaction*> react_sorter;
    for (auto& react_ptr : scenario_->reactions()) {
        Reaction* react = react_ptr.get();
        if (
            !params_->is_gene_indicated(react) &&
            !react->is_selected() && 
            react->obj_coeff() < params_->max_react_cost
        ) {
            react_sorter.add (react, react->obj_coeff());
        }
    }
    react_sorter.doSort();
    ReactList reacts_to_add;
    for (size_t k = 0; k < react_sorter.size(); k++) {
        reacts_to_add.push_back (react_sorter[k]);
    }

    int add_iter = 0;
    int found_best = 0;
    size_t num_iters = reacts_to_add.size();
    while (
        reacts_to_add.size() > 0 &&
        (params_->time_limit == 0 || has_time_left()) 
    ) {
        add_iter++;
        DEBUG ('x', "Add Iteration " << add_iter);

        if (quiet_ == 2) {
            // 'quiet', but occasionally show prog
            if (add_iter % 100 == 0 && best_sol_ != nullptr)
                cout <<
                    limeFormat (
                        "Add iter %4d best obj %6.2lf"
                        " cost %6d tau %5.3lf error %6.2lf"
                        " mis-growth%2d runaways%2d nonunit %4d",
                        add_iter, best_sol_->obj(), (int) best_sol_->cost(),
                        best_sol_->tau(), best_sol_->error(),
                        best_sol_->growth_mismatch(),
                        best_sol_->num_runaways(),
                        best_sol_->non_unit_used()
                    ) << endl;
        }
        else {
            screen_pos (0, 0);
            on_screen (progname_);
            on_screen (limeFormat ("Add iteration %d/%d", add_iter, num_iters));
        }

        int num_to_add = 1;

        ReactList add_list;
        
        if (
            !select_cheap_reactions_to_add (reacts_to_add, num_to_add, add_list)
        ) {
            // No more reactions to add...
            break;
        }
        set_unit_react_cost (add_list);
        enable_reactions (add_list);

        string str = "   Add reaction";
        double cost = 0.0f;
        for (auto react : add_list) {
            str += " " + react->name();
            cost += react->obj_coeff();
        }
        on_screen (str);
        str = "           Cost " + limeFormat ("%7.2lf", cost);
        on_screen (str);
        
        enum Outcome {
            FAIL, NEW_BEST, NEW_BEST_BY_SECOND,
            NEW_INCUMB, NEW_INCUMB_BY_ADD, NEW_INCUMB_BY_SECOND, 
            IGNORE, RESTORE_BEST
        };
        vector<const char*> outcomeName = {
            "Fail", "New_best", "New_best_by_2nd",
            "New_incumb", "New_incumb_by_add", "New_incumb_by_2nd", 
            "Ignore", "Restore_best"
        };
        
        bool make_incumbant = false;
        Outcome outcome = FAIL;
        bool is_pareto = false;
        double best_obj = LIME_BIG_INT;
        double best_obj2 = LIME_BIG_INT;
        if (best_sol_ != nullptr) {
            best_obj = best_sol_->obj();
            best_obj2 = best_sol_->obj2();
        }
        double incumb_obj = LIME_BIG_INT;
        double incumb_obj2 = LIME_BIG_INT;
        if (incumb_sol_ != nullptr) {
            incumb_obj = incumb_sol_->obj();
            incumb_obj2 = incumb_sol_->obj2();
        }
        
        MultiSolPtr sol = do_solve();
        SolAttributesPtr curr_sol = nullptr;
        if (sol == nullptr) {
            DEBUG ('x', "do_solve failed");
            if (
                best_sol_ == nullptr &&
                why_failed_ == NO_GROWTH0
            ) {
                // Save it anyway until we find a solution
                make_incumbant = true;
                outcome = NEW_INCUMB_BY_ADD;
            }
        }
        else {
            // We found a new, feasible solution
            outcome = IGNORE;
            curr_sol =
                make_shared<SolAttributes> (this, sol, non_unit_react_);
            
            // "Equal" here is within .1%
            bool bestObj1Equal =
                rel_equal (curr_sol->obj(), best_obj, 0.001f);
            bool incumbObj1Equal =
                rel_equal (curr_sol->obj(), incumb_obj, 0.001f);
            
            DEBUG ('x', "New sol " << curr_sol->to_string());
            DEBUG (
                'x', "     best obj1 equal: " << bestObj1Equal <<
                " " << fabs (curr_sol->obj() - best_obj) <<
                " " << curr_sol->obj() << " / " << best_obj
            );
            DEBUG (
                'x', "   incumb obj1 equal: " << incumbObj1Equal <<
                " " << fabs (curr_sol->obj() - incumb_obj) <<
                " " << curr_sol->obj() << " / " << incumb_obj
            );

            is_pareto = pareto_front_.add (curr_sol->multi_obj(), curr_sol);
            DEBUG ('x', "Sol is on the pareto front");
            
            if (best_sol_ == nullptr || curr_sol->obj() < best_obj) {
                DEBUG ('x', "Sol is new best :-)");
                save_best (curr_sol);
                found_best = add_iter;
                make_incumbant = true;
                outcome = NEW_BEST;
                iters_no_improve = 0;
            }
            else if (bestObj1Equal && curr_sol->obj2() < best_obj2) {
                DEBUG ('x', "New best by second");
                save_best (curr_sol);
                found_best = add_iter;
                make_incumbant = true;
                outcome = NEW_BEST_BY_SECOND;
            }
            else if (bestObj1Equal) {
                // Has more reactions enabled - save it
                DEBUG ('x', "New best by add (called by-2nd)");
                save_best (curr_sol);
                found_best = add_iter;
                make_incumbant = true;
                outcome = NEW_BEST_BY_SECOND;
            }
            else if (
                incumbObj1Equal &&
                curr_sol->obj2() < incumb_obj2
            ) {
                DEBUG ('x', "Sol has same obj but better obj2");
                make_incumbant = true;
                outcome = NEW_INCUMB_BY_SECOND;
            }
        }

        if (sol == nullptr) {
            on_screen ("");
            on_screen (
                limeFormat (
                    "Iter %4d Incumb %8.2lf %s Obj %8s %s",
                    add_iter, incumb_obj, obj_name_[which_obj_],
                    "---", outcomeName[outcome]
                )
            );
            on_screen (
                limeFormat (
                    "          error %8s tau %5s cost %8s "
                    "mis-growth%4s runaways%3s nonunit %4s ", 
                    "---", "---", "---", "---", "--", "---"
                )
            );
            on_screen (
                limeFormat (
                    "          num unsel %4d improve %4d noimprove %4d ",
                    num_unselected, num_improves_, iters_no_improve
                )
            );
        }
        else {
            on_screen ("");
            on_screen (
                limeFormat (
                    "Iter %4d Incumb %8.2lf %s Obj %8.2lf %s",
                    add_iter, incumb_obj,
                    obj_name_[which_obj_],
                    curr_sol->obj(), outcomeName[outcome]
                )
            );
            on_screen (
                limeFormat (
                    "          error %8.2lf tau %5.3lf cost %8d "
                    "mis-growth%4d runaways%3d nonunit %4d ", 
                    curr_sol->error(), curr_sol->tau(),
                    (int) curr_sol->cost(), curr_sol->growth_mismatch(),
                    curr_sol->num_runaways(), curr_sol->non_unit_used()
                )
            );
            on_screen (
                limeFormat (
                    "          num unsel %4d improve %4d noimprove %4d",
                    num_unselected, num_improves_, iters_no_improve
                )
            );
        }
        if (best_sol_ != nullptr) {
            on_screen ("");
            on_screen (
                limeFormat (
                    "Iter %4d Incumb %8.2lf %s Obj %8.2lf %s",
                    found_best, best_sol_->obj(),
                    obj_name_[which_obj_],
                    best_sol_->obj(), "Best"
                )
            );
            on_screen (
                limeFormat (
                    "          error %8.2lf tau %5.3lf cost %8d "
                    "mis-growth%4d runaways%3d nonunit %4d ", 
                    best_sol_->error(), best_sol_->tau(),
                    (int) best_sol_->cost(), best_sol_->growth_mismatch(),
                    best_sol_->num_runaways(), best_sol_->non_unit_used()
                )
            );
            clear_screen (false);
        }
        
        if (progress_ != nullptr || Debug::doDebug('I')) {
            string name = "";
            double cost = 0.0f;
            for (auto react : add_list) {
                name += react->name() + " ";
                cost += react->obj_coeff();
            }
            string sol_string =
                sol == nullptr ? string("FAILED") : curr_sol->to_string();
            stringstream prog;
            prog <<
                "Iter " << add_iter <<
                " outcome " << outcomeName[outcome] << 
                " added " << add_list.size() << 
                " name " << name <<
                " cost " << cost <<
                " rank0flux " << rank0_flux_ <<
                " ok " << (sol != nullptr) << 
                " whyfailed " << why_failed_ << 
                " " << sol_string << 
                " make_incumb " << make_incumbant <<
                " num_improve " << num_improves_ <<
                " simanneal " << num_simanneal_ <<
                " iters_no_improve " << iters_no_improve <<
                " target_flux " << target_flux_[scenario_->biolog_rank0()] <<
                " " << simanneal_;
            
            DEBUG ('I', prog.str());
            if (progress_ != nullptr) 
                (*progress_) << prog.str() << endl;
        }
            
        if (make_incumbant) {
            incumb_sol_ = curr_sol;
            
            // Don't consider re-adding/re-removing for a while
            DEBUG ('x', "Keep new sol as incumbant");
            num_unselected--;
        }
        else {
            DEBUG ('x', "Revert to incumbant");
            disable_reactions (add_list);
            // num_unselected stays the same
        }
    }

    if (best_sol_ == nullptr) {
        on_screen ("No solution found");
    }
    else {
        // Reset target flux to that which gave the best result
        scenario_->calc_target_flux (target_flux_, best_sol_->rank0_target());

        on_screen ("Best sol: " + best_sol_->to_string());
        for (size_t ik = 0; ik < lp_solver_.size(); ik++) {
            size_t k = scenario_->biolog_rank(ik)->index();
            on_screen (
                limeFormat (
                    "  Experiment %2d %10s target %6.4lf achieved %6.4lf "
                    " reqd growth %d is growth %d is runaway %d", 
                    k, experiment(k)->name().c_str(), 
                    target_flux_[k], best_sol_->sol(k)->biomass_flux(),
                    experiment(k)->is_growth(),
                    params_->is_growth_flux(best_sol_->sol()->biomass_flux(k)),
                    best_sol_->sol(k)->is_runaway()
                )
            );
        }
    }
    if (best_sol_ == nullptr)
        return nullptr;

    return improve();
}

void
LnsMxSolver::unselect_non_gene_ind(int& num_unselected)
{
    DEBUG ('x', "Disable all non-gene-indicated reacts");
    for (auto& react_ptr : scenario_->reactions()) {
        auto react = react_ptr.get();
        if (
            react->is_selected() &&
            !params_->is_gene_indicated (react) &&
            !react->is_biomass()
        ) {
            disable_reaction (react);
            num_unselected++;
        }
    }
}


void
LnsMxSolver::make_all_unit_cost()
{
    if (quiet_ == 2)
        cout << "Make unit" << endl;
    DEBUG ('x', "Make all unit, num non-unit is " << num_non_unit());
    int local_iter = 0;
    int num_to_change = 1;
    int num_to_remove = 1;
    int fails = 0;
    int double_fails = 0;
    while (
        num_non_unit() > 0 &&
        (params_->time_limit == 0 || has_time_left()) 
    ) {
        local_iter++;
        DEBUG ('x', "Local Iteration " << local_iter);

        // Only ONE list will have elements
        ReactList remove_list;
        ReactList unit_list;
 
        Method method = MAKE_UNIT;
        if (bad_reacts_.size() > 0) {
            // Try to get rid of the bad reactions
            method = BAD_SEL;
        }
        DEBUG ('x', "  Method is " << method_name_[method]);

        if (method == MAKE_UNIT) {
            if (
                select_unused_reactions_to_make_unit (
                    num_to_change, unit_list
                )
            ) {
                set_unit_react_cost (unit_list);
            }
            else {
                // No more reactions to try
                break;
            }
        }
        else {
            if (
                select_reactions_to_remove (
                    num_to_remove, method, remove_list
                )
            ) {
                disable_reactions (remove_list);
            }
        }
        
        if (
            remove_list.size() == 0 &&
            unit_list.size() == 0
        ) {
            on_screen ("Nothing selected - give up");
            break;
        }
        string str = "Remove reaction";
        for (auto react : remove_list)
            str += " " + react->name();
        on_screen (str);
        str = "      Make unit";
        for (auto react : unit_list)
            str += " " + react->name();
        on_screen (str);
        
        bool keep = false;

        MultiSolPtr sol = do_solve();
        SolAttributesPtr curr_sol = nullptr;

        if (sol == nullptr) {
            DEBUG ('x', "do_solve failed");
        }
        else {
            curr_sol =
                make_shared<SolAttributes> (this, sol, non_unit_react_);

            // Flux is too low if it drops to less than 90% - but less than
            // 90% is OK if it decreases the number of runaways
            bool flux_too_low = 
                (
                    curr_sol->biomass_flux() <
                    0.9f * best_sol_->biomass_flux()
                ) &&
                (curr_sol->num_runaways() >= best_sol_->num_runaways());
            bool flux_too_high =
                curr_sol->biomass_flux() > 1.1f * best_sol_->biomass_flux();
            
            keep = !flux_too_low && !flux_too_high;
            
            DEBUG ('x', "New sol " << curr_sol->to_string());
            DEBUG (
                'x', "     keep " << keep <<
                " " << curr_sol->biomass_flux() << " / " <<
                best_sol_->biomass_flux()
            );
        }

        string upto_str =
            limeFormat (
                "Iter %4d Solved %3s Non-unit %5d Keep %3s Fails %5d Double %5d",
                local_iter, (curr_sol == nullptr ? "No": "Yes"), num_non_unit(),
                (keep ? "Yes" : "No"),
                fails, double_fails
            );

        if (quiet_ == 2) {
            // 'quiet', but occasionally show prog
            if (local_iter % 100 == 0)
                cout << upto_str << endl;
        }
        else {
            screen_pos (1, 0);
            on_screen (upto_str);
            clear_screen (false);
        }
        
        if (progress_ != nullptr || Debug::doDebug('I')) {
            const Reaction* react = nullptr;
            if (remove_list.size() > 0) {
                react = remove_list.front();
            }
            else if (unit_list.size() > 0) {
                react = unit_list.front();
            }
            assert (react != nullptr);
            string name = react->name();
            double cost = react->obj_coeff();
            int tabu = tabu_list_[react->index()];
            double biomass_flux = 0.0f;
            if (curr_sol != nullptr) 
                biomass_flux = curr_sol->biomass_flux();
            
            string sol_string =
                sol == nullptr ? string("FAILED") : curr_sol->to_string();
        
            stringstream prog;
            prog <<
                "UnitIter " << local_iter <<
                " method " << method_name_[method] << 
                " keep " << keep << 
                " added " << 0 << 
                " removed " << remove_list.size() <<
                " makeunit " << unit_list.size() <<
                " name " << name <<
                " cost " << cost <<
                " tabu " << tabu <<
                " biomass_flux " << biomass_flux <<
                " best_flux " << best_sol_->biomass_flux() <<
                " ok " << (sol != nullptr) << 
                " " << sol_string;

            DEBUG ('I', prog.str());
            if (progress_ != nullptr) 
                (*progress_) << prog.str() << endl;
        }
            
        if (keep) {
            DEBUG ('x', "Keep new sol as best");
            best_sol_ = curr_sol;
          
            if (bad_reacts_.size() > 0) {
                DEBUG ('x', "    Successfully removed bad reaction");
                unit_fails_.insert (
                    unit_fails_.end(), remove_list.begin(), remove_list.end()
                );
                // if it comes back, make sure it comes back as unit cost
                set_unit_react_cost (remove_list);
            }
        }
        else {
            DEBUG ('x', "Revert");
            enable_reactions (remove_list);
            reset_react_cost (unit_list);
            if (bad_reacts_.size() > 0) {
                DEBUG (
                    'x', "    Failed to remove bad reactions"
                );
                reset_react_cost (remove_list);
            }
            if (bad_reacts_.size() > 0) {
                double_fails++;
            }
            if (unit_list.size() > 0) {
                fails++;
                // If I failed to make a reaction unit cost,
                // try to get rid of it
                DEBUG (
                    'x', "  Failed to make react unit cost - mark as bad"
                );
                bad_reacts_.insert (
                    bad_reacts_.begin(), unit_list.begin(), unit_list.end()
                );
                // Don't try to remove again (or when we start normal iters)
                make_tabu (unit_list, iter_ + params_->incumb_tabu_tenure);
            }
        }
        if (bad_reacts_.size() > 0) {
            // Remove reacts in 'remove_list' from bad_reacts
            for (auto react : remove_list) {
                auto it =
                    find (bad_reacts_.begin(), bad_reacts_.end(), react);
                if (it != bad_reacts_.end())
                    bad_reacts_.erase(it);
            }
        }
    }
    incumb_sol_ = best_sol_;
    if (best_sol_ != nullptr && best_sol_->sol() != nullptr) {
        // Make sure it is not zero
        base_cost_ = best_sol_->sol()->abs_obj_value() + 1.0f; 
        DEBUG ('x', "New base cost is " << base_cost_);
    }
    // Reset non_unit
    DEBUG ('x', "Finished make-unit");
    for (size_t k = 0; k < num_reactions(); k++) {
        if (non_unit_react_[k]) {
            auto react = scenario_->reaction(k);
            DEBUG (
                'x', "  " << react->name() <<
                " cost " << react->obj_coeff() << 
                " is still non-unit"
            );
        }
    }
}

void
LnsMxSolver::init_solve()
{
    DEBUG ('x', "Init solve");

    clear_screen ();
    on_screen (progname_);
    on_screen (string("Objective is ") + obj_name_[which_obj_]);
    on_screen ("");
    on_screen ("Initial solve");

    auto sol = solve_all(true);
    
    // Calc target flux for each exp
    double target_flux0 = scenario_->biolog_rank(0)->base_flux();
    if (!limeIsZero(params_->target_flux))
        target_flux0 = params_->target_flux;
    DEBUG (
        'x', "Using rank 0 target flux " << target_flux0 <<
        " params target_flux is " << params_->target_flux
    );
    scenario_->calc_target_flux (target_flux_, target_flux0);
    for (size_t k = 0; k < num_exp(); k++) {
        lp_solver_[k]->set_target_flux(target_flux_[k]);
    }

    // What is the delta to get down to lower_target_mult * the current max flux
    // in mx_num_intervals steps.
    target_delta_ =
        target_flux_[scenario_->biolog_rank0()] *
        (1.0 - params_->lower_target_mult) /
        params_->num_intervals;

    // Note all non-unit reactions (of reasonable cost)
    for (auto& react_ptr : scenario_->reactions()) {
        auto react = react_ptr.get();
        if (
            react->obj_coeff() > unit_cost_ &&
            !params_->is_gene_indicated (react) &&
            !params_->exceeds_max_cost(react) 
        ) {
            non_unit_react_[react->index()] = true;
        }
    }

    // Don't do any searching
    params_->max_biomass_search_iters = 0;
    // Double the obj mult, to account for more expensive reactions
    // potentially being used during search.
    biomass_obj_mult_ = 0.0f;
    for (size_t k = 0; k < num_exp(); k++) {
        auto exp = experiment(k);

        if (exp->is_growth()) {
            SolutionPtr this_sol = sol->sol(k);

            limeAssert (!limeIsZero(target_flux_[k]));
            double this_mult =
                lp_solver_[k]->biomass_obj_mult() *
                params_->mx_biomass_inflation;
            
            DEBUG (
                'x', "  Exp " << k <<
                " raising mult to " << this_mult
            );
            lp_solver_[k]->set_biomass_obj_mult (
                k, this_mult
            );

            // Save max mult for no-growth experiments
            if (this_mult > biomass_obj_mult_)
                biomass_obj_mult_ = this_mult;
        }
    }
    // Set no-growth mults to generic biomass_mult_
    for (size_t k = 0; k < num_exp(); k++) {
        auto exp = experiment(k);

        if (!exp->is_growth()) {
            SolutionPtr this_sol = sol->sol(k);
            DEBUG (
                'x', "  Exp " << k <<
                " gets default biomass mult " << biomass_obj_mult_
            );
            lp_solver_[k]->set_biomass_obj_mult (k, biomass_obj_mult_);
        }
    }
    
    // Set base cost - used in objective calculations
    base_cost_ = sol->abs_obj_value() + 1.0f; // Make sure it is not zero
    DEBUG ('x', "Base cost is " << base_cost_);
    
    // Set up our current and best solutions
    incumb_sol_ = make_shared<SolAttributes> (this, sol, non_unit_react_);

    DEBUG (
        'x', "Init obj is " << incumb_sol_->obj() <<
        " error " << incumb_sol_->error() << " tau " << incumb_sol_->tau() <<
        " cost " << incumb_sol_->cost()
    );
    
    // Save initial solution as best
    save_best (incumb_sol_);

    // Set up simulated annealing
    simanneal_.init (incumb_sol_->obj(), max_iters_, params_->sa_restarts);
        
    if (progress_ != nullptr) {
        (*progress_) <<
            "Iter 0 " << 
            incumb_sol_->to_string() << 
            " adaptchoice " << adapt_choice_ <<
            " " << simanneal_ <<
            endl;
    }
}

MultiSolPtr
LnsMxSolver::solve_all (bool first_time)
{
    auto sol = std::make_shared<MultiSol> (scenario_, params_);

    bool ok = true;
    for (size_t k = 0; k < lp_solver_.size(); k++) {
        DEBUG (
            'Y', "        Solve for " << experiment(k)->name()
        );
        on_screen_nl (
            limeFormat ("Solve for %10s ", experiment(k)->name().c_str())
        );
        lp_solver_[k]->set_quiet(true);
        SolutionPtr this_sol = nullptr;
        if (first_time)
            this_sol = lp_solver_[k]->solve();
        else
            this_sol = lp_solver_[k]->do_solve();
        
        sol->set_sol(k, this_sol);
        
        if (this_sol == nullptr) {
            limeCrash ("No sol for experiment " << experiment(k)->name());
        }
        double flux = this_sol->biomass_flux();
        DEBUG (
            'Y', "          Found sol, flux " << flux <<
            " dummy flux " << this_sol->dummy_flux()
        );
        string msg = 
            limeFormat (
                "Found sol, biomass flux %7.3lf biolog target %7.3lf",
                flux, experiment(k)->biolog_score()
            );
        
        // Check if this is runaway flux
        if (params_->is_runaway_flux (flux)) {
            DEBUG ('Y', "            Runaway flux");
            msg += " Runaway reaction";
        }
        else if (first_time) {
            // Set up the base flux for this experiment
            experiment(k)->set_base_flux(flux);
        }
        
        on_screen (msg);
    }
    return sol;
}

// Run all solvers
// 
MultiSolPtr
LnsMxSolver::do_solve ()
{
    MultiSolPtr sol = std::make_shared<MultiSol> (scenario_, params_);
    
    DEBUG ('x', "Do solve");
    bool ok = true;
    int mismatch_count = 0;
    int runaway_count = 0;
    why_failed_ = NO_FAIL;
    for (size_t ik = 0; ik < lp_solver_.size(); ik++) {
        size_t k = scenario_->biolog_rank(ik)->index();
        
        DEBUG ('x', "  Inner iter, exp " << k << " " << experiment(k)->name());
        
        double best_flux = 0.0f;
        if (best_sol_ != nullptr && best_sol_->sol(k) != nullptr)
            best_flux = best_sol_->sol(k)->biomass_flux();
        on_screen_nl (
            limeFormat (
                "%10s biomass target %7.3lf best %7.3lf ",
                experiment(k)->name().c_str(), target_flux_[k],
                best_flux
            )
        );
        
        SolutionPtr this_sol = nullptr;
        string message2 = "";
        if (!ok) {
            DEBUG ('x', "    Skip");
            // Make empyty sol
            this_sol = make_shared<Solution> (scenario_, params_);
            message2 = " (skip)";
        }
        else {
            inner_iters_++;
            this_sol = lp_solver_[k]->do_solve();
        }
        
        sol->set_sol(k, this_sol);
        if (this_sol == nullptr) {
            DEBUG ('x', "    Solve failed");
            why_failed_ = NO_SOL;
            ok = false;
        }
        else {
            DEBUG (
                'x', "          Found sol, flux " << this_sol->biomass_flux() <<
                " dummy flux " << this_sol->dummy_flux()
            );
        }
        
        double flux = 0.0f;
        if (ok) {
            flux = this_sol->biomass_flux();
        }
       if (ik == 0)
            rank0_flux_ = flux;
            
        string message1 = limeFormat ("flux %7.3lf ", flux);
        
        // Check if this is runaway flux
        if (params_->is_runaway_flux (flux)) {
            message2 += " Runaway";
            runaway_count++;
            if (best_sol_ != nullptr) {
                DEBUG (
                    'x', "            Runaway flux " << runaway_count <<
                    " / " << best_sol_->num_runaways()
                );
                if (runaway_count > best_sol_->num_runaways()) {
                    why_failed_ = RUNAWAYS;
                    ok = false;
                }
            }
        }
        if (ik == 0 && !params_->is_growth_flux (flux)) {
            DEBUG ('x', "  No growth in rank0 experiment - give up");
            ok = false;
            why_failed_ = NO_GROWTH0;
            message2 += " No growth 0";
        }
        else if (
            ik == 0 &&
            (flux > params_->excess_growth_mult * target_flux_[k]) &&
            (
                best_sol_ == nullptr ||
                flux >
                params_->excess_growth_mult * best_sol_->sol(k)->biomass_flux()
            )
        ) {
            DEBUG (
                'x', "  Excess growth in rank0 experiment: flux " <<
                flux << " target " << target_flux_[k] <<
                " max " << (params_->excess_growth_mult * target_flux_[k]) <<
                " or " <<
                (best_sol_ == nullptr
                 ? 9999.0f
                 : params_->excess_growth_mult *
                   best_sol_->sol(k)->biomass_flux()
                )
            );
            ok = false;
            why_failed_ = EXCESS_GROWTH0;
            message2 += " Excess growth 0";
        }
        else if (
            ok &&
            experiment(k)->is_growth() &&
            !params_->is_growth_flux (flux)
        ) {
            mismatch_count++;
            if (
                best_sol_ != nullptr &&
                mismatch_count > best_sol_->growth_mismatch()
            ) {
                DEBUG ('x', "  Too many growth mismatches - give up");
                why_failed_ = GROWTH_MISMATCH;
                ok = false;
                message2 += " Growth mismatch";
            }
        }
        on_screen (message1 + message2);
    }
    return ok ? sol : nullptr;
}

/** Calculate the Kendall Tau rank score, compared to
    the biolog ranking
*/
// Return sign of x, or sign of y if x == 0
double sign (double x, double y)
{
    if (limeIsZero (x)) {
        return sign (y, 1.0);
    }
    return x < 0.0 ? -1.0f : 1.0f; 
}

double
LnsMxSolver::calc_kendall_tau (MultiSolPtr sol) const
{
    int n = num_exp();
    
    double sum = 0.0f;
    for (int i = 0; i < n; i++) {
        DEBUG (
            'x', "  Sol " << i << " " << experiment(i)->name() <<
            " target " << target_flux_[i] <<
            " got " << sol->sol(i)->biomass_flux()
        );
        
        for (int j = i+1; j < n; j++) {
            double x =
                experiment(i)->biolog_score() - experiment(j)->biolog_score();
            double y =
                sol->biomass_flux(i) - sol->biomass_flux(j);
            sum += sign (x, y) * sign (y, x);
        }
    }
    double val = 2.0f * sum / (n * (n-1));
    DEBUG ('x', "Kendall Tau is " << val);
    return val;
}

double
LnsMxSolver::mean_squared_error (MultiSolPtr sol) const
{
    double sum = 0.0f;
    for (size_t k = 0; k < num_exp(); k++) {
        double diff = 10.0f * (sol->biomass_flux(k) - target_flux_[k]);
        sum += diff * diff;
        DEBUG (
            'y', "  Sol " << k << " " << experiment(k)->name() <<
            " diff " << diff
        );
    }
    return sum / num_exp();
}

int 
LnsMxSolver::num_used () const
{
    int count = 0;
    if (best_sol_ == nullptr)
        return 0;
    for (size_t k = 0; k < num_reactions(); k++) 
        if (best_sol_->sol()->uses_react(k) )
            count++;
    
    return count;
}

void
LnsMxSolver::enable_reaction (Reaction* react)
{
    DEBUG ('y', "  Enable reaction " << *react);
    react->set_selected(true);
    for (auto& solver : lp_solver_) {
        solver->enable_reaction (react);
    }
}

void
LnsMxSolver::disable_reaction (Reaction* react)
{
    react->set_selected(false);
    for (auto& solver : lp_solver_) {
        solver->disable_reaction (react);
    }
}

void
LnsMxSolver::set_react_cost (Reaction* react, double cost)
{
    for (size_t k = 0; k < num_exp(); k++) {
        lp_solver_[k]->set_react_cost (k, react, cost);
    }
}

void
LnsMxSolver::save_best (SolAttributesPtr sol)
{
    DEBUG ('x', "New best found at iter " << iter_ << " obj " << sol->obj());
    if (best_sol_ != nullptr)
        DEBUG ('x', "  Old best " << best_sol_->obj());
    
    best_sol_ = sol;
    
    DEBUG ('x', "  " << sol->to_string());

    if (Debug::doDebug('X')) {
        DEBUG ('X', "Used reactions: ");
        int count = 0;
        for (size_t k = 0; k < num_reactions(); k++) {
            if (sol->sol()->uses_react(k)) {
                count++;
                auto react = scenario_->reaction(k);
                DEBUG_NL (
                    'X', react->name() <<
                    " (" << sol->sol()->max_flux (react) << ") ";
                );
            }
        }
        DEBUG ('X', "\n  Count is " << count);
        DEBUG ('X', "Unselected:");
        for (auto react : sol->unselected()) {
            DEBUG_NL ('X', react->name() << " ");
        }
        DEBUG ('X', "\n  Count is " << sol->num_unselected());
        lp_solver_[0]->write_model("best.lp");
    }


    iter_found_best_ = iter_;
    num_improves_++;
    
    if (isFilename (save_best_fn_)) {
        // Rewrite unselected filename
        ofstream out (save_best_fn_);
        out << "# Produced by " << progname_ << " on " << todayString() << endl;
        write_best_unselected (out);
    }
    if (isFilename (save_best_sol_fn_)) {
        // Append to best sol file
        ofstream out (save_best_sol_fn_, ios::app);
        write_best_sol (out);
    }
}

void
LnsMxSolver::write_best_unselected (std::ostream& out)
{
    if (best_sol_ == nullptr)
        return;
    for (auto react : best_sol_->unselected())
        out << react->name() << endl;
}


void
LnsMxSolver::write_best_sol (std::ostream& out)
{
    if (best_sol_ == nullptr)
        return;
    
    out << "# -------------------------------------------------------" << endl;
    out << "# " << best_sol_->to_string() << endl;
    out << "# Target biomass" << endl;
    out << "# ";
    for (size_t k = 0; k < num_exp(); k++)
        out << " " << limeFormat ("%6.4lf", target_flux(k));
    out << endl;
    best_sol_->sol()->write_flux (out);
}


bool
LnsMxSolver::select_reactions_to_remove (
    int num_to_remove, Method method, ReactList& remove_list
)
{
    size_t target_exp = exp_chooser_.choose(); 
    if (method == BIAS_SEL) {
        // Set up for experiment seleceion in biased choice
        // Select an experiment to target.
        // Bias based on excess biomass
        DEBUG ('x', "Choosing experiment to target");
        exp_chooser_.clear();
        for (size_t k = 0; k < num_exp(); k++) {
            if (is_over_target (k, incumb_sol_->sol())) {
                DEBUG (
                    'x', "  Add choice " << k <<
                    " " << scenario_->experiment(k)->name() <<
                    " with val " <<
                    incumb_sol_->sol()->biomass_flux(k) - target_flux(k) <<
                    " flux " << incumb_sol_->sol()->biomass_flux(k) <<
                    " target " << target_flux(k)
                );
                exp_chooser_.addChoice (
                    k, incumb_sol_->sol()->biomass_flux(k) - target_flux(k)
                );
            }
        }
        if (exp_chooser_.numChoices() == 0) {
            DEBUG ('x', "  Everyone is underperforming. Use add");
            return false;
        }
    }
    
    react_chooser_.clear();
    DEBUG ('x', "  Using select method " << method_name_[method]);
    switch (method) {
    case BIAS_SEL: {
        size_t target_exp = exp_chooser_.choose(); 
        DEBUG ('x', "    Target experiment " << target_exp);
        remove_reactions_biased (target_exp);
        break;
    }
    case COST_SEL:
        remove_reactions_cost ();
        break;
    case FLUX_SEL:
        remove_reactions_flux ();
        break;
    case COST_FLUX_SEL:
        remove_reactions_cost_flux ();
        break;
    case RAND_SEL:
        remove_reactions_rand ();
        break;
    case BAD_SEL:
        remove_reactions_bad ();
        break;
    }
        
    DEBUG ('x', "Make del choices");
    
    while (
        react_chooser_.numChoices() > 0 &&
        remove_list.size() < num_to_remove
    ) {
        auto react = react_chooser_.choose();
        DEBUG (
            'x', "  Choose " << react->name() <<
            " cost " << react->obj_coeff() <<
            " score " << react_chooser_.weightFor (react)
        );
        react_chooser_.removeChoice (react);
        remove_list.push_back (react);
    }
    DEBUG (
        'x', "  Removing " << remove_list.size() << " reacts"
    );
    return remove_list.size() > 0;
}

bool
LnsMxSolver::select_reactions_to_add (
    int num_to_add, Method method, ReactList& add_list
)
{
    react_chooser_.clear();
    int tabu_count = 0;
    int unsel_count = 0;
    int maxcost_count = 0;
    
    for (auto& react_ptr : scenario_->reactions()) {
        Reaction* react = react_ptr.get();
        if (react->is_selected()) 
            continue;
        unsel_count++;
        
        if (is_tabu (react)) {
            tabu_count++;
            continue;
        }
        if (params_->exceeds_max_cost(react)) {
            maxcost_count++;
            continue;
        }
        if (params_->is_gene_indicated(react)) // Can't remove gene-indicated
            continue;

        assert (!react->is_biomass());
        
        double choice_val = 0.0f;
        switch (method) {
        case ADD_COST:
            choice_val = 1.0f/react->obj_coeff();
            break;
        case ADD_RAND:
            choice_val = 1.0f;
            break;
        }
        react_chooser_.addChoice (react, choice_val);
    }
    DEBUG ('x', "Make add choices");
    DEBUG (
        'x', "  " << tabu_count << " of " << unsel_count <<
        " unselected were tabu. " << maxcost_count <<
        " exceeded maxcost. Num choices is " <<
        react_chooser_.numChoices()
    );
    
    while (
        react_chooser_.numChoices() > 0 &&
        add_list.size() < num_to_add
    ) {
        auto react = react_chooser_.choose();
        DEBUG (
            'x', "  Choose " << react->name() <<
            " cost " << react->obj_coeff() <<
            " score " << react_chooser_.weightFor (react)
        );
        react_chooser_.removeChoice (react);
        add_list.push_back (react);
    }
    DEBUG (
        'x', "  Adding " << add_list.size() << " reacts"
    );
    return add_list.size() > 0;
}

bool
LnsMxSolver::select_cheap_reactions_to_add (
    ReactList& reacts_to_add,
    int num_to_add, ReactList& add_list
)
{
    DEBUG ('x', "Choosing cheap reactions from " << reacts_to_add.size() << " choices");
    react_chooser_.clear();

    if (reacts_to_add.empty())
        return false;
    double cheapest = reacts_to_add.front()->obj_coeff();
    DEBUG ('x', "  Cheapest react cost is " << cheapest);
    // Now add all reacts that are within 5% of the cheapest
    cheapest *= 1.05;
    DEBUG ('x', "         Selecting below " << cheapest);
    
    for (auto react : reacts_to_add) {
        if (react->obj_coeff() > cheapest)
            break;
        // Add with equal probability
        react_chooser_.addChoice (react, 1.0f);
    }
    DEBUG ('x', "  Make add choices");
    DEBUG ('x', "  Num choices is " << react_chooser_.numChoices());
    
    while (
        react_chooser_.numChoices() > 0 &&
        add_list.size() < num_to_add
    ) {
        auto react = react_chooser_.choose();
        DEBUG (
            'x', "  Choose " << react->name() <<
            " cost " << react->obj_coeff()
        );
        react_chooser_.removeChoice (react);
        add_list.push_back (react);
        reacts_to_add.remove (react);
    }
    DEBUG (
        'x', "  Adding " << add_list.size() << " reacts. Now have " <<
        reacts_to_add.size() << " reacts left to add"
    );
    return add_list.size() > 0;
}

bool
LnsMxSolver::select_reactions_to_make_unit (
    int num_to_change, Method method, ReactList& unit_list
)
{
    react_chooser_.clear();
    
    for (auto& react_ptr : scenario_->reactions()) {
        auto react = react_ptr.get();
        if (
            !react->is_selected() ||
            !non_unit_react_[react->index()] ||
            is_tabu (react) ||
            !incumb_sol_->sol()->uses_react(react)
        ) 
            continue;

        assert (!react->is_biomass());

        // Bias 
        double choice_val = 1.0f / react->obj_coeff();
        react_chooser_.addChoice (react, choice_val);
    }
    DEBUG ('x', "Make unit choices");
    DEBUG ('x', "   Num choices is " << react_chooser_.numChoices());
    
    while (
        react_chooser_.numChoices() > 0 &&
        unit_list.size() < num_to_change
    ) {
        auto react = react_chooser_.choose();
        DEBUG (
            'x', "  Choose " << react->name() <<
            " cost " << react->obj_coeff() <<
            " score " << react_chooser_.weightFor (react)
        );
        react_chooser_.removeChoice (react);
        unit_list.push_back (react);
    }
    DEBUG (
        'x', "  Make " << unit_list.size() << " reacts unit cost"
    );
    return unit_list.size() > 0;
}

bool
LnsMxSolver::select_unused_reactions_to_make_unit (
    int num_to_change, ReactList& unit_list
)
{
    react_chooser_.clear();
    
    for (auto& react_ptr : scenario_->reactions()) {
        auto react = react_ptr.get();
        if (
            !react->is_selected() ||
            !non_unit_react_[react->index()] ||
            is_tabu (react)
        ) 
            continue;
        
        assert (!react->is_biomass());

        // Bias 
        double choice_val = 1.0f / react->obj_coeff();
        react_chooser_.addChoice (react, choice_val);
    }
    DEBUG ('x', "Make unused react unit choices");
    DEBUG ('x', "   Num choices is " << react_chooser_.numChoices());
    
    while (
        react_chooser_.numChoices() > 0 &&
        unit_list.size() < num_to_change
    ) {
        auto react = react_chooser_.choose();
        DEBUG (
            'x', "  Choose " << react->name() <<
            " cost " << react->obj_coeff() <<
            " score " << react_chooser_.weightFor (react)
        );
        react_chooser_.removeChoice (react);
        unit_list.push_back (react);
    }
    DEBUG (
        'x', "  Make " << unit_list.size() << " reacts unit cost"
    );
    return unit_list.size() > 0;
}

/** Set up react_chooser_ with biased weights */
void
LnsMxSolver::remove_reactions_biased (size_t target_exp)
{
    DEBUG (
        'x', "  Make biased selection targetting " << target_exp <<
        " " << experiment(target_exp)->name()
    );
    // Choose reaction(s)
    // - Reaction should be used by chosen experiment
    // - Reduce bias (a lot) if used by underperforming experiments
    // - Reduce bias (a little) if used by on-target experiments
    // - Increase bias (a lot) if used by over-target experiments
    // - Bias manipulation based on proportion of experiments in
    //   those categories

    static constexpr double used_by_under_discount = 0.5f;
    static constexpr double used_by_on_discount = 0.25f;
    static constexpr double used_by_over_penalty = 2.0f;
    
    vector<bool> under_target (num_exp(), false);
    vector<bool> on_target (num_exp(), false);
    vector<bool> over_target (num_exp(), false);
    for (size_t k = 0; k < num_exp(); k++) {
        if (is_over_target (k, incumb_sol_->sol())) {
            DEBUG (
                'x', "  Experiment " << k <<
                " " << experiment(k)->name() << " is  over target" 
            );
            over_target[k] = true;
        }
        else if (is_on_target (k, incumb_sol_->sol())) {
            DEBUG (
                'x', "  Experiment " << k <<
                " " << experiment(k)->name() << " is    on target" 
            );
            on_target[k] = true;
        }
        else if (is_under_target (k, incumb_sol_->sol())) {
            DEBUG (
                'x', "  Experiment " << k <<
                " " << experiment(k)->name() << " is under target" 
            );
            under_target[k] = true;
        }
    }

    int react_count = 0;
    int tabu_count = 0;
    for (auto& react_ptr : scenario_->reactions()) {
        Reaction* react = react_ptr.get();

        if (incumb_sol_->sol()->uses_react(react))
            react_count++;
        if (incumb_sol_->sol()->uses_react(react) && is_tabu (react))
            tabu_count++;
        if (
            !incumb_sol_->sol(target_exp)->uses_react(react) ||
            is_tabu (react) ||
            react->is_biomass() ||
            params_->is_gene_indicated(react)
        )
            continue;

        int count_used_by_over_target = 0;
        int count_used_by_on_target = 0;
        int count_used_by_under_target = 0;
        
        for (size_t k = 0; k < num_exp(); k++) {
            if (!incumb_sol_->sol(k)->uses_react(react))
                continue;
            if (over_target[k])
                count_used_by_over_target++;
            else if (on_target[k])
                count_used_by_on_target++;
            else if (under_target[k])
                count_used_by_under_target++;
        }

        /*
        double prop_used_by_over_target =
            (double) count_used_by_over_target / num_exp();
        double prop_used_by_on_target =
            (double) count_used_by_on_target / num_exp();
        double prop_used_by_under_target =
            (double) count_used_by_under_target / num_exp();

        DEBUG ('y', "  Reaction " << react->name());
        DEBUG (
            'y', "     flux: " << incumb_sol_->sol(target_exp)->flux(react)
        );
        DEBUG (
            'y', "    Under: " << count_used_by_under_target << 
            " = prop " << prop_used_by_under_target 
        );
        DEBUG (
            'y', "       On: " << count_used_by_on_target << 
            " = prop " << prop_used_by_on_target 
        );
        DEBUG (
            'y', "     Over: " << count_used_by_over_target << 
            " = prop " << prop_used_by_over_target 
        );
        */
        
        /*
        double orig_react_score = react->obj_coeff();
        double react_score = orig_react_score;
        react_score -= 
            prop_used_by_under_target * used_by_under_discount *
            orig_react_score;
        react_score -= 
            prop_used_by_on_target * used_by_on_discount *
            orig_react_score;
        react_score += 
            prop_used_by_over_target * used_by_over_penalty *
            orig_react_score;
        
        DEBUG (
            'y', "    Orig score " << orig_react_score <<
            " used score " << react_score
        );
        */
        DEBUG (
            'y', "  Reaction " << react->name() <<
            " under " << count_used_by_under_target <<
            " on " << count_used_by_on_target <<
            " over " << count_used_by_over_target
        );
        if (count_used_by_over_target == 0) {
            DEBUG (
                'y', "    Not used by over-target - skip"
            );
        }
        else {
            double react_score = 1.0;
            if (count_used_by_under_target + count_used_by_on_target > 0) {
                react_score =
                    ((double) count_used_by_over_target) /
                    (count_used_by_under_target + count_used_by_on_target + 1);
            }
            else {
                react_score = (double) count_used_by_over_target;
            }
            DEBUG (
                'y', "    score is " << react_score
            );
            react_chooser_.addChoice (react, react_score);
        }
    }
    DEBUG (
        'x', "  " << tabu_count << " of " << react_count <<
        " used reactions are tabu. Num choices " << react_chooser_.numChoices()
    );
}

/** Set up react_chooser_ with weights based only on cost */
void
LnsMxSolver::remove_reactions_cost ()
{
    DEBUG ('x', "  Using remove by cost");
    int react_count = 0;
    int tabu_count = 0;
    for (auto& react_ptr : scenario_->reactions()) {
        Reaction* react = react_ptr.get();
        if (incumb_sol_->sol()->uses_react(react))
            react_count++;
        if (incumb_sol_->sol()->uses_react(react) && is_tabu (react))
            tabu_count++;
        if (
            !incumb_sol_->sol()->uses_react(react) ||
            is_tabu (react) ||
            react->is_biomass() ||
            params_->is_gene_indicated(react)
        )
            continue;
        react_chooser_.addChoice (react, react->obj_coeff());
    }
    DEBUG (
        'x', "  " << tabu_count << " of " << react_count <<
        " used reactions are tabu. Num choices " << react_chooser_.numChoices()
    );
}

/** Set up react_chooser_ with weights based on flux*/
void
LnsMxSolver::remove_reactions_flux ()
{
    DEBUG ('x', "  Using remove by flux");
    int react_count = 0;
    int tabu_count = 0;
    for (auto& react_ptr : scenario_->reactions()) {
        Reaction* react = react_ptr.get();
        if (incumb_sol_->sol()->uses_react(react))
            react_count++;
        if (incumb_sol_->sol()->uses_react(react) && is_tabu (react))
            tabu_count++;
        if (
            !incumb_sol_->sol()->uses_react(react) ||
            is_tabu (react) ||
            react->is_biomass() || 
            params_->is_gene_indicated(react)
        )
            continue;
        react_chooser_.addChoice (react, incumb_sol_->sol()->max_flux(react));
    }
    DEBUG (
        'x', "  " << tabu_count << " of " << react_count <<
        " used reactions are tabu. Num choices " << react_chooser_.numChoices()
    );
}

/** Set up react_chooser_ with weights based on cost * flux
     - we want the most expensive reaction that is used a lot
       this is useful particularly when we have runaway reactions.
 */
void
LnsMxSolver::remove_reactions_cost_flux ()
{
    DEBUG ('x', "  Using remove by cost * flux");
    int react_count = 0;
    int tabu_count = 0;
    for (auto& react_ptr : scenario_->reactions()) {
        Reaction* react = react_ptr.get();
        if (incumb_sol_->sol()->uses_react(react))
            react_count++;
        if (incumb_sol_->sol()->uses_react(react) && is_tabu (react))
            tabu_count++;
        if (
            !incumb_sol_->sol()->uses_react(react) ||
            is_tabu (react) ||
            react->is_biomass() ||
            params_->is_gene_indicated(react)
        )
            continue;
        
        react_chooser_.addChoice (
            react,
            incumb_sol_->sol()->max_flux(react) * react->obj_coeff() /
            scenario_->max_react_cost()
        );
    }
    DEBUG (
        'x', "  " << tabu_count << " of " << react_count <<
        " used reactions are tabu. Num choices " << react_chooser_.numChoices()
    );
}

/** Set up react_chooser_ with uniform weights */
void
LnsMxSolver::remove_reactions_rand ()
{
    DEBUG ('x', "  Using remove randomly");
    int react_count = 0;
    int tabu_count = 0;
    for (auto& react_ptr : scenario_->reactions()) {
        Reaction* react = react_ptr.get();
        if (incumb_sol_->sol()->uses_react(react))
            react_count++;
        if (incumb_sol_->sol()->uses_react(react) && is_tabu (react))
            tabu_count++;
        if (
            !incumb_sol_->sol()->uses_react(react) ||
            is_tabu (react) ||
            react->is_biomass() ||
            params_->is_gene_indicated(react)
        )
            continue;
        react_chooser_.addChoice (react, 1.0f);
    }
    DEBUG (
        'x', "  " << tabu_count << " of " << react_count <<
        " used reactions are tabu. Num choices " << react_chooser_.numChoices()
    );
}

/** Set up react_chooser_ from bad_reacts list */
void
LnsMxSolver::remove_reactions_bad ()
{
    DEBUG ('x', "  Using remove bad reacts");
    int react_count = 0;
    for (auto& react : bad_reacts_) {
        react_chooser_.addChoice (react, 1.0f);
    }
    DEBUG (
        'x', "  Num choices " << react_chooser_.numChoices()
    );
}

/** Remove the oldest tabu  */
bool
LnsMxSolver::remove_tabu ()
{
    DEBUG ('x', "Clear tabu status for a reaction");
    size_t max_tabu = 0;
    size_t max_at = 0;
    for (size_t k = 0; k < tabu_list_.size(); k++) {
        if (tabu_list_[k] > max_tabu && tabu_list_[k] != never) {
            max_tabu = tabu_list_[k];
            max_at = k;
        }
    }
    if (max_tabu == 0) { // no tabu entries found
        DEBUG ('x', "No tabu entries found");
        return false;
    }
    DEBUG (
        'x', "  Clearing tabu for react " << max_at <<
        " " << scenario_->reaction(max_at)->name() <<
        " was " << tabu_list_[max_at] <<
        " curr iter is " << iter_
    );
    tabu_list_[max_at] = never;
    return true;
}


void
LnsMxSolver::write_model(std::ostream& out, MultiSolPtr sol)
{
    // List of used reactions
    list<const Reaction*> reacts;
    // Unique mets
    set<const Metabolite*> mets;

    for (size_t k = 0; k < num_reactions(); k++) {
        auto react = scenario_->reaction(k);
        
        if (
            !sol->uses_react(k) &&
            !params_->is_gene_indicated (react)
        )
            continue; // No flux, and not gene-indicated. Not interested
        
        reacts.push_back (react);

        for (auto metidx : react->mets()) {
            auto met = scenario_->metabolite(metidx);

            mets.insert (met);
        }
    }

    for (auto met : mets) {
        met->write_to (out);
    }

    for (auto react : reacts) {
        react->write_to (out, react->obj_coeff());
    }
}

void
LnsMxSolver::write_model(std::ostream& out)
{
    if (best_sol_ == nullptr)
        return;
    write_model (out, best_sol_->sol());
}

void
LnsMxSolver::write_pareto (ostream& out, string pareto_fn, string header)
{
    for (size_t k = 0; k < pareto_front_.size(); k++) {
        auto& objvec = pareto_front_.obj(k);
        auto sol = pareto_front_.solution(k);
        size_t id = (k + 1);
        out <<
            "id " << id << 
            " " << sol->to_string() << 
            endl;
        
        
        if (isFilename (pareto_fn)) {
            string model_fn = pareto_fn + "-" + to_string(id) + ".pld";
            ofstream model (model_fn);
            model << header << endl;
            model << "# Pareto front element " << id << endl;
            write_model (model, sol->sol());
        }
    }
}

void
LnsMxSolver::write_descr(std::ostream& out)
{
    out << "Summary stats" << endl;
    if (best_sol_ == nullptr) {
        out << "NO SOLUTION" << endl;
        return;
    }
    out << limeFormat ("%20s: %8.2lf", "Best obj", best_sol_->obj()) << endl;
    out << limeFormat (
        "%20s: %8d", "Best growth mismatch", best_sol_->growth_mismatch()
    ) << endl;
    out << limeFormat (
        "%20s: %8.2lf", "Best cost", best_sol_->cost()
    ) << endl;
    out << limeFormat (
        "%20s: %8.2lf", "Best max cost", best_sol_->max_cost()
    ) << endl;
    out << limeFormat (
        "%20s: %8.2lf", "Best tau", best_sol_->tau()
    ) << endl;
    out << limeFormat (
        "%20s: %8.2lf", "Best error", best_sol_->error()
    ) <<
        endl;

    out << endl;
    out << "Experiment stats" << endl;
    for (size_t ik = 0; ik < lp_solver_.size(); ik++) {
        size_t k = scenario_->biolog_rank(ik)->index();
        out <<
            limeFormat (
                "  Experiment %2d %10s biolog %6.4lf target %6.4lf"
                " achieved %6.4lf reqd growth %d is growth %d biomass mult %10.2lf", 
                k, experiment(k)->name().c_str(),
                experiment(k)->biolog_score(),
                target_flux_[k], best_sol_->sol(k)->biomass_flux(),
                experiment(k)->is_growth(),
                params_->is_growth_flux(best_sol_->sol()->biomass_flux(k)),
                lp_solver_[k]->biomass_obj_mult()
            ) << endl;
    }

    out << endl;
    out << "Reactions that could not be made unit cost and and were initially removed" << endl;
    for (auto& react : unit_fails_) {
        out << "  " << react->name() <<
            " cost " << react->obj_coeff() <<
            " max flux " << best_sol_->sol()->max_flux (react) <<
            " full name \"" << react->full_name() <<
            "\" " << react->nice_formula() << endl;
    }

    out << endl;
    out << "Reactions that could not be made unit cost" << endl;
    for (auto& react_ptr : scenario_->reactions()) {
        Reaction* react = react_ptr.get();
        if (
            non_unit_react_[react->index()] &&
            best_sol_->sol()->uses_react(react)
        ) {
            out << "  " << react->name() <<
                " cost " << react->obj_coeff() <<
                " max flux " << best_sol_->sol()->max_flux (react) <<
                " full name \"" << react->full_name() <<
                "\" " << react->nice_formula() << endl;
        }
    }

    out << endl;
    out << "Used reactions" << endl;
    SortPairT<Reaction*> used;
    for (size_t k = 0; k < num_reactions(); k++) {
        if (!best_sol_->sol()->uses_react(k))
            continue;
        auto react = scenario_->reaction(k);
        used.add (react, react->obj_coeff());
    }
    used.doSort();
    for (size_t k = 0; k < used.size(); k++) {
        auto react = used.get(k);
        out << "  " << react->name() <<
            " cost " << react->obj_coeff() <<
            " max flux " << best_sol_->sol()->max_flux (react) <<
            " full name \"" << react->full_name() <<
            "\" " << react->nice_formula() << endl;
    }
}

