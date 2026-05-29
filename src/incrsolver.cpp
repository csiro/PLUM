#include <sstream>
#include <list>
#include <iomanip>

#include "lime/debug.h"
#include "lime/numutil.h"
#include "lime/error.h"
#include "lime/biaschoice.h"
#include "lime/sortpairt.h"

#include "mosh/constants.h"
#include "mosh/incrsolver.h"
#include "mosh/intsolver.h"

using namespace std;
using namespace lime;
using namespace mosh;

void
IncrSolver::init()
{
    // We want to propagate in cost order
    SortPairT<Reaction*> to_propagate;
    
    if (selected_react_file_.is_open())
        selected_react_file_ << "# Init" << endl;
    for (auto& react : scenario_->reactions()) {
        // Enable (but don't propagate) dummy reactions
        if (react->is_dummy()) {
            react->set_selected(true);
        }
        else if (react->obj_coeff() <= min_cost_) {
            // Bootstrap with reactions up to min cost
            to_propagate.add (react.get(), react->obj_coeff());
        }
        else {
            react->set_selected(false);
            costs_.insert (react->obj_coeff());
        }
    }
    DEBUG ('i', "Propagating init");
    for (size_t k = 0; k < to_propagate.size(); k++) {
        auto react = to_propagate[k];
        DEBUG (
            'i', "  Propagate " << react->name() <<
            " val " << to_propagate.getKey(k)
        );
        react->set_selected(true);
        propagate (react);
    }
        
    DEBUG ('i', "Reaction costs are ");
    for (auto val : costs_) {
        DEBUG ('i', "  " << val);
    }
    scenario_->calc_sources();
    if (selected_react_file_.is_open())
        selected_react_file_ << "# Sources" << endl;
    for (auto& met : scenario_->metabolites()) {
        if (met->is_source())  {
            DEBUG ('i', "    Met " << *met << " is source");
            available_[met->index()] = true;
            if (selected_react_file_.is_open())
                selected_react_file_ << "MET " << met->name() << endl;
        }
    }
}

SolutionPtr
IncrSolver::solve ()
{
    on_screen ("Solving with incremental solver");
    DEBUG ('i', "Solve with incremental solver");

    if (num_experiments() != 1)
        limeCrash ("Incr solver only works on one experiment");
        
    status_ = INFEASIBLE_;
    
    lp_solver_.set_quiet(true);
    
    // Set up the formulation with only initial reactions selected
    lp_solver_.formulate (scenario_->experiment(0));

    // Check sol
    SolutionPtr sol = do_solve();
    
    while (sol == nullptr && !costs_.empty()) {
        cost_iters_++;
        // Add reactions with the next-lowest cost
        auto first = costs_.begin();
        max_used_cost_ = *first;
        costs_.erase(first);

        sol = cost_iter (max_used_cost_);
    }
    
#ifndef NDEBUG
    DEBUG ('i', "Available mets:");
    for (size_t k = 0; k < scenario_->num_metabolites(); k++) {
        if (available_[k]) {
            DEBUG ('i', "  " << scenario_->metabolite(k)->name());
        }
    }
    DEBUG ('i', "The following reactions are NOT enabled");
    for (auto& react : scenario_->reactions())
        if (!react->is_selected())
            DEBUG (
                'i', "  " << react->name() << " cost " << react->obj_coeff()
            );
#endif
    if (return_int_sol_ && sol != nullptr) {
        assert (imp2_ != nullptr);
        MultiSolPtr cts_sol = make_shared<MultiSol> (sol);
        // This will only use the selected reactions
        IntSolver int_solver (
            scenario_, params_, rand_.generateSeed(), 0, imp2_, cts_sol
        );
        sol = int_solver.solve();
        summary_ = int_solver.summary();
        status_ = int_solver.status();
    }
    
    return sol;
}

SolutionPtr
IncrSolver::cost_iter (double val)
{
    DEBUG ('i', "Adding reactions with obj-coeff " << val);
    on_screen ("  Add reactions with obj-coeff " + to_string(val));
    if (selected_react_file_.is_open())
        selected_react_file_ << "# Cost " << val << endl;

    int level = 0;
    
    while (true) {
        if (selected_react_file_.is_open())
            selected_react_file_ << "# Level " << level << endl;
        if (level_iter (val, level)) { // Add reactions at this level
            auto sol = do_solve();
            if (sol != nullptr) {
                DEBUG ('i', "Found sol at level: " << level);
                return sol;
            }
        }
        else if (!add_cycle_met (val)) { // Try adding required mets
            break; // No waiting mets left
        }
        level++;
        level_iters_++;
    }
    DEBUG ('i', "No sol for cost " << val);
    on_screen ("  No sol for obj-coeff " + to_string(val));
    return nullptr;
}

SolutionPtr
IncrSolver::do_solve ()
{
    DEBUG ('i', "Solve LP");

    DEBUG ('i', "  Reset biomass mult");
    lp_solver_.set_biomass_obj_mult (0, params_->init_biomass_obj_mult);
    
    auto sol = lp_solver_.do_solve();
    if (
        sol != nullptr &&
        sol->biomass_flux() > 0 &&
        limeIsZero (sol->dummy_flux()) 
    ) {
        DEBUG ('i', "Found sol!");
        status_ = lp_solver_.status();
        return sol;
    }
    return nullptr;
}

// Are the reactants for this reaction available?
bool
IncrSolver::is_enabled (const Reaction* react)
{
    for (auto met : react->in_mets()) {
        if (!is_available(met))
            return false;
    }
    return true;
}

// Do the products of this reaction include a carbon source?
bool
IncrSolver::makes_carbon (const Reaction* react)
{
    // Allowed to produce carbon if it uses one of the current
    // carbon sources to do it
    bool makes_c = false;
    for (auto met : react->out_mets()) {
        if (met->is_c_source() && !met->is_source()) {
            DEBUG (
                'i', "    React " << *react <<
                " makes c-src " << met->name()
            );
            makes_c = true;
        }
    }
    if (makes_c) {
        // Check if it gets its carbon from a valid source
        for (auto met : react->in_mets()) {
            if (met->is_source())  {
                DEBUG ('i', "    Input met " << *met << " is source, so OK");
                makes_c = false;
            }
        }
    }
    return makes_c;
}

// Make the products of this reaction available
void
IncrSolver::propagate (const Reaction* react)
{
    if (selected_react_file_.is_open())
        selected_react_file_ << "REACT " << react->name() <<
            " " << react->obj_coeff() << endl;
    for (auto met : react->out_mets()) {
        if (!is_available(met)) {
            DEBUG ('i', "    Met " << *met << " is available");
            available_[met->index()] = true;
            if (selected_react_file_.is_open())
                selected_react_file_ << "MET " << met->name() << endl;
        }
    }
}


bool
IncrSolver::level_iter (double cost, int level)
{
    DEBUG ('i', "  Adding reactions at level " << level);

    int count = 0;
    for (auto& react_ptr : scenario_->reactions()) {
        auto react = react_ptr.get();
        if (!react->is_selected()) {
            if (react->obj_coeff() <= cost) {
                DEBUG (
                    'i', "    Consider reaction " << *react <<
                    " cost " << react->obj_coeff()
                );
                if (
                    is_enabled (react) &&
                    !makes_carbon (react)
                ) {
                    DEBUG (
                        'i', "      " << *react <<
                        " enabled, cost " << react->obj_coeff()
                    );
                    count++;
                    react_ptr->set_selected(true);
                    lp_solver_.enable_reaction (react);
                    // Set up for next iter
                    propagate (react);
                }
            }
        }
    }
    on_screen (
        "    Adding " + to_string (count) +
        " reactions at level " + to_string (level) +
        " cost " + to_string(cost)
    );

    DEBUG ('i', "    Added " << count << " reactions");
    return count > 0;
}

// Add a metabolite that someone is waiting for, that can be supplied
// by someone else (i.e. that is in a cycle)
bool
IncrSolver::add_cycle_met (double cost)
{
    DEBUG ('i', "    Choose waiting");
    set<size_t> waiting;  // Indices of mets we are waiting for
    set<size_t> supplied; // Indices of mets that can be supplied at curr cost

    // See what is waiting, and what can be supplied
    BiasChoice choice (rand_.generateSeed());
    for (auto& react_ptr : scenario_->reactions()) {
        auto react = react_ptr.get();
        if (!react->is_selected()) {
            if (react->obj_coeff() <= cost) {
                for (auto met : react->in_mets()) {
                    if (!is_available (met)) {
                        DEBUG (
                            'j', "        React " << react->name() <<
                            " is waiting for " << *met
                        );
                        choice.addChoice(met->index());
                    }
                }
                for (auto met : react->out_mets()) {
                    if (!is_available (met))
                        supplied.insert(met->index());
                }
            }
        }
    }

    DEBUG ('i', "    Choose from " << choice.numChoices() << " mets");
    DEBUG ('i', "    Supplied has " << supplied.size() << " entries");
    DEBUG ('J', "Supplied:");
    for (auto k : supplied) {
        DEBUG('J', "  " << scenario_->metabolite(k)->name());
    }
    if (selected_react_file_.is_open())
        selected_react_file_ << "# Add cycle mets" << endl;
    bool ok = false;
    while(choice.numChoices() > 0) {
        int idx = choice.choose();
        auto met = scenario_->metabolite(idx);
        if (
            met->is_cycle_met() &&
            supplied.find(idx) != supplied.end()
        ) {
            // It can be supplied - we have our winner
            DEBUG (
                'i', "    Choose met " << *met
            );
            available_[idx] = true;
            ok = true;
            on_screen ("      Add cycle-break met " + met->name());
            if (selected_react_file_.is_open())
                selected_react_file_ << "MET " << met->name() << endl;
            cycle_mets_++;
        }
        if (!met->is_cycle_met()) { 
            DEBUG (
                'j', "        Not allowed to start a cycle " << *met <<
                " - try again"
            );
        }
        else { 
            DEBUG (
                'j', "        Can't supply " << *met << " - try again"
            );
        }
        // Can't be supplied - delete it
        choice.removeChoice (idx);
    }
    return ok;
}


string
IncrSolver::summary()
{
    stringstream str;
    str << " mincost " << min_cost_ <<
        " max_used " << max_used_cost_ <<
        " cost_iters " << cost_iters_ <<
        " level_iters " << level_iters_ <<
        " cycle_mets " << cycle_mets_;
    return str.str() + lp_solver_.summary() + summary_;
}
