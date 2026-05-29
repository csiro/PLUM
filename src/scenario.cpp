
#include <sstream>
#include <algorithm>
#include <list>
#include <math.h>

#include "lime/debug.h"
#include "lime/linereader.h"
#include "lime/limetok.h"
#include "lime/numutil.h"
#include "lime/strutil.h"
#include "lime/sortpair.h"
#include "lime/constants.h"
#include "lime/error.h"
#include "lime/line.h"

#include "mosh/scenario.h"
#include "mosh/solution.h"
#include "mosh/multisol.h"

using namespace std;
using namespace lime;

namespace mosh
{
    bool 
    isDontCare (double val)
    {
        return (int) val == DONT_CARE;
    }
    
    bool
    isIgnoreVal (double val)
    {
        return (int) val == IGNORE;
    }
}
    
using namespace mosh;

void
Scenario::read_data (std::string data_fn)
{
    cout << "Reading " << data_fn << endl;
    DEBUG ('d', "Reading pld data " << data_fn);
    LineReader reader (data_fn);

    string line;
    while (reader.getLine (line)) {
        DEBUG ('d', "Line " << line);
        LimeTok tok (line);

        string tag = tok.nextString ();
        DEBUG ('d', "  Tag " << tag);

        if (tag == "MET") {
            bool error = false;
            string name = tok.nextString (error);
            if (error) 
                reader.error ("Bad format in MET line");
            string full_name = tok.nextString (error);
            
            // Make sure it is not already there...
            if (find_metabolite (name) != nullptr)
                cout << "Skip duplicate metabolite " << name << endl;
            else {
                add_metabolite (std::make_shared<Metabolite> (name, full_name));
            }
        }
        else if (tag == "REACTION") {
            bool error = false;
            string name = tok.nextString (error);
            double obj_coeff = tok.nextDouble (error);
            if (error) 
                reader.error ("Bad format in REACTION line");
            double flux_ub = tok.nextDouble (error);
            if (error) {
                // Probably old-format file. Use 1000 as default
                flux_ub = 1000.0f;
            }

            // Make sure it is not already there...
            if (find_reaction (name) != nullptr) {
                cout << "Skip duplicate reaction " << name << endl;
                continue;
            }
            
            if (obj_coeff < 0.0f) {
                // We have a neg value reaction - use flux maximisation
                cout << "Use flux max for reaction " << name << endl;
            }
            // Select for LP if it has zero cost
            bool is_active = !isIgnoreVal (obj_coeff);
            bool is_selected = is_active;
            auto reaction =
                make_shared<Reaction> (
                    name, obj_coeff, flux_ub, is_active, is_selected
                );
            add_reaction (reaction);
            DEBUG ('d', "    Add react " << reaction->name());
            
            string met_name = tok.nextString (error);
            while (!error) {
                if (met_name.compare ("__fullname__") == 0) {
                    string fullname = tok.nextString (error);
                    if (error) 
                        reader.error (
                            string("Error in REACTION (fullname) for react ") +
                            name
                        );
                    reaction->set_full_name (fullname);
                }
                else {
                    double met_coeff = tok.nextDouble (error);
                    if (error) 
                        reader.error (
                            string("Error in REACTION (met coeff) for react ") +
                            name
                        );
                    
                    auto met = find_metabolite (met_name);
                    if (met == nullptr) {
                        // Quietly add a new metabolite 
                        DEBUG ('d', "      Add unseen met " << met_name);
                        auto met_ptr = 
                            std::make_shared<Metabolite> (met_name, met_name);
                        add_metabolite (met_ptr);
                        met = met_ptr.get();
                    }
                    else {
                        DEBUG ('d', "      Found met " << met_name << " " << *met);
                    }
                    DEBUG (
                        'd', "      Add coeff " << met_coeff <<
                        " for " << met->name()
                    );
                    reaction->set_coeff (met, met_coeff);
                }
                met_name = tok.nextString (error);
            }
        }
        else
            reader.error (string ("Unknown tag: ") + tag);
    }
    for (auto& met : metabolite_)
        DEBUG ('d', "MET " << *met.get());
            
    cout << "Read " << metabolite_.size() << " metabolites" << endl;
    cout << "Read " << reaction_.size() << " reactions" << endl;

    orig_num_reactions_ = reaction_.size();
    select_all();
}

void
Scenario::finalise (Params* params)
{
    if (params->use_dummy)
        add_dummy_reactions (params->dummy_cost);
    else if (
        params->use_dummy_biomass_react_dm ||
        params->use_dummy_biomass_react_ex ||
        params->use_dummy_biomass_prod_dm
    )
        add_dummy_biomass_dm (params);
    
    if (params->max_react_cost > 0)
        unselect_above_cost (params->max_react_cost);
    
    if (params->preserve_dummies > 0)
        preserve_dummies ();
    
    if (params->unit_cost)
        set_react_cost (params->gene_ind_cost);

    if (biomass_react_ != nullptr && !limeIsZero (params->biomass_ub))
        biomass_react_->set_flux_ub (params->biomass_ub);
}

// Add dummay reactions for all mets.
void
Scenario::add_dummy_reactions (double dummy_cost)
{
    cout << "Adding dummy reactions with obj value " << dummy_cost << endl;
    // Add expensive reactions with singleton input for each metabolic 
    // residual, just so eqns can balance out
    for (size_t k = 0; k < num_metabolites(); k++) {
        auto met = metabolite(k);

        // Add a dummy with +ve coeff - dummy supply
        string name = met->name() + "-dummy-supply";
        double obj_coeff = dummy_cost;
        bool is_active = true;
        // Select col for LP (always make dummies available)
        bool is_selected = true;
        bool is_dummy = true;
        auto reaction =
            make_shared<Reaction> (
                name, obj_coeff, DUMMY_FLUX_UB,
                is_active, is_selected, is_dummy
            );
        add_reaction (reaction);
        
        double met_coeff = 1.0f;
        reaction->set_coeff (met, met_coeff);
        
        // Add a dummy with -ve coeff - dummy demand
        name = met->name() + "-dummy-demand";
        obj_coeff = dummy_cost;
        is_active = true;
        is_selected = true;
        is_dummy = true;
        reaction =
            make_shared<Reaction> (
                name, obj_coeff, DUMMY_FLUX_UB,
                is_active, is_selected, is_dummy
            );
        add_reaction (reaction);
        
        met_coeff = -1.0f;
        reaction->set_coeff (met, met_coeff);
    }
    cout << " Now " << reaction_.size() << " reactions with dummies" << endl;
}

// Select existing dummies, identified by their name
void
Scenario::preserve_dummies ()
{
    int count = 0;
    for (auto& react : reaction_) {
        if (
            !react->is_selected() &&
            (react->name().find("-dummy-") != string::npos)
        ) {
            react->set_selected(true);
            count++;
        }
    }

    cout << "Preserved " << count << " '-dummy-' reactions" << endl;
}

// Add dummy DM reactions for biomass reactants 
void
Scenario::add_dummy_biomass_dm (Params* params)
{
    if (biomass_react_ == nullptr)
        return;
    
    cout << "Adding dummy DM reactions for biomass (obj value " <<
        params->dummy_cost << ")" << endl;
    DEBUG (
        'a', "Adding dummy DM reactions for biomass, obj value " <<
        params->dummy_cost <<
        " react-dm: " << params->use_dummy_biomass_react_dm <<
        " react-ex: " << params->use_dummy_biomass_react_ex <<
        " prod-dm: " << params->use_dummy_biomass_prod_dm
    );

    if (params->use_dummy_biomass_react_dm) {
        for (auto reactant : biomass_react_->in_mets()) {
            // Add a dummy with -ve coeff - dummy demand - to mop up
            // any excess biomass reactants
            string name = reactant->name() + "-dummy-demand";
            double obj_coeff = params->dummy_cost;
            bool is_active = true;
            // Select col for LP (always make dummies available)
            bool is_selected = true;
            bool is_dummy = true;
            DEBUG (
                'a', "  Adding " << name << " for reactant " <<
                reactant->name()
            );
            auto reaction =
                make_shared<Reaction> (
                    name, obj_coeff, DUMMY_FLUX_UB,
                    is_active, is_selected, is_dummy
                );
            add_reaction (reaction);
        
            double met_coeff = -1.0f;
            reaction->set_coeff (reactant, met_coeff);
        }
    }
    if (params->use_dummy_biomass_react_ex) {
        for (auto reactant : biomass_react_->in_mets()) {
            // Add a dummy with +ve coeff - dummy supply - to supply 
            // biomass reactants
            string name = reactant->name() + "-dummy-supply";
            double obj_coeff = params->dummy_cost;
            bool is_active = true;
            // Select col for LP (always make dummies available)
            bool is_selected = true;
            bool is_dummy = true;
            DEBUG (
                'a', "  Adding " << name << " for reactant " <<
                reactant->name()
            );
            auto reaction =
                make_shared<Reaction> (
                    name, obj_coeff, DUMMY_FLUX_UB,
                    is_active, is_selected, is_dummy
                );
            add_reaction (reaction);
        
            double met_coeff = 1.0f;
            reaction->set_coeff (reactant, met_coeff);
        }
    }
    if (params->use_dummy_biomass_prod_dm) {
        for (auto prod : biomass_react_->out_mets()) {
            // Add a dummy with -ve coeff - dummy demand - to mop up
            // any excess biomass products
            string name = prod->name() + "-dummy-demand";
            double obj_coeff = params->dummy_cost;
            bool is_active = true;
            // Select col for LP (always make dummies available)
            bool is_selected = true;
            bool is_dummy = true;
            DEBUG (
                'a', "  Adding " << name << " for product " <<
                prod->name()
            );
            auto reaction =
                make_shared<Reaction> (
                    name, obj_coeff, DUMMY_FLUX_UB,
                    is_active, is_selected, is_dummy
                );
            add_reaction (reaction);
        
            double met_coeff = -1.0f;
            reaction->set_coeff (prod, met_coeff);
        }
    }

    cout << " Now " << reaction_.size() << " reactions with bio-dummies" <<
        endl;
}

void
Scenario::read_flux (std::string flux_fn)
{
    cout << "Reading flux " << flux_fn << endl;
    DEBUG ('d', "Reading flux " << flux_fn);
    LineReader reader (flux_fn);

    // When using a sol, make 0-flux reactions inactive
    for (auto& react : reaction_)
        react->set_active(false);

    // Each line is of the form
    // reaction-id flux [error-bound]
    // where
    //   reaction-id appears in reaction list
    //   flux is the estimated flux
    //  error-bound is a positive real, indicating the error bounds

    string line;
    while (reader.getLine (line)) {
        LimeTok tok (line);

        bool error = false;
        string react_str = tok.nextString (error);
        double flux = tok.nextDouble (error);
        if (error)
            reader.error (string ("Bad format in flux file"));
        if (flux < 0) {
            // Biomass reaction - skip
            continue;
        }
        double flux_error = tok.nextDouble (error);
        if (error) 
            flux_error = 0.0;
        
        if (isIgnoreVal (flux_error)) {
            DEBUG (
                'f', "  Reaction " << react_str <<
                " has uncertain flux - skip"
            );
            continue;
        }

        string name = react_str;
        auto react = find_reaction (name);
        if (react == nullptr) {
            limeWarning (string ("Lost reaction ") + name + " in flux file");
            continue;
        }
        
        DEBUG (
            'f', "  Set known flux for reaction " << react->name() <<
            " to " <<  flux
        );
        react->set_known_flux (flux, flux_error);
        react->set_active(true);
    }
}

// Read, updating according to policy.
// Note: those wil gene-indicated cost remain gene-indicated
void
Scenario::read_react_cost (
    string react_cost_fn, const Params* params, ReactCostPolicy policy
)
{
    cout << "Reading reaction cost file " << react_cost_fn << endl;
    DEBUG ('d', "Reading react cost " << react_cost_fn);
    LineReader reader (react_cost_fn);

    // Each line is of the form
    // reaction-id flux [error-bound]
    // where
    //   reaction-id appears in reaction list
    //   flux is the estimated flux
    //  error-bound is a positive real, indicating the error bounds

    string line;
    int count = 0;
    while (reader.getLine (line)) {
        LimeTok tok (line);

        bool error = false;
        string react_str = tok.nextString (error);
        double new_cost = tok.nextDouble (error);
        if (error)
            reader.error (string ("Bad format in react-cost file"));

        auto react = find_reaction (react_str);
        if (react == nullptr) {
            DEBUG (
                'f', "    Ignore reaction " << react_str <<
                " cost " <<  new_cost
            );
            // quietly ignore
            continue;
        }

        if (limeDblEqual (react->obj_coeff(), params->gene_ind_cost)) {
            DEBUG (
                'f', "    Ignore gene-indicated reaction " << react->name()
            );
            continue;
        }
        if (react->is_biomass()) {
            DEBUG ('f', "    Ignore biomass reaction " << react->name());
            continue;
        }

        double old_cost = react->obj_coeff();
        bool update = false;
        switch (policy)
        {
        case USEMIN:
            update = new_cost < old_cost;
            break; 
        case USEMAX:
            update = new_cost > old_cost;
            break; 
        case REPLACE:
            update = new_cost != old_cost;
            break;
        }
        if (update) {
            DEBUG (
                'f', "  Set obj coeff for reaction " << react->name() <<
                " to " <<  new_cost << " was " << old_cost
            );
            react->set_obj_coeff (new_cost);
            count++;
        }
        else {
            DEBUG (
                'f', "  Leave obj coeff for reaction " << react->name() <<
                " at " <<  old_cost << " ignore " << new_cost
            );
        }
    }
    cout << "Updated cost on " << count << " reactions" << endl;
}

void
Scenario::read_c_sources (string c_source_fn)
{
    DEBUG ('d', "Reading c sources " << c_source_fn);
    if (
        c_source_fn.length() == 0 ||
        c_source_fn.compare ("none") == 0
    )
        return;
            
    cout << "Reading carbon sources " << c_source_fn << endl;
    LineReader reader (c_source_fn);
            
    string line;
    while (reader.getLine (line)) {
        LimeTok tok (line);
                
        string name = tok.nextString ();
        auto met = (Metabolite*)find_metabolite (name);
        if (met == nullptr) {
            limeWarning (string ("Lost metabolite ") + name + " in carbon sources");
            continue;
        }
        met->set_c_source(true);
    }
}

void
Scenario::read_cycle_mets (string cycle_met_fn)
{
    DEBUG ('d', "Reading cycle mets " << cycle_met_fn);
    if (
        cycle_met_fn.length() == 0 ||
        cycle_met_fn.compare ("none") == 0
    )
        return;
            
    cout << "Reading cycle mets " << cycle_met_fn << endl;
    LineReader reader (cycle_met_fn);
            
    string line;
    while (reader.getLine (line)) {
        LimeTok tok (line);
                
        string name = tok.nextString ();
        auto met = (Metabolite*)find_metabolite (name);
        if (met == nullptr) {
            limeWarning (string ("Lost metabolite ") + name + " in cycle mets");
            continue;
        }
        met->set_cycle_met(true);
    }
}

void
Scenario::read_supply_demand (string supply_demand_fn, const Params* params)
{
    cout << "Reading supply/demand " << supply_demand_fn << endl;
    DEBUG ('d', "Reading supply/demand " << supply_demand_fn);
    LineReader reader (supply_demand_fn);

    // Each line is of the form
    // metabolite-id supply_demand 
    // where
    //   metabolite-id appears in metabolite list
    //   supply_demand is the estimated supply_demand -
    //      +ve if supplied in substrate
    //      -ve if produced by the reactions,

    DEBUG (
        'm', "Creating exp " << basename (supply_demand_fn) <<
        " from " << supply_demand_fn
    );
    ExperimentPtr experiment =
        make_shared<Experiment> (
            basename (supply_demand_fn), num_experiments(), num_metabolites()
        );
    experiment_.push_back (experiment);
    biolog_rank_.push_back (NULL);
    
    string line;
    while (reader.getLine (line)) {
        LimeTok tok (line);

        bool error = false;
        string met_str = tok.nextString (error);

        if (met_str == "biolog") {
            // This is the biolog reading for the experiment
            double biolog = tok.nextDouble (error);
            if (error)
                reader.error ("Bad format in biolog in supply_demand file");
            experiment->set_biolog_score (
                biolog, params->is_growth_biolog (biolog)
            );
            DEBUG (
                'm', "Set biolog reading for experiment " <<
                experiment->index() << " to " <<  biolog <<
                " growth " << params->is_growth_biolog (biolog)
            );
            continue;
        }
        double lb = tok.nextDouble (error);
        double ub = tok.nextDouble (error);
        if (error)
            reader.error (string ("Bad format in supply_demand file"));
        
        bool is_carbon_source = false;
        bool val = tok.nextBool (error); // Allow blank to mean false
        if (!error) {
            is_carbon_source = val;
        }
        
        auto met = find_metabolite (met_str);
        if (met == nullptr) {
            DEBUG (
                'A', "Lost metabolite " << met_str <<
                " in sd file " << supply_demand_fn
            );
            continue;
        }

        experiment->set_supply_resid (met, lb, ub);
        if (is_carbon_source)
            experiment->add_carbon_source (met);
        
        DEBUG (
            'm', "Set supply/resid range of " << met->name() <<
            " to [" << -lb << "," << -ub << "] = [" <<
            experiment->lb(met) << "," << experiment->ub(met) << "]" <<
            " max1 " << limeMax (lb,ub) <<
            " max2 " << limeMax (-lb,-ub) <<
            " c-src " << is_carbon_source
        );
    }
    // Sort as we read them
    SortPairT<Experiment*> sorter;
    for (auto& exp : experiment_) 
        sorter.add (exp.get(), -exp->biolog_score());
    for (size_t k = 0; k < sorter.size(); k++) {
        biolog_rank_[k] = sorter[k];
        sorter[k]->set_biolog_rank(k);
    }
}

void
Scenario::read_base_flux (string base_flux_fn)
{
    cout << "Reading base flux " << base_flux_fn << endl;
    DEBUG ('d', "Reading base_flux " << base_flux_fn);
    
    has_base_flux_ = true;
    LineReader reader (base_flux_fn);

    // Each line is of the form
    // sd_fn base_flux 
    // where
    //   sd_fn is the name of an experiment file
    //   base_flux is the base_flux to use for rank constraints

    string line;
    while (reader.getLine (line)) {
        LimeTok tok (line);

        bool error = false;
        string sd_fn = tok.nextString (error);
        double flux = tok.nextDouble (error);
        if (error)
            reader.error (string ("Bad format in base_flux file"));
        
        auto exp = find_experiment (basename (sd_fn));
        if (exp == nullptr) {
            DEBUG (
                'A', "Lost experiment " << basename (sd_fn) << 
                " for sd file " << sd_fn <<
                " in base flux file " << base_flux_fn
            );
            continue;
        }
        exp->set_base_flux (flux);
        DEBUG (
            'D', "Set base flux for " << exp->name() <<
            " to " << flux
        );
    }
}

void
Scenario::calc_target_flux (vector<double>& target, double rank0_flux)
{
    DEBUG ('E', "  Calc target flux for rank 0 flux " << rank0_flux);
    size_t biolog0 = biolog_rank0();
    double rank0_biolog = experiment(biolog0)->biolog_score();
    
    target.reserve(num_experiments());
    for (size_t exp = 0 ; exp < num_experiments(); exp++) {
        double exp_biolog = experiment(exp)->biolog_score();
        double flux = rank0_flux * exp_biolog / rank0_biolog;
        DEBUG ('E', "    Exp " << exp << " target flux " << flux);
        if (
            false && //??
            experiment(exp)->has_base_flux() &&
            flux > experiment(exp)->base_flux()
        ) {
            flux = experiment(exp)->base_flux();
            DEBUG ('E', "    Clamp at base flux " << flux);
        }
        target[exp] = flux;
    }
}

int
Scenario::depth ()
{
    vector<int> react_depth;
    vector<bool> avail;
    vector<int> enabled_by;
    int residual_depth;
    int never;

    // Use experiment 0 as indicator
    int depth =
        calc_reachability (
            0, react_depth, avail, enabled_by, residual_depth, never
        );
    DEBUG ('P', "Depth is " << depth);
    
    return depth;
}

int
Scenario::sol_depth (Solution* sol)
{
    vector<int> react_depth;
    vector<bool> avail;
    vector<int> enabled_by;
    int residual_depth;
    int never;

    int depth =
        calc_reachability (
            0, react_depth, avail, enabled_by, residual_depth, never, sol
        );
    DEBUG ('P', "Sol Depth is " << depth);
    
    return depth;
}

/**
   Counts the number of reactions, metabolites and metabolites with
   residuals that *cannot* be reached
*/
void
Scenario::count_reachability (
    size_t exp_idx, int& reactions, int& mets, int& residuals
)
{
    vector<int> react_depth;
    vector<bool> avail;
    vector<int> enabled_by;
    int residual_depth;
    int never;
    auto exp = experiment(exp_idx);

    int depth =
        calc_reachability (
            exp_idx, react_depth, avail, enabled_by, residual_depth, never
        );
    
    reactions = 0;
    for (size_t r = 0; r < num_reactions(); r++) 
        if (react_depth[r] >= never) // Not enabled
            reactions++;

    mets = 0;
    residuals = 0;
    for (size_t k = 0; k < num_metabolites(); k++) 
        if (!avail[k]) {
            mets++;
            if (exp->is_residual(k))
                residuals++;
        }
    
    DEBUG ('P', "Reachable count is " << reactions);
}

int
Scenario::calc_reachability (
    size_t exp_idx, 
    vector<int>& react_depth, vector<bool>& avail, vector<int>& enabled_by,
    int& residual_depth, int& never, Solution* sol // can be null
)
{
    auto exp = experiment(exp_idx);
    // Never is a const that means not yet available
    never = num_reactions()+1;
    residual_depth = never;
    // react_depth is the depth that the reaction became active
    // i.e. all inputs satisfied
    // -1 means reaction not used, so don't care
    // never means hasn't been enabled yet
    react_depth.resize (num_reactions());
    fill (react_depth.begin(), react_depth.end(), never);
        
    // avail is whether metabolite available
    avail.resize(num_metabolites());
    fill (avail.begin(), avail.end(), false);
    
    // amount is only used if solution is provided, and shows
    // the total amount of a metabolite available through enabled
    // reactions
    vector<double> amount (num_metabolites(), 0.0f);
    
    // enabled_by is the reaction that enabled a metabolite
    enabled_by.resize(num_metabolites());
    fill (enabled_by.begin(), enabled_by.end(), -2);
    
    // Make all those with a supply available
    for (size_t k = 0; k < num_metabolites(); k++)
        if (exp->is_supply(k)) {
            DEBUG (
                'P', "Metabolite " << k << " available at depth 0"
            );
            avail[k] = true;
            if (sol != nullptr) 
                amount[k] = exp->supply(k);
            enabled_by[k] = -1;
        }
    int curr_depth = 0;
    bool finished = false;
    bool found_one = true;
    while (!finished && found_one && curr_depth <= num_reactions()) {
        DEBUG ('P', "Depth " << curr_depth);
        finished = true;
        found_one = false;
        curr_depth++;
        for (size_t r = 0; r < num_reactions(); r++) {
            auto react = reaction(r);
            if (react->is_dummy() || react_depth[r] < never) // already active
                continue;
            if (sol != nullptr && limeIsZero (sol->flux(r))) 
                continue; // Not in sol
            
            bool inputs_satisfied = true;
            for (auto k : react->mets()) {
                if (react->met_coeff (k) < 0.0f) {
                    // the met is an input
                    if (!avail[k]) { // not available yet
                        DEBUG ('P', "  Reaction " << r << " needs met " << k);
                        inputs_satisfied = false;
                        break;
                    }
                    else if (
                        sol != nullptr &&
                        amount[k] < sol->flux(r) * -react->met_coeff (k)
                    ) {
                        // Avail, but check enough in amount
                        DEBUG (
                            'P', "    Avail but not enough " <<
                            amount[k] << " < " <<
                            (sol->flux(r) * -react->met_coeff (k))
                        );
                        inputs_satisfied = false;
                        break;
                    }
                }
            }
            if (inputs_satisfied) {
                DEBUG (
                    'P', "  Reaction " << r <<
                    " is available at depth " << curr_depth
                );
                // This reaction has become active
                react_depth[r] = curr_depth;
                found_one = true;
            }
            else {
                DEBUG (
                    'P', "    Reaction " << r <<
                    " is still not available at depth " << curr_depth
                );
                finished = false;
            }
        }
        if (residual_depth == never) {
            // See if all residuals satisfied
            DEBUG ('P', "Check residuals");
            bool inputs_satisfied = true;
            for (size_t k = 0; k < num_metabolites(); k++) {
                if (exp->is_residual(k)) {
                    bool enough = avail[k];
                    if (enough && sol != nullptr)
                        enough =
                            (amount[k] + LIME_EPSILON) >=
                            exp->residual(k);
                    if (!enough) {
                        DEBUG ('P', "  Met " << k << " not enough");
                        inputs_satisfied = false;
                        break;
                    }
                }
            }
            if (inputs_satisfied) {
                DEBUG ('P', "All residuals satisfied at depth " << curr_depth);
                residual_depth = curr_depth;
            }
        }
        for (size_t r = 0; r < num_reactions(); r++) {
            auto react = reaction(r);
            if (react_depth[r] == curr_depth) {
                // Make all products available
                for (auto k : react->mets()) {
                    if (
                        react->met_coeff (k) > 0.0f  // the met is an output
                    ) {
                        if (!avail[k]) {
                            avail[k] = true;
                            if (sol != nullptr)
                                amount[k] +=
                                    sol->flux(k) * react->met_coeff (k);
                            enabled_by[k] = r;
                            DEBUG (
                                'P', "    Metabolite " << k <<
                                " enabled by react " << r
                            );
                        }
                    }
                }
            }
        }
        DEBUG (
            'P', "  End of loop - finished " << finished <<
            " found one " << found_one <<
            " curr_depth " << curr_depth
        );
    }
    if (!finished && !found_one) {
        DEBUG ('P', "No more reactions enabled at depth " << curr_depth);
        // Set to 'infinite' depth
        curr_depth = num_reactions() + 1;
    }
    DEBUG ('P', "Reachability depth is " << curr_depth);
    return curr_depth;
}

size_t
Scenario::make_dijkstra (
    lime::Dijkstra<double>& graph, const Params* params
)
{
    //vector<bool> avail(num_metabolites(), false);
    // Cycles mess us up, so just pretend everything is available
    vector<bool> avail(num_metabolites(), true);
    vector<const Reaction*> supplied_by (num_metabolites(), nullptr);

    return make_dijkstra (graph, params, avail, supplied_by);
}

size_t
Scenario::make_dijkstra (
    lime::Dijkstra<double>& graph, const Params* params,
    vector<bool>& avail, vector<const Reaction*>& supplied_by
)
{
    DEBUG ('J', "Making dijkstra graph");
    size_t dummy_source = num_metabolites();
    graph.setSize(dummy_source + 1);
    
    for (auto& react : reaction_)
        react->set_selected(false);
    
    // Ensure all supply mets are avail
    for (size_t k = 0; k < num_metabolites(); k++) {
        if (experiment_[0]->is_supply(k)) {
            avail[k] = true;
            DEBUG (
                'J', "      Edge Source -> " << *metabolite(k) << " (supply)"
            );
            graph.addEdge (dummy_source, k, 1, num_reactions());
        }
    }
    // Link all sources to dummy source
    for (auto& react : reaction_) {
        if (react->is_active() && react->is_export()) {
            for (auto k : react->mets()) {
                assert (react->met_coeff(k) > 0);
                avail[k] = true;
                supplied_by[k] = react.get();
                graph.addEdge (dummy_source, k, 1, react->index());
                DEBUG (
                    'J', "      Edge Source -> " << *metabolite(k) <<
                    " (exported by " << react->name() << ")"
                );
            }
        }
    }

    bool change_made = true;
    int num_selected = 0;
    while (change_made) {
        change_made = false;

        // See if any reactions are triggered
        for (size_t r = 0; r < num_reactions(); r++) {
            auto react = reaction(r);
            // Skip if already selected, or ignorable
            if (
                react->is_selected() ||
                !react->is_active() ||
                react->is_dummy() ||
                react->is_biomass() ||
                params->exceeds_max_cost (react)
            )
                continue;
            
            DEBUG ('K', "    Consider " << *react);
            bool triggered = true;
            for (auto met : react->in_mets()) {
                if (!avail[met->index()]) {
                    DEBUG ('K', "        Met " << *met << " not avail");
                    triggered = false;
                }
            }
            if (!triggered)
                continue;
            DEBUG ('J', "  Reaction " << *react << " is triggered");
            change_made = true;
            react->set_selected(true);
            num_selected++;

            for (auto met_i : react->in_mets()) {
                double coeff_i = react->met_coeff (met_i);
                size_t i = met_i->index();
                for (auto met_j : react->out_mets()) {
                    size_t j = met_j->index();
                    double coeff_j = react->met_coeff (met_j);
                    graph.addEdge (i, j, react->obj_coeff(), r);
                    if (!avail[j]) {
                        DEBUG ('J', "    Met " << *met_j << " is available");
                    }
                    avail[j] = true;
                    supplied_by[j] = react;
                    DEBUG ('J', "      Edge " << *met_i << " -> " << *met_j);
                }
            }
        }
    }
    
    cout << "Selected " << num_selected << " / " << num_reactions() <<
        " reactions" << endl;
    cout << "Graph has " << graph.num_edges() << " edges" << endl;

    size_t num_not_supplied = 0;
    for (size_t k = 0; k < num_metabolites(); k++) {
        auto met = metabolite(k);
        if (!experiment_[0]->is_supply(k) && supplied_by[k] == nullptr) {
            num_not_supplied++;
        }
    }
    cout << "Num mets not supplied " << num_not_supplied << endl;
    
    return dummy_source;
}

void
Scenario::draw_reachability (Dig* dig, Solution* sol)
{
    vector<int> react_depth;
    vector<bool> avail;
    vector<int> enabled_by;
    int residual_depth;
    int never;

    // Use experiment 0
    auto exp = experiment(0);

    int depth =
        calc_reachability (
            0, react_depth, avail, enabled_by, residual_depth, never, sol
        );
    
    SortPair srt;
    for (int r = 0; r < (int)num_reactions(); r++) {
        if (!reaction(r)->is_dummy()) {
            double key = (double)react_depth[r];
            srt.add (r, key);
        }
    }
    srt.doSort();
    vector<int> revsrt(num_reactions(), 0);
    for (size_t k = 0; k < srt.size(); k++)
        revsrt[srt[k]] = k;

    dig->title ("Reachability");
    dig->xTic (0.0f, "Supply");
    for (size_t i = 0; i < srt.size(); i++) {
        int r = srt[i];
        double x = i + 1.0f;
        dig->xTic (
            x, reaction(r)->name() + "(" + to_string(react_depth[r]) + ")"
        );
    }
    double residual_x = srt.size() + 1.0f;
    dig->xTic (
        residual_x,
        string("Residual(") + to_string(residual_depth) + ")"
    );

    for (size_t j = num_metabolites(); j > 0; j--) {
        double y = num_metabolites() - (j-1);
        dig->yTic (y, metabolite_[j-1]->name());
    }

    int diam = 10;
        
    DEBUG ('P', "Draw reachability");
    
    // Dots for supply and residual
    list<Line> green_lines;
    for (size_t j = 0; j < num_metabolites(); j++) {
        auto met = metabolite(j);
        double y = num_metabolites() - j;
        int colour = 0;
        if (exp->is_supply(met)) {
            DEBUG ('P', "Supply " << j);
            double x = 0.0f;
            dig->style (Dig::StyleColours::BLUE, 0);
            dig->circle (x, y, diam, true);
        }
        if (exp->is_residual(met)) {
            DEBUG ('P', "Residual " << j);
            double x = residual_x;
            dig->style (Dig::StyleColours::RED, 0);
            dig->circle (x, y, diam, true);

            // Enabled-by line
            double enabled_x = 0;
            if (!avail[j]) {
                // Not enabled
                DEBUG ('P', "      Not available");
                colour = Dig::StyleColours::DARK_GREEN;
                int stroke = Dig::StyleStroke::DASHED;
                enabled_x = x - 1.0f;
                dig->style (colour, 0, stroke);
                dig->draw (x, y, enabled_x, y);
            }
            if (enabled_by[j] == -1) {
                enabled_x = 0.0f;
                green_lines.push_back (Line (x, y, enabled_x, y));
            }
            else if (enabled_by[j] >= 0) {
                DEBUG (
                    'P', "      Enabled by " << enabled_by[j] <<
                    " srt " << revsrt[enabled_by[j]]
                );
                enabled_x = revsrt[enabled_by[j]] + 1.0f;
                green_lines.push_back (Line (x, y, enabled_x, y));
            }
        }
    }
    int prev_depth = -1;
    for (size_t i = 0; i < srt.size(); i++) {
        size_t r = srt[i];
        auto react = reaction(r);
        
        DEBUG (
            'P', "  Reaction " << r << " at pos " << i <<
            " depth " << react_depth[r]
        );
        
        double x = i + 1.0f;

        if (prev_depth != react_depth[r]) {
            prev_depth = react_depth[r];
            dig->style (Dig::StyleColours::CORAL, 0);
            dig->label ("Level " + to_string (prev_depth));
            dig->draw (x-0.5f, 0.0f, x-0.5f, num_metabolites() + 1.0f);
        }

        // One row for each metabolite
        for (size_t j : react->mets()) {
            auto met = metabolite(j);
            double y = num_metabolites() - j;
            DEBUG (
                'P', "    Met " << j << " at y " << y 
            );

            int colour = 0;
            if (react->met_coeff(met) < 0.0f) {
                colour = Dig::StyleColours::RED;
            }
            else {
                colour = Dig::StyleColours::BLUE;
            }
            dig->style (colour, 0);
            dig->circle (x, y, diam, true);
            
            double enabled_x = 0;
            if (react->met_coeff(met) < 0.0f) {
                // Only draw lines for consumers
                if (!avail[j]) {
                    // Not enabled
                    DEBUG ('P', "      Not available");
                    colour = Dig::StyleColours::DARK_GREEN;
                    int stroke = Dig::StyleStroke::DASHED;
                    enabled_x = x - 1.0f;
                    dig->style (colour, 0, stroke);
                    dig->draw (x, y, enabled_x, y);
                }
                else {
                    if (enabled_by[j] == -1) {
                        DEBUG ('P', "      Enabled by supply");
                        enabled_x = 0.0f;
                        green_lines.push_back (Line (x, y, enabled_x, y));
                    }
                    else if (enabled_by[j] >= 0) {
                        DEBUG (
                            'P', "      Enabled by " << enabled_by[j] <<
                            " srt " << revsrt[enabled_by[j]]
                        );
                        enabled_x = revsrt[enabled_by[j]] + 1.0f;
                        green_lines.push_back (Line (x, y, enabled_x, y));
                    }
                }
            }
        }
    }
    dig->wait();
    int colour = Dig::StyleColours::DARK_GREEN;
    int stroke = Dig::StyleStroke::THIN;
    dig->style (colour, 0, stroke);
    for (auto& line : green_lines)
        dig->draw (&line);
    
    DEBUG ('P', "Done draw reachability");
}
