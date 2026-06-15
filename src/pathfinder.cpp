/**
 * @file pathfinder.cpp
 * @brief Implementation of PathFinder class for metabolic pathway analysis
 *
 * This file implements the PathFinder class which performs shortest path finding,
 * loop detection, and connectivity analysis in metabolic networks. It supports
 * gap-filling analysis by identifying biomass precursor accessibility and
 * detecting metabolic cycles with flux constraints.
 */

#include <sstream>
#include <algorithm>
#include <list>
#include <iomanip>
#include <math.h>

#include "lime/debug.h"
#include "lime/dijkstra.h"

#include "mosh/pathfinder.h"
#include "mosh/scenario.h"

using namespace std;
using namespace lime;
using namespace mosh;

#define write(X) {cout << X; if (out_ != nullptr) *(out_) << X;}

/**
 * @brief Find shortest path from source to a target metabolite
 *
 * Resets graph costs and computes the shortest path from the source node
 * to the specified metabolite using Dijkstra's algorithm.
 *
 * @param met Target metabolite to find path to
 * @return true if a valid path exists, false otherwise
 */
bool
PathFinder::find_sp_to (const Metabolite* met)
{
    graph_.reset_costs();
    return graph_.findPath (source_idx_, met->index());
}

/**
 * @brief Check connectivity of biomass reaction precursors
 *
 * Analyzes whether metabolites required by the biomass reaction can be
 * reached from source metabolites. Traces back paths through reactions,
 * records maximum reaction costs, and builds a graph representation for
 * visualization. Outputs connectivity statistics and path information.
 *
 * @param[in,out] count Total number of metabolites checked (incremented)
 * @param[in,out] count_ok Number of metabolites with valid paths (incremented)
 */
void
PathFinder::check_biomass (int& count, int& count_ok)
{
    const Reaction* biomass = scenario_->biomass_react();
    int biomass_mets = 0;
    int biomass_mets_ok = 0;
    double max_react_cost = 0.0f;
    for (auto met : biomass->in_mets()) {
        count++;
        biomass_mets++;
        if (find_sp_to (met)) {
            write ("Path to " << *met << endl);
            count_ok++;
            biomass_mets_ok++;
            write ("  ");
            auto upto = met;
            const Reaction* prev_react = nullptr;
            DEBUG ('P', "Path to " << *met);
            while (upto != nullptr) {
                dot_mets_.insert (upto);
                DEBUG ('P', "  Upto is " << *upto << " idx " << upto->index());
                write (upto->name());
                if (prev_react != nullptr)
                    add_dot_edge (
                        upto->name(), prev_react->name(), prev_react
                    );
                auto react = react_for (upto);
                DEBUG ('P', "    React idx " << react_idx_for (upto));
                write (" [ ");
                if (react_idx_for (upto) >= scenario_->num_reactions()) {
                    // It's a source met
                    write ("Source");
                    add_dot_edge ("Source", upto->name());
                    prev_react = nullptr;
                }
                else {
                    DEBUG ('P', "    React " << react->name());
                    write (
                        react->name() << " " <<
                        limeFormat ("%.3lf", react->obj_coeff())
                    );
                    dot_reacts_.insert (react);
                    if (react-> is_export()) {
                        write (" (export)");
                        add_dot_edge ("Export", react->name());
                    }
                    add_dot_edge (react->name(), upto->name(), react);
                    if (react->obj_coeff() > max_react_cost) 
                        max_react_cost = react->obj_coeff();
                    prev_react = react;
                }
                write (" ]" << endl << "  ");
                upto = parent (upto);
            }
            write (endl);
        }
        else {
            write ("NO Path to " << *met << endl);
            dot_mets_.insert (met);
            add_dot_edge ("NoPath", met->name());
        }
    }
    write (
        "Found paths to " << biomass_mets_ok << " out of " <<
        biomass_mets << " reactants in biomass reaction - " <<
        limeFormat (
            "%5.1lf", (double) 100.0f * biomass_mets_ok / biomass_mets
        ) << "%" << endl
    );
    write (
        "Max react cost " << limeFormat ("%.3lf", max_react_cost) << endl
    );
}

/**
 * @brief Perform comprehensive biomass reaction connectivity check
 *
 * Similar to check_biomass but with extended analysis. Examines all biomass
 * precursor metabolites, traces paths back to sources, and builds complete
 * graph structures including export reactions. Used for detailed gap analysis.
 *
 * @param[in,out] count Total number of metabolites checked (incremented)
 * @param[in,out] count_ok Number of metabolites with valid paths (incremented)
 */
void
PathFinder::full_check_biomass (int& count, int& count_ok)
{
    const Reaction* biomass = scenario_->biomass_react();
    list<const Metabolite*> target;
    set<const Reaction*> exclude;
    
    for (auto met : biomass->in_mets()) {
        if (find_sp_to (met)) {
            auto upto = met;
            const Reaction* prev_react = nullptr;
            DEBUG ('P', "Path to " << *met);
            while (upto != nullptr) {
                DEBUG ('P', "  Upto is " << *upto << " idx " << upto->index());
                auto react = react_for (upto);
                DEBUG ('P', "    React idx " << react_idx_for (upto));
                if (react_idx_for (upto) >= scenario_->num_reactions()) {
                    // It's a source met
                    prev_react = nullptr;
                }
                else {
                    DEBUG ('P', "    React " << react->name());
                    write (
                        react->name() << " " <<
                        limeFormat ("%.3lf", react->obj_coeff())
                    );
                    dot_reacts_.insert (react);
                    if (react-> is_export()) {
                        write (" (export)");
                        add_dot_edge ("Export", react->name());
                    }
                    add_dot_edge (react->name(), upto->name(), react);
                    prev_react = react;
                }
                write (" ]" << endl << "  ");
                upto = parent (upto);
            }
            write (endl);
        }
        else {
            write ("NO Path to " << *met << endl);
            dot_mets_.insert (met);
            add_dot_edge ("NoPath", met->name());
        }
    }
}

/**
 * @brief Build metabolic graph reachable from a carbon source
 *
 * Performs forward propagation from a specified carbon source metabolite
 * to determine which reactions can be activated and which metabolites can
 * be produced. Uses iterative flood-fill algorithm limited by carbon_depth
 * to simulate metabolic network expansion from the carbon source. Marks
 * reactions and metabolites for graph visualization.
 *
 * @param carbon_source Starting carbon source metabolite
 * @param carbon_depth Maximum number of reaction layers to propagate
 * @param dot Output stream for DOT graph format (currently unused in implementation)
 */
void
PathFinder::carbon_source_graph (
    const Metabolite* carbon_source, int carbon_depth, ostream& dot
)
{
    // Marked mets are directly supplied by carbon source, or can track back 
    vector<bool> marked (scenario_->num_metabolites(), false);
    marked[carbon_source->index()] = true;
    
    DEBUG ('f', "Check supply");
    fill (supplied_by_.begin(), supplied_by_.end(), nullptr);
    dot_reacts_.insert (source_.get());
    supplied_by_[carbon_source->index()] = source_.get();
    dot_mets_.insert (carbon_source);
    bool found_carbon = false;
    for (auto& met : scenario_->metabolites()) {
        if (scenario_->experiment(0)->is_supply(met->index())) {
            DEBUG ('f', "  " << *met << " is source (supply)");
            met->set_is_source(true);
            if (met.get() == carbon_source)
                found_carbon = true;
            supplied_by_[met->index()] = source_.get();
            dot_mets_.insert (met.get());
            add_dot_edge (source_->name(), met->name());
        }
    }
    if (!found_carbon) {
        add_dot_edge (source_->name(), carbon_source->name());
    }
    
    // List of not-yet-activated reactions
    list<Reaction*> waiting;
    for (auto& react : scenario_->reactions()) {
        waiting.push_back(react.get());
        react->set_active(false);
        react->set_selected(false);
    }
        
    bool something_changed = true;
    int iter = 0;
    while (something_changed && iter++ < carbon_depth) {
        DEBUG ('f', "Next iter");
        something_changed = false;
        auto iter = waiting.begin();
        while (iter != waiting.end()) {
            auto react = *iter;
            DEBUG ('F', "  Check "<< react->name());
            bool ok = true;
            bool from_marked = false;
            for (auto met : react->in_mets()) {
                if (!is_supplied(met)) {
                    ok = false;
                    DEBUG ('F', "    Waiting for "<< met->name());
                }
                if (marked[met->index()])
                    from_marked = true;
            }
            if (ok) {
                // Reaction is enabled
                DEBUG ('f', "  Enabled " << react->name());
                something_changed = true;
                react->set_active(true);
                
                if (from_marked) {
                    dot_reacts_.insert (react);
                    react->set_selected(true);
                    for (auto met : react->in_mets()) {
                        dot_mets_.insert (met);
                        add_dot_edge (met->name(), react->name());
                    }
                }
                iter = waiting.erase (iter);
                
                for (auto met : react->out_mets()) {
                    DEBUG ('f', "      Activate " << met->name());
                    supplied_by_[met->index()] = react;
                    if (from_marked) {
                        dot_mets_.insert (met);
                        add_dot_edge (react->name(), met->name());
                        marked[met->index()] = true;
                    }
                }
            }
            else {
                ++iter;
            }
        }
    }
    for (auto& react : scenario_->reactions()) {
        DEBUG (
            'f', "Reaction " << react->name() <<
            " active? " << react->is_active()
        );
    }
    for (auto& met : scenario_->metabolites()) {
        DEBUG (
            'f', "Met  " << met->name() <<
            " supplied? " << is_supplied(met.get())
        );
    }
    
}

/**
 * @brief Calculate maximum flux in a metabolic loop
 *
 * Traverses a detected loop starting from the given metabolite index and
 * determines the maximum flux value among all reactions participating in
 * the loop. The maximum flux represents the loop's capacity.
 *
 * @param start_idx Index of metabolite where loop detection starts
 * @return Maximum flux value found in the loop
 */
double
PathFinder::calc_flux (size_t start_idx)
{
    DEBUG ('f', "Calc flux for " << start_idx);
    double flux = 0.0f;
    size_t upto_idx = start_idx;
    const Reaction* react = nullptr;
    int count = 0;
    do {
        DEBUG ('f', "  upto " << upto_idx);
        auto react_idx = graph_.parent_edge (upto_idx);
        DEBUG ('f', "  react idx  " << react_idx);
        assert (react_idx < scenario_->num_reactions());
        auto react = scenario_->reaction (react_idx);
        if (react->known_flux() > flux) 
            flux = react->known_flux();
        upto_idx = graph_.parent (upto_idx);
    }
    while (upto_idx != start_idx && count++ < scenario_->num_reactions());
    return flux;
}

/**
 * @brief Detect and analyze a metabolic loop from a given metabolite
 *
 * Searches for a cycle in the metabolic graph starting from the specified
 * metabolite. If found, calculates the loop's flux and filters out low-flux
 * loops below the minimum threshold. Records all metabolites and reactions
 * participating in the loop for visualization and counts unique metabolites.
 *
 * @param met Starting metabolite for loop detection
 * @param min_flux Minimum flux threshold for valid loops
 * @param[in,out] num_mets Counter for unique metabolites found in loops (incremented)
 * @return true if a valid high-flux loop is found, false otherwise
 */
bool
PathFinder::find_loop (const Metabolite* met, double min_flux, int& num_mets)
{
    DEBUG ('l', "  Look for loop from " << *met);
    write ("Check for loop from " << *met << endl);
    if (graph_.findLoop (met->index())) {
        write ("  Found loop" << endl);
        DEBUG ('l', "    Found loop involving " << *met);
        auto upto_idx = met->index();
        double flux = calc_flux (upto_idx);
        DEBUG ('l', "      Flux " << flux);
        if (flux < min_flux) {
            write ("  Low-flux loop" << endl);
            DEBUG ('l', "    Low-flux loop for " << *met << ": " << flux);
            return false;
        }
        const Reaction* react = nullptr;
        int count = 0;
        double max_flux = 0.0f;
        do {
            auto upto = scenario_->metabolite(upto_idx);
            DEBUG ('l', "        Met " << upto_idx << " " << *upto);
            if (react != nullptr)
                add_dot_edge (upto->name(), react->name(), react);
            auto react_idx = graph_.parent_edge (upto_idx);
            assert (react_idx < scenario_->num_reactions());
            react = scenario_->reaction (react_idx);
            DEBUG ('l', "        React " << react_idx << " " << *react);
            write ("    " << *upto << " via " << *react);
            count++;
            if (react->known_flux() > max_flux)
                max_flux = react->known_flux();
            if (react->known_flux() > 0.0f) {
                write (" (flux " <<
                    fixed << setprecision(2) << react->known_flux() <<
                       ")");
            }
            write (endl);
            DEBUG (
                'l', "      Loop met " << *upto <<
                " via " << *react
            );
            dot_mets_.insert (upto);
            dot_reacts_.insert (react);
            add_dot_edge (react->name(), upto->name(), react);
            if (!in_loop_[upto_idx]) {
                num_mets++;
                in_loop_[upto_idx] = true;
            }
            upto_idx = graph_.parent (upto_idx);
        }
        while (upto_idx != met->index());
                
        if (react != nullptr) {
            DEBUG (
                'l', "      Close loop to " << *met <<
                " via " << *react
            );
            write ("    " << *met << " via " << *react << endl);
            count++;
            if (react->known_flux() > max_flux)
                max_flux = react->known_flux();
            write (
                "      Loop len " << count << " max flux " << max_flux << endl
            );
            
            add_dot_edge (met->name(), react->name(), react);
        }
        return true;
    }
    write ("  No loop" << endl);
    DEBUG ('l', "    No loop for " << *met);
    return false;
}

/**
 * @brief Find all metabolic loops in the network above flux threshold
 *
 * Systematically searches for loops starting from each non-source metabolite
 * that is not already part of a detected loop. Counts total loops found and
 * tracks total number of metabolites participating in any loop. Used for
 * identifying futile cycles and metabolic modules.
 *
 * @param min_flux Minimum flux threshold to consider a loop significant
 * @param[in,out] num_loops Counter for total loops detected (incremented)
 * @param[in,out] num_mets Counter for total unique metabolites in loops (incremented)
 */
void
PathFinder::find_loops (double min_flux, int& num_loops, int& num_mets)
{
    DEBUG ('l', "Find loops");
#ifndef NDEBUG
    if (Debug::doDebug('L')) {
        graph_.show (Debug::debugFile());
    }
#endif
    for (auto& metptr : scenario_->metabolites()) {
        graph_.reset_costs();
        auto met = metptr.get();
        if (!met->is_source() && !in_loop_[met->index()]) {
            if (find_loop (met, min_flux, num_mets)) {
                num_loops++;
            }
        }
    }
}

/**
 * @brief Perform flood-fill to determine metabolite supply connectivity
 *
 * Identifies all metabolites that can be synthesized from available source
 * metabolites by iteratively activating reactions whose substrates are all
 * available. Continues until no new reactions can be activated. This forward
 * propagation analysis determines network connectivity and identifies gaps.
 */
void
PathFinder::flood()
{
    DEBUG ('f', "Check supply");
    fill (supplied_by_.begin(), supplied_by_.end(), nullptr);
    for (auto& met : scenario_->metabolites()) {
        if (scenario_->experiment(0)->is_supply(met->index())) {
            DEBUG ('f', "  " << *met << " is source (supply)");
            met->set_is_source(true);
            supplied_by_[met->index()] = source_.get();
        }
    }
    
    // List of not-yet-activated reactions
    list<const Reaction*> waiting;
    for (auto& react : scenario_->reactions()) {
        waiting.push_back(react.get());
        react->set_active(false);
    }
        
    bool something_changed = true;
    while (something_changed) {
        DEBUG ('f', "Next iter");
        something_changed = false;
        auto iter = waiting.begin();
        while (iter != waiting.end()) {
            auto react = *iter;
            DEBUG ('f', "  Check "<< react->name());
            bool ok = true;
            for (auto met : react->in_mets()) {
                if (!is_supplied(met)) {
                    ok = false;
                    DEBUG ('f', "    Waiting for "<< met->name());
                }
            }
            if (ok) {
                // Reaction is enabled
                DEBUG ('f', "    Enabled");
                something_changed = true;
                ((Reaction*)react)->set_active(true);
                iter = waiting.erase (iter);
                for (auto met : react->out_mets()) {
                    DEBUG ('f', "      Activate " << met->name());
                    supplied_by_[met->index()] = react;
                }
            }
            else {
                ++iter;
            }
        }
    }
    for (auto& react : scenario_->reactions()) {
        DEBUG (
            'f', "Reaction " << react->name() <<
            " active? " << react->is_active()
        );
    }
    for (auto& met : scenario_->metabolites()) {
        DEBUG (
            'f', "Met  " << met->name() <<
            " supplied? " << is_supplied(met.get())
        );
    }
}

/**
 * @brief Generate DOT format graph for visualization
 *
 * Outputs metabolic network graph in GraphViz DOT format including all
 * metabolites and reactions collected during analysis. Supports optional
 * explosion modes that expand edges to show all reaction substrates/products
 * or all reactions using specific metabolites.
 *
 * @param dot Output stream to write DOT format graph
 * @param explode_react If true, add explicit edges for all reaction inputs/outputs
 * @param explode_met If true, add all reactions using each metabolite in dot_mets_
 */
void
PathFinder::write_dot (ostream& dot, bool explode_react, bool explode_met)
{
    if (explode_react) {
        for (auto react : dot_reacts_) {
            for (auto met : react->in_mets()) {
                add_dot_edge (met->name(), react->name(), react);
            }
            for (auto met : react->out_mets()) {
                add_dot_edge (react->name(), met->name(), react);
            }
        }
    }
    if (explode_met) {
        for (auto met : dot_mets_) {
            for (auto& react_ptr : scenario_->reactions()) {
                auto react = react_ptr.get();
                if (!react->uses (met))
                    continue;
                dot_reacts_.insert (react);
                if (react->met_coeff(met) < 0) {
                    // It's an input met
                    add_dot_edge (met->name(), react->name(), react);
                }
                else {
                    add_dot_edge (react->name(), met->name(), react);
                }
            }
        }
    }
    
    dot << "digraph G {" << endl;
    for (auto met : dot_mets_) {
        dot << "  \"" << met->name() <<
            "\" [color=blue tooltip=\"" << met->full_name() << "\"]" << endl;
    }
    for (auto react : dot_reacts_) {
        string colour = "red";
        if (react->is_selected())
            colour="orange";
        dot << "  \"" << react->name() <<
            "\" [shape=box color=" << colour <<
            " tooltip=\"" << react->formula() << "\"]" <<
            endl;
    }
    for (auto edge : dot_edges_) {
        dot << "  " << edge << endl;
    }
    dot << "}" << endl;
}


