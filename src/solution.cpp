/**
 * @file solution.cpp
 * @brief Implementation of Solution class for metabolic network flux analysis
 *
 * This file contains the implementation of the Solution class which represents
 * a solution to a flux balance analysis problem. It provides methods for:
 * - Computing objective function values (absolute and relative)
 * - Analyzing biomass and dummy reaction fluxes
 * - Calculating metabolite balances and production rates
 * - Detecting runaway reactions based on carbon source constraints
 * - Generating various output formats (text, DOT graphs, visualization)
 * - Path analysis using Dijkstra's algorithm for metabolic pathways
 */

#include <sstream>
#include <algorithm>
#include <list>
#include <math.h>

#include "lime/debug.h"
#include "lime/linereader.h"
#include "lime/limetok.h"
#include "lime/numutil.h"
#include "lime/sortpair.h"
#include "lime/constants.h"
#include "lime/line.h"
#include "lime/dijkstra.h"

#include "mosh/solution.h"
#include "mosh/metabolite.h"

using namespace std;
using namespace lime;
using namespace mosh;

/**
 * @brief Computes the total objective function value for the solution
 *
 * Calculates either absolute or relative objective value based on parameters,
 * and adds the biomass contribution multiplied by the given factor.
 *
 * @param biomass_mult Multiplier for biomass reaction contribution
 * @return Total objective function value
 */
double
Solution::obj_value (double biomass_mult) const
{
    if (params_->use_abs_obj)
        return abs_obj_value () + biomass_obj(biomass_mult);
    return rel_obj_value() + biomass_obj(biomass_mult);
}

/**
 * @brief Computes the relative objective function value
 *
 * Sums the weighted flux values for all non-biomass reactions, where each
 * reaction's contribution is its flux multiplied by its objective coefficient.
 *
 * @return Relative objective function value
 */
double
Solution::rel_obj_value () const
{
    double sum = 0.0f;
    for (size_t k = 0; k < scenario_->num_reactions(); k++) {
        if (!scenario_->reaction(k)->is_biomass())
            sum += scenario_->reaction(k)->obj_coeff() * flux(k);
    }
    DEBUG ('o', "Rel obj val " << sum);
    return sum;
}

/**
 * @brief Computes the absolute objective function value
 *
 * Sums the objective coefficients for all used non-biomass reactions,
 * regardless of their flux magnitude. This gives a count-like metric
 * weighted by reaction costs.
 *
 * @return Absolute objective function value
 */
double
Solution::abs_obj_value () const
{
    double sum = 0.0f;
    for (size_t k = 0; k < scenario_->num_reactions(); k++) {
        if (uses_react(k) && !scenario_->reaction(k)->is_biomass()) {
            DEBUG (
                'o', "  Reaction " << scenario_->reaction(k)->name() <<
                " flux " << flux(k) <<
                " obj " << scenario_->reaction(k)->obj_coeff()
            );
            sum += scenario_->reaction(k)->obj_coeff();
        }
    }
    DEBUG ('o', "Abs obj val " << sum);
    return sum;
}

/**
 * @brief Finds the maximum objective coefficient among used reactions
 *
 * @return Maximum objective coefficient value, or 0.0 if no reactions used
 */
double
Solution::max_cost () const
{
    double the_max = 0.0f;
    for (size_t k = 0; k < scenario_->num_reactions(); k++) {
        if (uses_react(k) && scenario_->reaction(k)->obj_coeff() > the_max)
            the_max = scenario_->reaction(k)->obj_coeff();
    }
    return the_max;
}

/**
 * @brief Counts the number of reactions with non-zero flux
 *
 * @return Number of active reactions in the solution
 */
size_t
Solution::num_reactions_used () const
{
    size_t count = 0;
    for (size_t k = 0; k < scenario_->num_reactions(); k++) {
        if (uses_react(k))
            count++;
    }
    return count;
}

/**
 * @brief Computes the sum of all reaction fluxes
 *
 * @return Total flux across all reactions
 */
double
Solution::sum_flux () const
{
    double sum = 0.0f;
    for (size_t k = 0; k < scenario_->num_reactions(); k++) 
        sum += flux(k);
    return sum;
}

/**
 * @brief Gets the flux value of the biomass reaction
 *
 * The biomass reaction represents the growth rate or biomass production
 * in the metabolic network.
 *
 * @return Flux through the biomass reaction
 */
double
Solution::biomass_flux () const
{
    size_t idx = scenario_->biomass_react()->index();
    return flux(idx);
}

/**
 * @brief Computes the biomass contribution to the objective function
 *
 * Calculates biomass flux multiplied by the given multiplier and the
 * biomass reaction's objective coefficient.
 *
 * @param biomass_mult Multiplier for biomass contribution
 * @return Biomass objective function contribution
 */
double
Solution::biomass_obj (double biomass_mult) const
{
    DEBUG (
        'o', "Biomass obj = flux " << biomass_flux () <<
        " * mult " << biomass_mult <<
        " = " << biomass_flux() * biomass_mult
    );
    return biomass_flux() * biomass_mult *
        scenario_->biomass_react()->obj_coeff();
}

/**
 * @brief Computes the sum of fluxes through dummy reactions
 *
 * Dummy reactions are artificial reactions added for gap-filling or
 * completing metabolic pathways.
 *
 * @return Total flux through all dummy reactions
 */
double
Solution::dummy_flux () const
{
    double sum = 0.0f;
    for (size_t k = 0; k < scenario_->num_reactions(); k++)
        if (scenario_->reaction(k)->is_dummy())
            sum += flux(k);
    return sum;
}

/**
 * @brief Computes the cost contribution of dummy reactions
 *
 * Sums the objective coefficients of all used dummy reactions to
 * quantify the penalty for using artificial reactions.
 *
 * @return Total objective cost of dummy reactions
 */
double
Solution::dummy_obj_val () const
{
    double sum = 0.0f;
    for (size_t k = 0; k < scenario_->num_reactions(); k++)
        if (scenario_->reaction(k)->is_dummy() && uses_react(k))
            sum += scenario_->reaction(k)->obj_coeff();
    return sum;
}

/**
 * @brief Counts the number of dummy reactions with non-zero flux
 *
 * @return Number of active dummy reactions
 */
int
Solution::num_dummy() const
{
    int count = 0;
    for (size_t k = 0; k < scenario_->num_reactions(); k++)
        if (scenario_->reaction(k)->is_dummy() && uses_react(k))
            count++;
    return count;
}

/**
 * @brief Writes flux values for all used reactions to output stream
 *
 * Outputs reaction name, flux value, objective coefficient, and participating
 * metabolites with their stoichiometric coefficients. Includes known flux
 * constraints if specified.
 *
 * @param out Output stream to write flux data
 */
void
Solution::write_flux (std::ostream& out)
{
    out << "# <reaction> <flux> <obj-coeff> {<metabolite> <coeff>}..." << endl;
    for (size_t k = 0; k < scenario_->num_reactions(); k++) {
        if (uses_react(k)) {
            auto react = scenario_->reaction(k);
            out << react->name() << " " << flux(k) <<
                "  " << react->obj_coeff() << " ";
            for (size_t m : react->mets()) {
                out << " " << scenario_->metabolite(m)->name() <<
                    " " << react->met_coeff(m);
            }
            if (react->has_known_flux())
                out << " # Known flux " << react->known_flux();
            out << endl;
        }
    }

    /*
    out << "# Unused reactions with reduced_cost" << endl;
    out << "# <reaction> <col-val> <obj-coeff> {<metabolite> <coeff>}..." << endl;
    SortPair srt;
    for (size_t i = 0; i < scenario_->num_reactions(); i++) {
        // If used, sort to the front, in order
        if (!scenario_->reaction(i)->is_dummy() && limeIsZero(flux(i))) 
            srt.add (i, scenario_->reaction(i)->reduced_cost());
    }
    srt.doSort();
    for (size_t k = 0; k < srt.size(); k++) {
        auto react = scenario_->reaction(srt[k]);
        out << "### " << react->name() << " " << react->reduced_cost() <<
            "  " << react->obj_coeff() << " ";
        for (size_t m : react->mets()) {
            out << " " << scenario_->metabolite(m)->name() <<
                " " << react->met_coeff(m);
        }
        out << endl;
    }
    */
}

/**
 * @brief Checks if the solution has runaway reactions in any experiment
 *
 * A runaway solution occurs when the biomass flux exceeds the available
 * carbon sources, indicating an unrealistic or infeasible solution.
 *
 * @return true if runaway detected in any experiment, false otherwise
 */
bool
Solution::is_runaway() const
{
    for (size_t k = 0; k < scenario_->num_experiments(); k++) {
        if (is_runaway(scenario_->experiment(k)))
            return true;
    }
    return false;
}

/**
 * @brief Checks if the solution has runaway reactions for a specific experiment
 *
 * Compares the total biomass flux against the sum of available carbon sources.
 * A solution is considered runaway if biomass production exceeds carbon availability.
 *
 * @param exp Pointer to the experiment to check
 * @return true if biomass flux exceeds available carbon sources, false otherwise
 */
bool
Solution::is_runaway(const Experiment* exp) const
{
    DEBUG ('r', "      Runaway check for " << exp->name());
    if (exp->carbon_sources().size() == 0) {
        // We don't have any data on carbon sources, so we can't check
        return false;
    }
    double sum_carbon = 0.0f;
    for (auto met : exp->carbon_sources()) {
        DEBUG (
            'r', "        C source " << met->name() <<
            " ub " << exp->ub (met) << " lb " << exp->lb (met)
        );
        sum_carbon += fabs (exp->lb (met));
    }
    DEBUG (
        'r', "      sum c " << sum_carbon << " biomass flux " << biomass_flux()
    );
    return biomass_flux() > sum_carbon;
}


/**
 * @brief Counts the total number of reactions with non-zero flux
 *
 * Similar to num_reactions_used() but includes debug output for each
 * counted reaction.
 *
 * @return Number of active reactions
 */
int
Solution::reaction_count () const
{
    int count = 0;
    for (auto& react_ptr : scenario_->reactions()) {
        const Reaction* react = react_ptr.get();
        if (uses_react(react)) {
            DEBUG (
                'g', " Counting react " << *react << " val " << flux(react)
            ); 
            count++;
        }
    }
    return count;
}

/**
 * @brief Writes the metabolic model including only used reactions
 *
 * Outputs all metabolites involved in active reactions, followed by the
 * reaction definitions with optionally normalized costs.
 *
 * @param out Output stream for model data
 * @param cost_one If true, normalize high costs to a standard value
 */
void
Solution::write_model (std::ostream& out, bool cost_one)
{
    // Work out the mets we need
    set<const Metabolite*> mets;
    for (size_t k = 0; k < scenario_->num_reactions(); k++) {
        auto react = scenario_->reaction(k);
        if (uses_react(k)) {
            for (auto k : react->mets())
                mets.insert (scenario_->metabolite(k));
        }
    }
    for (auto met : mets) {
        met->write_to(out);
    }
    for (size_t k = 0; k < scenario_->num_reactions(); k++) {
        if (uses_react(k)) {
            auto react = scenario_->reaction(k);
            double cost = react->obj_coeff();
            double one = params_->gene_ind_cost + 0.01f;
            if (cost_one && cost > one)
                cost = one;
            react->write_to(out, cost);
        }
    }
}

/**
 * @brief Generates a DOT format graph of metabolite connections
 *
 * Creates a directed graph showing connections between reactions through
 * shared metabolites. Edges represent metabolite flow between reactions.
 *
 * @param out Output stream for DOT graph data
 */
void
Solution::write_met_dot (std::ostream& out)
{
    double max_flux = 0.0f;
    for (size_t k = 0; k < scenario_->num_reactions(); k++) 
        if (flux(k) > max_flux)
            max_flux = flux(k);

    out << "digraph G {" << endl;
    for (size_t i = 0; i < scenario_->num_reactions(); i++) {
        if (limeIsZero (flux(i)))
            continue;
        auto react_i = scenario_->reaction(i);
        for (size_t j = 0; j < scenario_->num_reactions(); j++) {
            if (i == j || limeIsZero (flux(j)))
                continue;
            auto react_j = scenario_->reaction(j);
            
            for (size_t k : react_i->mets()) {
                if (
                    react_i->met_coeff(k) < 0 &&
                    react_j->met_coeff(k) < 0
                ) {
                    // There's a link
                    out << "  " << react_i->name() <<
                        " -> " << react_j->name() <<
                        //scenario_->metabolite(i)->name() <<
                        ";" << endl;
                    break;
                }
            }
        }
    }
    out << "}" << endl;
}

/**
 * @brief Generates a DOT graph showing paths from carbon sources to biomass
 *
 * Traces metabolic pathways from the primary carbon source to biomass
 * precursors using Dijkstra's shortest path algorithm.
 *
 * @param out Output stream for DOT graph data
 */
void
Solution::write_carbon_source_dot (std::ostream& out)
{
    set<const Metabolite*> from_mets;
    set<const Metabolite*> to_mets;

    DEBUG ('t', "Write carbon source dot");
    const Metabolite* carbon_source = nullptr;
    auto exp = scenario_->experiment(0);
    for (size_t k = 0; k < scenario_->num_metabolites(); k++)
        if (exp->is_supply(k) && exp->supply(k) < 100.0) {
            carbon_source = scenario_->metabolite(k);
            break;
        }
    if (carbon_source == nullptr) {
        DEBUG ('t', "No carbon source!");
        return;
    }
    DEBUG ('t', "Carbon source is " << *carbon_source);
    from_mets.insert (carbon_source);

    // Find the biomass reactants
    for (size_t k : scenario_->biomass_react()->mets()) {
        if (scenario_->biomass_react()->met_coeff(k) < 0) {
            to_mets.insert (scenario_->metabolite(k));
        }
    }

    plot_paths (out, from_mets, to_mets);
}

/**
 * @brief Generates a DOT graph showing paths from all supplies to demands
 *
 * Traces metabolic pathways from all supplied metabolites to all residual
 * (demanded) metabolites using shortest path analysis.
 *
 * @param out Output stream for DOT graph data
 */
void
Solution::write_supply_demand_dot (std::ostream& out)
{
    set<const Metabolite*> from_mets;
    set<const Metabolite*> to_mets;

    DEBUG ('t', "Write sd dot");
    const Metabolite* carbon_source = nullptr;
    auto exp = scenario_->experiment(0);
    for (size_t k = 0; k < scenario_->num_metabolites(); k++) {
        if (exp->is_supply(k)) {
            from_mets.insert (scenario_->metabolite(k));
        }
        else if (exp->is_residual(k)) {
            to_mets.insert (scenario_->metabolite(k));
        }
    }

    plot_paths (out, from_mets, to_mets);
}


/**
 * @brief Generates a DOT graph with shortest paths between metabolite sets
 *
 * Uses Dijkstra's algorithm to find paths from each source metabolite to each
 * target metabolite through the reaction network. Visualizes reactions as boxes,
 * metabolites as circles, and includes flux information on edges.
 *
 * @param out Output stream for DOT graph data
 * @param from_mets Set of source metabolites
 * @param to_mets Set of target metabolites
 */
void
Solution::plot_paths (
    std::ostream& out,
    set<const Metabolite*>& from_mets,
    set<const Metabolite*>& to_mets
)
{
    DEBUG (
        't', "Plot paths from " << from_mets.size() << " mets to " <<
        to_mets.size() << " mets"
    );
    DEBUG ('t', "Make Dijkstra graph");
    Dijkstra<int> graph (scenario_->num_metabolites());
    
    for (size_t k = 0; k < scenario_->num_reactions(); k++) {
        if (limeIsZero (flux(k)))
            continue;
        auto react = scenario_->reaction(k);

        for (size_t i : react->mets()) {
            if (react->met_coeff(i) < 0) {
                for (size_t j : react->mets()) {
                    if (i == j)
                        continue;
                    if (react->met_coeff(j) > 0) {
                        // Add an edge, and use the reaction idx as
                        // the edge index
                        graph.addEdge (i, j, 1, react->index());
                        DEBUG (
                            'T', "  Edge " << *scenario_->metabolite(i) <<
                            " " << *scenario_->metabolite(j) << " 1.0 " <<
                            react->name()
                        );
                    }
                }
            }
        }
    }
    DEBUG ('t', "  Added " << graph.num_edges() << " edges");
#ifndef NDEBUG
    if (Debug::doDebug('T')) {
        graph.show (Debug::debugFile());
    }
#endif
    // Find a path for each metabolite used in the biomass reaction
    set<const Reaction*> reacts;
    set<const Metabolite*> mets;
    /**
     * @struct Edge
     * @brief Represents a directed edge in the metabolic pathway graph
     *
     * Used for storing pathway connections with associated flux values
     * for visualization purposes.
     */
    struct Edge
    {
        std::string from; /**< Source node name (reaction or metabolite) */
        std::string to; /**< Destination node name (reaction or metabolite) */
        double flux; /**< Flux value through this edge */
        /**
         * @brief Constructs an Edge with specified source, destination, and flux
         *
         * @param from_ Source node name
         * @param to_ Destination node name
         * @param flux_ Flux value through the edge
         */
        Edge (string from_, string to_, double flux_) :
            from(from_), to(to_), flux(flux_) {}
    };
    list<Edge> edges;
    
    out << "digraph G {" << endl;
    for (auto from : from_mets) {
        for (auto to : to_mets) {
            DEBUG ('t', "Find a path from " << *from << " to " << *to);
            if (graph.findPath (from->index(), to->index())) {
                DEBUG ('t', "  Found a path!");
                
                mets.insert (from);
                size_t curr = to->index();
                const Reaction* prev = nullptr;
                while (curr != from->index()) {
                    DEBUG (
                        't', "  to " << curr <<
                        " = " << 
                        *(scenario_->metabolite(curr)) <<
                        " cost " << graph.costTo(curr) <<
                        " via " << graph.parent_edge(curr) <<
                        " = " <<
                        *(scenario_->reaction(graph.parent_edge(curr)))
                    );
                    auto met = scenario_->metabolite(curr);
                    auto react = scenario_->reaction(graph.parent_edge(curr));
                    reacts.insert (react);
                    mets.insert (met);
                    edges.push_back (
                        Edge (
                            react->name(), met->name(), flux_[react->index()]
                        )
                    );
                    if (prev != nullptr)
                        edges.push_back (
                            Edge (
                                met->name(), prev->name(), flux_[prev->index()]
                            )
                        );
                    prev = react;
                    
                    curr = graph.parent(curr);
                }
                if (prev != nullptr) {
                    edges.push_back (
                        Edge (
                            from->name(), prev->name(),
                            flux_[prev->index()]
                        )
                    );
                }
            }
            else {
                DEBUG ('t', "  No path");
            }
            graph.reset_costs();
        }
    }
    for (auto react : reacts) {
        out << "  " << react->name() <<
            " [shape=box color=red tooltip=\"" << react->formula() << "\"]" <<
            endl;
    }
    for (auto met : mets) {
        out << "  " << met->name() <<
            " [color=blue tooltip=\"" << met->full_name() << "\"]" << endl;
    }
    // Find max flux
    double max_flux = 1.0f;
    for (auto edge : edges) 
        if (edge.flux > max_flux) 
            max_flux = edge.flux;
    DEBUG ('t', "Max flux is " << max_flux);
    
    for (auto edge : edges) {
        double width = 1.0f + (edge.flux * 8.0f / max_flux);
        out << "  " << edge.from << " -> " << edge.to <<
            " [edgetooltip=\"Flux " << limeFormat("%g", edge.flux) <<
            "\" penwidth=" + to_string(width) + "]" <<
            endl;
    }
    out << "}" << endl;
}

/**
 * @brief Writes metabolite balance for all metabolites
 *
 * Outputs the net balance (supply - demand) for each metabolite, with
 * positive values indicating required supply and negative values indicating
 * residuals. Metabolites with zero balance are listed last.
 *
 * @param out Output stream for metabolite balance data
 */
void
Solution::write_metbal (std::ostream& out)
{
    out << "# <metabolite> <reqd supply(+ve)/resid(-ve)> lb ub" << endl;
    vector<double> sum (scenario_->num_metabolites());
    calc_metbal (sum);

    // List of mets with sum == 0
    list<string> zero_sum;
    for (size_t i = 0; i < scenario_->num_metabolites(); i++) {
        auto met = scenario_->metabolite(i);
        // negate everything because they are stored in reverse sign
        sum[i] = 0.0f - sum[i];

        if (sum[i] == 0.0f)
            zero_sum.push_back (met->name());
        else
            out << met->name() << " " <<
                sum[i] << 
                endl;
    }
    for (auto name : zero_sum) {
            out << name << " " <<
                0.0f << 
                endl;
    }
}

/**
 * @brief Calculates metabolite balance across all reactions
 *
 * Computes the net production/consumption for each metabolite by summing
 * stoichiometric coefficients times flux for all active non-dummy reactions.
 *
 * @param metbal Vector to store metabolite balance values (resized if needed)
 */
void
Solution::calc_metbal (std::vector<double>& metbal)
{
    if (metbal.size() < scenario_->num_metabolites())
        metbal.resize (scenario_->num_metabolites());
    std::fill (metbal.begin(), metbal.end(), 0.0f);
    
    for (size_t k = 0; k < scenario_->num_reactions(); k++) {
        auto react = scenario_->reaction(k);
        if (uses_react(k) && !react->is_dummy()) {
            for (size_t m : react->mets()) {
                metbal[m] += react->met_coeff(m) * flux(k);
            }
        }
    }
}

/**
 * @brief Calculates total production amount for each metabolite
 *
 * Computes the sum of positive (output) flux contributions for each metabolite
 * across all active non-dummy reactions.
 *
 * @param metprod Vector to store metabolite production values (resized if needed)
 */
void
Solution::calc_metprod (std::vector<double>& metprod)
{
    if (metprod.size() < scenario_->num_metabolites())
        metprod.resize (scenario_->num_metabolites());
    std::fill (metprod.begin(), metprod.end(), 0.0f);
    
    for (size_t k = 0; k < scenario_->num_reactions(); k++) {
        auto react = scenario_->reaction(k);
        if (uses_react(k) && !react->is_dummy()) {
            for (auto met : react->out_mets()) {
                metprod[met->index()] += react->met_coeff(met) * flux(k);
            }
        }
    }
}

/**
 * @brief Writes total metabolite usage (supply + production)
 *
 * Outputs the sum of supplied amounts and produced amounts for each metabolite.
 * Includes both external supply and internal production through reactions.
 *
 * @param out Output stream for metabolite usage data
 */
void
Solution::write_met_use (std::ostream& out)
{
    out << "# <metabolite> <total flux>" << endl;
    vector<double> sum (scenario_->num_metabolites(), 0.0f);
    
    // Calc metbal
    calc_metbal (sum);
    // But only keep supplied stuff.
    // Note signs are reversed, so supplied is -ve
    for (size_t i = 0; i < scenario_->num_metabolites(); i++) {
        if (sum[i] < 0.0f)
            sum[i] = 0.0f - sum[i];
        else
            sum[i] = 0.0f;
    }

    // Now add in +ve flows
    for (size_t k = 0; k < scenario_->num_reactions(); k++) {
        auto react = scenario_->reaction(k);
        if (uses_react(k) && !react->is_dummy()) {
            for (auto met : react->out_mets()) {
                sum[met->index()] += react->met_coeff(met) * flux(k);
            }
        }
    }

    // List of mets with sum == 0
    list<string> zero_sum;
    for (size_t i = 0; i < scenario_->num_metabolites(); i++) {
        auto met = scenario_->metabolite(i);

        if (sum[i] == 0.0f)
            zero_sum.push_back (met->name());
        else
            out << met->name() << " " <<
                sum[i] << 
                endl;
    }
    for (auto name : zero_sum) {
            out << name << " " <<
                0.0f << 
                endl;
    }
}

/**
 * @brief Draws a visual representation of the gap-fill solution
 *
 * Creates a 2D plot with reactions on x-axis and metabolites on y-axis.
 * Circle sizes represent flux magnitudes, colors indicate production (blue)
 * vs consumption (red). Includes supply and residual columns.
 *
 * @param dig Pointer to Dig visualization object
 */
void
Solution::draw (Dig* dig)
{
    // Use experiment 0 as bassis for supply/demand
    auto exp = scenario_->experiment(0);
    
    dig->title ("Gap Fill Sol");
    dig->xTic (0.0f, "Supply");
    for (size_t i = 0; i < scenario_->num_reactions(); i++) {
        double x = i + 1.0f;
        dig->xTic (x, scenario_->reaction(i)->name());
    }
    dig->xTic (scenario_->num_reactions()+1.0f, "Residual");

    for (size_t j = scenario_->num_metabolites(); j > 0; j--) {
        double y = scenario_->num_metabolites() - (j-1);
        dig->yTic (y, scenario_->metabolite(j-1)->name());
    }
    double y = scenario_->num_metabolites() + 1.0f;
    dig->yTic (y, "Flux");

    double base_diam = 10.0f;
    char buffer[100];
        
    // One col for each reaction
    for (size_t i = 0; i < scenario_->num_reactions(); i++) {
        auto react = scenario_->reaction(i);
        double x = i + 1.0f;
        
        if (!limeIsZero (flux(i))) {
            double y = scenario_->num_metabolites() + 1.0f;
            int diam = (int) (base_diam * sqrt(flux(i)));
            dig->style (Dig::DARK_GREEN, 0);
            dig->circle (x, y, diam, true);
            snprintf (buffer, 100, "%.4g", flux(i));
            if (!limeIsZero (flux(i)))
                dig->labelPoint (x, y, buffer);
        }
        
        // One row for each metabolite
        for (size_t j = 0; j < scenario_->num_metabolites(); j++) {
            auto met = scenario_->metabolite(j);
            double y = scenario_->num_metabolites() - j;

            double mult = 1;
            int colour = 0;
            if (limeIsZero (flux(i))) {
                if (react->met_coeff(met) < 0) {
                    mult = -react->met_coeff(met);
                    colour = Dig::StyleColours::PINK;
                }
                else {
                    mult = react->met_coeff(met);
                    colour = Dig::StyleColours::CADET_BLUE;
                }
            }
            else {
                if (react->met_coeff(met) < 0) {
                    mult = flux(i) * -react->met_coeff(met);
                    colour = Dig::StyleColours::RED;
                }
                else {
                    mult = flux(i) * react->met_coeff(met);
                    colour = Dig::StyleColours::BLUE;
                }
            }
            if (mult > 0) {
                dig->style (colour, 0);
                int diam = (int) (base_diam * sqrt(mult));
                dig->circle (x, y, diam, true);
                snprintf (buffer, 100, "%.4g", mult);
                if (!limeIsZero (flux(i)))
                    dig->labelPoint (x, y, buffer);
            }
        }
    }
        
    // Do supply/residual for each metabolite
    for (size_t j = 0; j < scenario_->num_metabolites(); j++) {
        auto met = scenario_->metabolite(j);
        double y = scenario_->num_metabolites() - j;
        if (exp->is_supply(met)) {
            double x = 0.0f;
            dig->style (Dig::StyleColours::BLUE, 0);
            int diam = (int) (base_diam * sqrt(exp->supply(met)));
            dig->circle (x, y, diam, true);
            snprintf (buffer, 100, "%.4g", exp->supply(met));
            dig->labelPoint (x, y, buffer);
        }
        if (exp->is_residual(met)) {
            double x = scenario_->num_reactions()+1.0f;
            dig->style (Dig::StyleColours::BLUE, 0);
            int diam = (int) (base_diam * sqrt(exp->residual(met)));
            dig->circle (x, y, diam, true);
            snprintf (buffer, 100, "%.4g", exp->residual(met));
            dig->labelPoint (x, y, buffer);
        }
    }
}

/**
 * @brief Draws metabolite balance visualization
 *
 * Creates a bar chart showing the net balance for each metabolite with
 * upper and lower bounds from the experiment. Positive values indicate
 * required supply, negative values indicate residuals.
 *
 * @param dig Pointer to Dig visualization object
 */
void
Solution::draw_metbal (Dig* dig)
{
    // Use experiment 0
    auto exp = scenario_->experiment(0);
    
    vector<double> metbal (scenario_->num_metabolites());
    calc_metbal (metbal);

    double min_y = 0.0;
    double max_y = 0.0;
    
    dig->title ("Metabolite Balance");
    for (size_t i = 0; i < scenario_->num_metabolites(); i++) {
        double x = i + 1.0f;
        dig->xTic (x, scenario_->metabolite(i)->name());
    }
    // Do upper and lower bounds
    dig->style(0, 0, 0);
    for (size_t i = 0; i < scenario_->num_metabolites(); i++) {
        auto met = scenario_->metabolite(i);
        
        double x = i + 1.0f;
        dig->draw (x-0.25, -exp->lb(met), x+0.25, -exp->lb(met));
        dig->draw (x-0.25, -exp->ub(met), x+0.25, -exp->ub(met));
    }
    dig->style(1);
    for (size_t i = 0; i < scenario_->num_metabolites(); i++) {
        auto met = scenario_->metabolite(i);
        
        double x = i + 1.0f;
        dig->draw (x, 0.0f, x, -metbal[i]);
        if (-metbal[i] < min_y)
            min_y = -metbal[i];
        if (-metbal[i] > max_y)
            max_y = -metbal[i];
    }
    dig->mark(0);
    if (-max_y < min_y)
        min_y = -max_y;
    dig->labelPoint (0, max_y + 1, "Supply");
    dig->labelPoint (0, min_y - 1, "Residual");
}

/**
 * @brief Draws reduced cost visualization for reactions
 *
 * Creates a chart showing objective costs for used reactions and reduced costs
 * for unused reactions, sorted by value. Helps identify which unused reactions
 * are closest to being included in the solution.
 *
 * @param dig Pointer to Dig visualization object
 */
void
Solution::draw_reduced_cost (Dig* dig)
{
    // First, sort by col val
    SortPair srt;
    for (size_t k = 0; k < scenario_->num_reactions(); k++) {
        // If used, sort to the front, in order
        if (uses_react(k)) {
            //double srt_val = ((int)i - (int)num_reactions());
            srt.add (k, -scenario_->reaction(k)->obj_coeff());
        }
        else
            srt.add (k, scenario_->reaction(k)->reduced_cost());
    }
    srt.doSort();
        
    // Now chart reduced cost (reduced_cost) of each reaction
        
    dig->title ("Column Cost - smaller is better");
    for (size_t i = 0; i < scenario_->num_reactions(); i++) {
        double x = i + 1.0f;
        dig->xTic (x, scenario_->reaction(srt[i])->name());
    }
        
    for (size_t i = 0; i < scenario_->num_reactions(); i++) {
        auto react = scenario_->reaction(srt[i]);
        double x = i + 1.0f;
        
        if (uses_react(srt[i])) {
            dig->style (Dig::MAGENTA, Dig::NO_MARK);
            dig->label ("Cost of used reaction");
            double y = react->obj_coeff();
            dig->draw (x, 0, x, y);
            //int diam = 10;
            //dig->circle (x, 0.0f, diam, true);
        }
        else {
            if (!react->is_dummy()) {
                dig->style (Dig::DARK_GREEN, Dig::NO_MARK, Dig::THICK);
                dig->label ("Reduced cost of unused reaction");
                double y = react->reduced_cost();
                dig->draw (x, 0, x, y);
            }
        }
    }
}






