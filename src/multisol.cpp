/**
 * @file multisol.cpp
 * @brief Implementation of MultiSol class for multi-experiment flux balance analysis solutions
 *
 * This file contains the implementation of methods for handling solutions across multiple
 * experimental conditions in metabolic network analysis. The MultiSol class aggregates
 * flux distributions from multiple experiments and provides unified metrics for gap-filling
 * optimization and metabolic flux analysis.
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

#include "mosh/multisol.h"
#include "mosh/metabolite.h"

using namespace std;
using namespace lime;
using namespace mosh;

/**
 * @brief Compute the total objective value across all experiments
 *
 * Sums the objective function values from all individual experiment solutions,
 * scaled by the biomass multiplier. This is used in gap-filling optimization to
 * evaluate the quality of a multi-experiment solution.
 *
 * @param biomass_mult Multiplier for biomass flux contribution to objective
 * @return Total objective value summed across all experiments
 */
double
MultiSol::obj_value (double biomass_mult) const
{
    double sum = 0.0f;
    for (size_t k = 0; k < scenario_->num_experiments(); k++)
        sum += sol_[k]->obj_value (biomass_mult);
    return sum;
}

/** Return the average rel obj across all runs */
double
MultiSol::rel_obj_value () const
{
    double sum = 0.0f;
    for (size_t k = 0; k < scenario_->num_experiments(); k++)
        sum += sol_[k]->rel_obj_value ();
    return sum / scenario_->num_experiments();
}

/**
 * @brief Calculate the absolute objective value based on reaction costs
 *
 * Computes the sum of objective coefficients for all reactions used in any experiment.
 * This provides a cost-based metric independent of flux magnitudes, useful for
 * evaluating the parsimony of gap-filling solutions.
 *
 * @return Sum of objective coefficients for all used reactions
 */
double
MultiSol::abs_obj_value () const
{
    double sum = 0.0f;
    for (size_t k = 0; k < scenario_->num_reactions(); k++)
        if (uses_react(k))
            sum += scenario_->reaction(k)->obj_coeff();
    return sum;
}

/**
 * @brief Find the maximum reaction cost among all used reactions
 *
 * Identifies the highest objective coefficient among reactions used in any experiment.
 * This can be used to identify the most costly reaction added during gap-filling.
 *
 * @return Maximum objective coefficient among used reactions
 */
double
MultiSol::max_cost () const
{
    double the_max = 0.0f;
    for (size_t k = 0; k < scenario_->num_reactions(); k++)
        if (uses_react(k) && scenario_->reaction(k)->obj_coeff() > the_max)
            the_max = scenario_->reaction(k)->obj_coeff();
    return the_max;
}

/**
 * @brief Count the total number of reactions with non-zero flux
 *
 * Counts reactions that are active (have non-zero flux) in any experiment.
 * This metric is useful for evaluating solution complexity and sparsity.
 *
 * @return Number of reactions used across all experiments
 */
size_t
MultiSol::num_reactions_used () const
{
    size_t count = 0;
    for (size_t k = 0; k < scenario_->num_reactions(); k++) {
        if (uses_react(k))
            count++;
    }
    return count;
}

/**
 * @brief Calculate the total sum of all flux magnitudes
 *
 * Sums the absolute flux values across all reactions and all experiments.
 * This provides a measure of total metabolic activity in the solution.
 *
 * @return Total flux sum across all experiments
 */
double
MultiSol::sum_flux () const
{
    double sum = 0.0f;
    for (size_t k = 0; k < scenario_->num_experiments(); k++)
        sum += sol_[k]->sum_flux ();
    return sum;
}

/**
 * @brief Calculate the total biomass production flux
 *
 * Sums the biomass reaction flux values across all experiments. This represents
 * the total growth or biomass production rate predicted by the model.
 *
 * @return Total biomass flux across all experiments
 */
double
MultiSol::biomass_flux () const
{
    double sum = 0.0f;
    for (size_t k = 0; k < scenario_->num_experiments(); k++)
        sum += sol_[k]->biomass_flux ();
    return sum;
}

/** Biomass flux times mult */
double
MultiSol::biomass_obj (double biomass_mult) const
{
    double sum = 0.0f;
    for (size_t k = 0; k < scenario_->num_experiments(); k++)
        sum += sol_[k]->biomass_obj (biomass_mult);
    return sum;
}

/** Sum of fluxes over dummy ractions */
double
MultiSol::dummy_flux () const
{
    double sum = 0.0f;
    for (size_t k = 0; k < scenario_->num_experiments(); k++)
        sum += sol_[k]->dummy_flux ();
    return sum;
}

/** Cost of dummy reactions */
double
MultiSol::dummy_obj_val () const
{
    double sum = 0.0f;
    for (size_t k = 0; k < scenario_->num_reactions(); k++)
        if (scenario_->reaction(k)->is_dummy() && uses_react(k))
            sum += scenario_->reaction(k)->obj_coeff();
    return sum;
}

/** Number of dummy reactions used */
int
MultiSol::num_dummy() const
{
    int count = 0;
    for (size_t k = 0; k < scenario_->num_reactions(); k++)
        if (scenario_->reaction(k)->is_dummy() && uses_react(k))
            count++;
    return count;
}

/**
 * @brief Count distinct reactions used in at least one experiment
 *
 * Counts the number of unique reactions that have non-zero flux in at least one
 * experiment. Unlike num_reactions_used(), this counts each reaction only once
 * regardless of how many experiments use it.
 *
 * @return Number of distinct reactions used across all experiments
 */
int
MultiSol::reaction_count () const
{
    int count = 0;
    for (auto& react_ptr : scenario_->reactions()) {
        const Reaction* react = react_ptr.get();
        for (size_t k = 0; k < scenario_->num_experiments(); k++) {
            if (sol_[k]->uses_react(react)) {
                count++;
                break;
            }
        }
    }
    return count;
}

/**
 * @brief Write flux distribution to output stream
 *
 * Outputs a formatted report of all reaction fluxes across experiments. For each
 * reaction used in any experiment, writes:
 * - Reaction name
 * - Sum of fluxes across experiments
 * - Individual flux values for each experiment
 * - Objective coefficient (cost)
 * - Metabolite stoichiometry
 * - Known flux constraint if applicable
 *
 * The output format is suitable for analysis and visualization of multi-experiment
 * flux balance analysis results.
 *
 * @param out Output stream to write flux data to
 */
void
MultiSol::write_flux (std::ostream& out) const
{
    out << "# <reaction> <flux>... <obj-coeff> {<metabolite> <coeff>}..." <<
        endl;
    out << "# Flux is the following order:" << endl;
    out << "# <sum>";
    for (size_t exp = 0; exp < scenario_->num_experiments(); exp++)
        out << " " << scenario_->experiment(exp)->name();
    out << endl;
    for (size_t k = 0; k < scenario_->num_reactions(); k++) {
        auto react = scenario_->reaction(k);
        if (uses_react(k)) {
            out << react->name() << " " << flux(react->index()) << " ";
            for (size_t exp = 0; exp < scenario_->num_experiments(); exp++) 
                out << " " << flux(react,exp);
                
            out << "  " << react->obj_coeff() << " ";
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









