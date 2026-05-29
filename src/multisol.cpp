
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

double
MultiSol::abs_obj_value () const
{
    double sum = 0.0f;
    for (size_t k = 0; k < scenario_->num_reactions(); k++)
        if (uses_react(k))
            sum += scenario_->reaction(k)->obj_coeff();
    return sum;
}

double
MultiSol::max_cost () const
{
    double the_max = 0.0f;
    for (size_t k = 0; k < scenario_->num_reactions(); k++)
        if (uses_react(k) && scenario_->reaction(k)->obj_coeff() > the_max)
            the_max = scenario_->reaction(k)->obj_coeff();
    return the_max;
}

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

double
MultiSol::sum_flux () const
{
    double sum = 0.0f;
    for (size_t k = 0; k < scenario_->num_experiments(); k++)
        sum += sol_[k]->sum_flux ();
    return sum;
}

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









