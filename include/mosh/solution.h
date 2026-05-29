#pragma once

#include <vector>
#include <memory>
#include <set>

#include "lime/numutil.h"

#include "mosh/scenario.h"
#include "mosh/reaction.h"
#include "mosh/params.h"

namespace mosh
{
    class Solution 
    {
    public:
        Solution(
            const Scenario* scenario, const Params* params
        ) :
            scenario_(scenario),
            params_(params),
            flux_(scenario->num_reactions(), 0.0f)
        {
        }
        Solution(Solution* other) :
            scenario_(other->scenario_),
            params_(other->params_),
            flux_(other->flux_)
        {
        }

        const Scenario* scenario() const {return scenario_;}
        const Params* params() const {return params_;}
        
        virtual double flux (size_t k) const {return flux_[k];}
        double flux (const Reaction* react) const
        {
            return flux (react->index());
        }
        virtual void set_flux (size_t k, double val) {flux_[k] = val;}
        void set_flux (const Reaction* react, double val) 
        {
            set_flux (react->index(), val);
        }
        virtual bool uses_react(size_t k) const
        {
            return flux(k) > 0.0f;
        }
        virtual bool uses_react(const Reaction* react) const
        {
            return this->uses_react(react->index());
        }

        virtual void fill (double value)
        {
            std::fill (flux_.begin(), flux_.end(), value);
        }

        virtual double obj_value (double biomass_mult) const;
        virtual double rel_obj_value () const;
        virtual double abs_obj_value () const;
        virtual size_t num_reactions_used() const;
        virtual double sum_flux () const;
        virtual double max_cost () const;
        virtual double biomass_flux () const;
        virtual double biomass_obj (double biomass_mult) const;
        virtual double dummy_flux () const;
        bool has_dummy_flux () const
        {
            return !limeIsZero(dummy_flux());
        }
        virtual double dummy_obj_val () const;
        virtual int num_dummy () const;
        virtual void write_flux (std::ostream& out);

        // Check if the solution has run-away reactions (any experiment)
        bool is_runaway() const;
        // Check if the solution has run-away reactions
        // with regard to the carbon-sources in the experiment
        bool is_runaway(const Experiment* exp) const;
        
        int reaction_count () const;
        void calc_metbal (std::vector<double>& metbal);
        void calc_metprod (std::vector<double>& metprod);

        void write_model (std::ostream& out, bool cost_one);
        void write_met_dot (std::ostream& out);
        void write_metbal (std::ostream& out);
        void write_met_use (std::ostream& out);
        void write_carbon_source_dot (std::ostream& out);
        void write_supply_demand_dot (std::ostream& out);

        void draw (lime::Dig* dig);
        void draw_metbal (lime::Dig* dig);
        void draw_reduced_cost (lime::Dig* dig);

    protected:
        void plot_paths (
            std::ostream& out,
            std::set<const Metabolite*>& from_mets,
            std::set<const Metabolite*>& to_mets
        );
        
        const Scenario* scenario_;
        const Params* params_;
        std::vector<double> flux_;
    };

    using SolutionPtr = std::shared_ptr<Solution>;
}
