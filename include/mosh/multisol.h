#pragma once

#include <vector>
#include <memory>
#include <set>

#include "lime/numutil.h"

#include "mosh/scenario.h"
#include "mosh/solution.h"
#include "mosh/reaction.h"
#include "mosh/params.h"

namespace mosh
{
    class MultiSol : public Solution
    {
    public:
        MultiSol(
            const Scenario* scenario, Params* params
        ) :
            Solution (scenario, params),
            sol_(scenario->num_experiments(), nullptr)
        {
            for (size_t exp = 0; exp < sol_.size(); exp++)
                sol_[exp] =
                    std::make_shared<Solution> (
                        scenario, params
                    );
        }

        MultiSol (const MultiSol& other) :
            Solution (other.scenario_, other.params_),
            sol_(other.scenario_->num_experiments(), nullptr)
        {
            for (size_t exp = 0; exp < sol_.size(); exp++)
                sol_[exp] = std::make_shared<Solution> (*(other.sol_[exp]));
        }
        /** Make a multi-sol form a single sol */
        MultiSol (const SolutionPtr other) :
            Solution (other->scenario(), other->params()),
            sol_(other->scenario()->num_experiments(), nullptr)
        {
            assert (num_exp() == 1);
            sol_[0] = other;
        }

        void copy (const MultiSol& other)
        {
            for (size_t k = 0; k < sol_.size(); k++)
                sol_[k] = other.sol_[k];
        }
        void clear ()
        {
            for (size_t k = 0; k < sol_.size(); k++)
                sol_[k] = nullptr;
        }
        
        SolutionPtr sol (size_t exp)
        {
            return sol_[exp];
        }
        SolutionPtr sol (size_t exp) const
        {
            return sol_[exp];
        }
        void set_sol (size_t exp, SolutionPtr sol)
        {
            sol_[exp] = sol;
        }

        size_t num_exp() const
        {
            return sol_.size();
        }

        double flux (size_t k) const override
        {
            double sum = 0.0f;
            for (size_t exp = 0; exp < num_exp(); exp++)
                sum += sol_[exp]->flux(k);
            return sum;
        }
        void set_flux (size_t k, double val) override
        {
            limeCrash ("Called set_flux on a multisol");
        }
        bool uses_react(size_t k) const override
        {
            for (size_t exp = 0; exp < num_exp(); exp++)
                if (sol_[exp]->uses_react (k))
                    return true;
            return false;
        }
        bool uses_react(const Reaction* react) const override
        {
            return this->uses_react(react->index());
        }
        void fill (double value) override
        {
            for (size_t exp = 0; exp < num_exp(); exp++)
                sol_[exp]->fill(value);
        }

        double flux (size_t k, size_t exp) const {
            return sol_[exp]->flux(k);
        }
        double flux (const Reaction* react, size_t exp) const
        {
            return flux (react->index(), exp);
        }
        double max_flux (const Reaction* react) const
        {
            double max_val = 0.0f;
            for (size_t exp = 0; exp < num_exp(); exp++)
                if (flux (react->index(), exp) > max_val)
                    max_val = flux (react->index(), exp);
            return max_val;
        }
        void set_flux (size_t k, size_t exp, double val)
        {
            sol_[exp]->set_flux (k, val);
        }
        void set_flux (const Reaction* react, size_t exp, double val) 
        {
            set_flux (react->index(), exp, val);
        }

        double biomass_flux (size_t exp) const
        {
            if (sol_[exp] == nullptr)
                return 0.0f;
            return sol_[exp]->biomass_flux();
        }
        
        double obj_value (double biomass_mult) const override;
        double rel_obj_value () const override;
        double abs_obj_value () const override;
        double max_cost () const override;
        size_t num_reactions_used() const override;
        double sum_flux () const override;
        double biomass_flux () const override;
        double biomass_obj (double biomass_mult) const override;
        double dummy_flux () const override;
        double dummy_obj_val () const override;
        int num_dummy () const override;

        int num_runaways() const
        {
            int count = 0;
            for (size_t exp = 0; exp < num_exp(); exp++)
                if (sol_[exp] != nullptr && sol_[exp]->is_runaway())
                    count++;
            return count;
        }
        
        int reaction_count () const;

        void write_flux (std::ostream& out) const;

    private:
        std::vector<SolutionPtr> sol_;
    };

    using MultiSolPtr = std::shared_ptr<MultiSol>;
}
