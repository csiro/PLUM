/**
 * @file multisol.h
 * @brief Multi-experiment solution container for metabolic flux balance analysis
 *
 * This file defines the MultiSol class which manages solutions across multiple
 * experimental conditions in flux balance analysis (FBA) scenarios. It provides
 * an aggregate view of fluxes and reactions across all experiments.
 */
#pragma once

#include <vector>
#include <memory>
#include <set>

#include "lime/numutil.h"

#include "mosh/scenario.h"
#include "mosh/solution.h"
#include "mosh/reaction.h"
#include "mosh/params.h"

/**
 * @brief Namespace for metabolic optimization and solution handling
 */
namespace mosh
{
    /**
     * @class MultiSol
     * @brief Container for managing solutions across multiple experimental conditions
     *
     * MultiSol extends Solution to handle metabolic flux distributions across multiple
     * experiments simultaneously. It aggregates flux values, reaction usage, and
     * objective values from individual solutions for each experimental condition.
     * This is particularly useful for gap-filling problems where multiple growth
     * conditions must be satisfied simultaneously.
     */
    class MultiSol : public Solution
    {
    public:
        /**
         * @brief Constructs a MultiSol with solutions for each experiment
         * @param scenario Pointer to the scenario defining experimental conditions
         * @param params Pointer to optimization parameters
         */
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

        /**
         * @brief Copy constructor - creates a deep copy of another MultiSol
         * @param other The MultiSol instance to copy
         */
        MultiSol (const MultiSol& other) :
            Solution (other.scenario_, other.params_),
            sol_(other.scenario_->num_experiments(), nullptr)
        {
            for (size_t exp = 0; exp < sol_.size(); exp++)
                sol_[exp] = std::make_shared<Solution> (*(other.sol_[exp]));
        }
        /**
         * @brief Constructs a MultiSol from a single solution
         * @param other Shared pointer to a single Solution to wrap
         * @note Asserts that the scenario has exactly one experiment
         */
        MultiSol (const SolutionPtr other) :
            Solution (other->scenario(), other->params()),
            sol_(other->scenario()->num_experiments(), nullptr)
        {
            assert (num_exp() == 1);
            sol_[0] = other;
        }

        /**
         * @brief Copies solution pointers from another MultiSol
         * @param other The MultiSol instance to copy from
         */
        void copy (const MultiSol& other)
        {
            for (size_t k = 0; k < sol_.size(); k++)
                sol_[k] = other.sol_[k];
        }
        /**
         * @brief Clears all solution pointers, setting them to nullptr
         */
        void clear ()
        {
            for (size_t k = 0; k < sol_.size(); k++)
                sol_[k] = nullptr;
        }


        /**
         * @brief Retrieves the solution for a specific experiment
         * @param exp Experiment index
         * @return Shared pointer to the Solution for the specified experiment
         */
        SolutionPtr sol (size_t exp)
        {
            return sol_[exp];
        }
        /**
         * @brief Retrieves the solution for a specific experiment (const version)
         * @param exp Experiment index
         * @return Shared pointer to the Solution for the specified experiment
         */
        SolutionPtr sol (size_t exp) const
        {
            return sol_[exp];
        }
        /**
         * @brief Sets the solution for a specific experiment
         * @param exp Experiment index
         * @param sol Shared pointer to the Solution to assign
         */
        void set_sol (size_t exp, SolutionPtr sol)
        {
            sol_[exp] = sol;
        }

        /**
         * @brief Returns the number of experiments
         * @return Number of experimental conditions in this MultiSol
         */
        size_t num_exp() const
        {
            return sol_.size();
        }

        /**
         * @brief Computes the total flux for a reaction across all experiments
         * @param k Reaction index
         * @return Sum of flux values for reaction k across all experiments
         */
        double flux (size_t k) const override
        {
            double sum = 0.0f;
            for (size_t exp = 0; exp < num_exp(); exp++)
                sum += sol_[exp]->flux(k);
            return sum;
        }
        /**
         * @brief Attempts to set flux directly on MultiSol (not supported)
         * @param k Reaction index
         * @param val Flux value
         * @throws Always crashes with error message - use experiment-specific set_flux instead
         */
        void set_flux (size_t k, double val) override
        {
            limeCrash ("Called set_flux on a multisol");
        }
        /**
         * @brief Checks if a reaction is used in any experiment
         * @param k Reaction index
         * @return True if reaction k has non-zero flux in any experiment
         */
        bool uses_react(size_t k) const override
        {
            for (size_t exp = 0; exp < num_exp(); exp++)
                if (sol_[exp]->uses_react (k))
                    return true;
            return false;
        }
        /**
         * @brief Checks if a reaction is used in any experiment
         * @param react Pointer to the Reaction object
         * @return True if the reaction has non-zero flux in any experiment
         */
        bool uses_react(const Reaction* react) const override
        {
            return this->uses_react(react->index());
        }
        /**
         * @brief Fills all fluxes in all experiments with a constant value
         * @param value The value to assign to all fluxes
         */
        void fill (double value) override
        {
            for (size_t exp = 0; exp < num_exp(); exp++)
                sol_[exp]->fill(value);
        }

        /**
         * @brief Gets the flux for a specific reaction in a specific experiment
         * @param k Reaction index
         * @param exp Experiment index
         * @return Flux value for reaction k in experiment exp
         */
        double flux (size_t k, size_t exp) const {
            return sol_[exp]->flux(k);
        }
        /**
         * @brief Gets the flux for a specific reaction in a specific experiment
         * @param react Pointer to the Reaction object
         * @param exp Experiment index
         * @return Flux value for the reaction in experiment exp
         */
        double flux (const Reaction* react, size_t exp) const
        {
            return flux (react->index(), exp);
        }

        /**
         * @brief Finds the maximum flux for a reaction across all experiments
         * @param react Pointer to the Reaction object
         * @return Maximum flux value for the reaction across all experiments
         */
        double max_flux (const Reaction* react) const
        {
            double max_val = 0.0f;
            for (size_t exp = 0; exp < num_exp(); exp++)
                if (flux (react->index(), exp) > max_val)
                    max_val = flux (react->index(), exp);
            return max_val;
        }

        /**
         * @brief Sets the flux for a specific reaction in a specific experiment
         * @param k Reaction index
         * @param exp Experiment index
         * @param val Flux value to assign
         */
        void set_flux (size_t k, size_t exp, double val)
        {
            sol_[exp]->set_flux (k, val);
        }
        /**
         * @brief Sets the flux for a specific reaction in a specific experiment
         * @param react Pointer to the Reaction object
         * @param exp Experiment index
         * @param val Flux value to assign
         */
        void set_flux (const Reaction* react, size_t exp, double val)
        {
            set_flux (react->index(), exp, val);
        }

        /**
         * @brief Gets the biomass production flux for a specific experiment
         * @param exp Experiment index
         * @return Biomass flux value for the specified experiment, or 0.0 if solution is null
         */
        double biomass_flux (size_t exp) const
        {
            if (sol_[exp] == nullptr)
                return 0.0f;
            return sol_[exp]->biomass_flux();
        }

        /**
         * @brief Computes the objective value with biomass multiplier
         * @param biomass_mult Multiplier for biomass contribution to objective
         * @return Total objective value across all experiments
         */
        double obj_value (double biomass_mult) const override;
        /**
         * @brief Computes the relative objective value across all experiments
         * @return Relative objective value (implementation defined in source file)
         */
        double rel_obj_value () const override;
        /**
         * @brief Computes the absolute objective value across all experiments
         * @return Absolute objective value (implementation defined in source file)
         */
        double abs_obj_value () const override;
        /**
         * @brief Computes the maximum cost across all experiments
         * @return Maximum cost value (implementation defined in source file)
         */
        double max_cost () const override;
        /**
         * @brief Counts the number of reactions used across all experiments
         * @return Total number of reactions with non-zero flux in any experiment
         */
        size_t num_reactions_used() const override;
        /**
         * @brief Computes the sum of all fluxes across all experiments
         * @return Total sum of flux magnitudes
         */
        double sum_flux () const override;
        /**
         * @brief Computes total biomass flux across all experiments
         * @return Sum of biomass production fluxes from all experiments
         */
        double biomass_flux () const override;
        /**
         * @brief Computes biomass contribution to objective across all experiments
         * @param biomass_mult Multiplier for biomass contribution
         * @return Total biomass objective value
         */
        double biomass_obj (double biomass_mult) const override;
        /**
         * @brief Computes total flux through dummy reactions across all experiments
         * @return Sum of dummy reaction fluxes (used for gap-filling)
         */
        double dummy_flux () const override;
        /**
         * @brief Computes objective contribution from dummy reactions
         * @return Total dummy reaction objective value across all experiments
         */
        double dummy_obj_val () const override;
        /**
         * @brief Counts the number of dummy reactions used across all experiments
         * @return Total count of active dummy reactions (used in gap-filling)
         */
        int num_dummy () const override;

        /**
         * @brief Counts experiments with runaway (unbounded) solutions
         * @return Number of experiments where the solution is unbounded
         */
        int num_runaways() const
        {
            int count = 0;
            for (size_t exp = 0; exp < num_exp(); exp++)
                if (sol_[exp] != nullptr && sol_[exp]->is_runaway())
                    count++;
            return count;
        }

        /**
         * @brief Counts total reaction usage across all experiments
         * @return Total reaction count (implementation defined in source file)
         */
        int reaction_count () const;

        /**
         * @brief Writes flux distribution to an output stream
         * @param out Output stream to write to
         */
        void write_flux (std::ostream& out) const;

    private:
        /** @brief Vector of solution pointers, one per experimental condition */
        std::vector<SolutionPtr> sol_;
    };

    /** @brief Shared pointer type for MultiSol instances */
    using MultiSolPtr = std::shared_ptr<MultiSol>;
}
