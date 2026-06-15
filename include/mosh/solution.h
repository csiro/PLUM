/**
 * @file solution.h
 * @brief Solution class for metabolic network flux balance analysis
 *
 * This file defines the Solution class which represents a flux distribution
 * solution for metabolic gap-filling and flux balance analysis problems.
 * It stores reaction fluxes and provides methods for objective value calculation,
 * metabolic balance analysis, and visualization.
 */

#pragma once

#include <vector>
#include <memory>
#include <set>

#include "lime/numutil.h"

#include "mosh/scenario.h"
#include "mosh/reaction.h"
#include "mosh/params.h"

/**
 * @brief Namespace for metabolic optimization and simulation handling
 */
namespace mosh
{
    /**
     * @class Solution
     * @brief Represents a flux distribution solution for metabolic networks
     *
     * The Solution class encapsulates a complete flux distribution across all reactions
     * in a metabolic network scenario. It provides methods for accessing and modifying
     * reaction fluxes, computing objective values, checking solution feasibility,
     * analyzing metabolic balance, and generating various output formats including
     * visualization and DOT graph representations.
     */
    class Solution 
    {
    public:
        /**
         * @brief Construct a new Solution object
         * @param scenario Pointer to the metabolic scenario containing network structure
         * @param params Pointer to algorithm parameters for optimization
         */
        Solution(
            const Scenario* scenario, const Params* params
        ) :
            scenario_(scenario),
            params_(params),
            flux_(scenario->num_reactions(), 0.0f)
        {
        }
        /**
         * @brief Copy constructor
         * @param other Pointer to another Solution object to copy from
         */
        Solution(Solution* other) :
            scenario_(other->scenario_),
            params_(other->params_),
            flux_(other->flux_)
        {
        }

        /**
         * @brief Get the metabolic scenario
         * @return const Scenario* Pointer to the scenario
         */
        const Scenario* scenario() const {return scenario_;}
        /**
         * @brief Get the algorithm parameters
         * @return const Params* Pointer to the parameters
         */
        const Params* params() const {return params_;}
        
        /**
         * @brief Get the flux value for a reaction by index
         * @param k Reaction index
         * @return double Flux value for the specified reaction
         */
        virtual double flux (size_t k) const {return flux_[k];}
        /**
         * @brief Get the flux value for a specific reaction
         * @param react Pointer to the reaction
         * @return double Flux value for the specified reaction
         */
        double flux (const Reaction* react) const
        {
            return flux (react->index());
        }
        /**
         * @brief Set the flux value for a reaction by index
         * @param k Reaction index
         * @param val New flux value
         */
        virtual void set_flux (size_t k, double val) {flux_[k] = val;}
        /**
         * @brief Set the flux value for a specific reaction
         * @param react Pointer to the reaction
         * @param val New flux value
         */
        void set_flux (const Reaction* react, double val)
        {
            set_flux (react->index(), val);
        }
        /**
         * @brief Check if a reaction is used (has positive flux) by index
         * @param k Reaction index
         * @return bool True if the reaction flux is greater than zero
         */
        virtual bool uses_react(size_t k) const
        {
            return flux(k) > 0.0f;
        }
        /**
         * @brief Check if a reaction is used (has positive flux)
         * @param react Pointer to the reaction
         * @return bool True if the reaction flux is greater than zero
         */
        virtual bool uses_react(const Reaction* react) const
        {
            return this->uses_react(react->index());
        }

        /**
         * @brief Fill all flux values with a constant value
         * @param value The value to set for all fluxes
         */
        virtual void fill (double value)
        {
            std::fill (flux_.begin(), flux_.end(), value);
        }

        /**
         * @brief Calculate the total objective value
         * @param biomass_mult Multiplier for biomass contribution to objective
         * @return double The computed objective value
         */
        virtual double obj_value (double biomass_mult) const;
        /**
         * @brief Calculate the relative objective value
         * @return double The relative objective value based on reaction costs
         */
        virtual double rel_obj_value () const;
        /**
         * @brief Calculate the absolute objective value
         * @return double The absolute objective value
         */
        virtual double abs_obj_value () const;
        /**
         * @brief Count the number of reactions with non-zero flux
         * @return size_t Number of active reactions in the solution
         */
        virtual size_t num_reactions_used() const;
        /**
         * @brief Calculate the sum of all flux values
         * @return double Total sum of fluxes across all reactions
         */
        virtual double sum_flux () const;
        /**
         * @brief Get the maximum cost among all reactions
         * @return double Maximum reaction cost in the solution
         */
        virtual double max_cost () const;
        /**
         * @brief Get the flux value through the biomass reaction
         * @return double Biomass reaction flux representing growth rate
         */
        virtual double biomass_flux () const;
        /**
         * @brief Calculate the biomass contribution to the objective function
         * @param biomass_mult Multiplier for biomass objective
         * @return double Biomass objective value
         */
        virtual double biomass_obj (double biomass_mult) const;
        /**
         * @brief Get the total flux through dummy reactions
         * @return double Sum of fluxes through dummy (artificial) reactions
         */
        virtual double dummy_flux () const;
        /**
         * @brief Check if the solution uses dummy reactions
         * @return bool True if dummy flux is non-zero
         */
        bool has_dummy_flux () const
        {
            return !limeIsZero(dummy_flux());
        }
        /**
         * @brief Calculate the objective value contribution from dummy reactions
         * @return double Dummy reaction objective value
         */
        virtual double dummy_obj_val () const;
        /**
         * @brief Count the number of active dummy reactions
         * @return int Number of dummy reactions with non-zero flux
         */
        virtual int num_dummy () const;
        /**
         * @brief Write flux values to an output stream
         * @param out Output stream to write to
         */
        virtual void write_flux (std::ostream& out);

        /**
         * @brief Check if the solution has run-away reactions in any experiment
         *
         * Run-away reactions are those with extremely high flux values that may
         * indicate unbounded or thermodynamically infeasible solutions.
         *
         * @return bool True if run-away reactions are detected
         */
        bool is_runaway() const;
        /**
         * @brief Check if the solution has run-away reactions for a specific experiment
         *
         * Checks for run-away reactions with regard to the carbon sources in the
         * specified experiment.
         *
         * @param exp Pointer to the experiment to check
         * @return bool True if run-away reactions are detected for this experiment
         */
        bool is_runaway(const Experiment* exp) const;
        
        /**
         * @brief Get the total number of reactions in the solution
         * @return int Number of reactions
         */
        int reaction_count () const;
        /**
         * @brief Calculate metabolic balance for all metabolites
         * @param metbal Output vector to store metabolite balance values
         */
        void calc_metbal (std::vector<double>& metbal);
        /**
         * @brief Calculate metabolite production rates
         * @param metprod Output vector to store metabolite production values
         */
        void calc_metprod (std::vector<double>& metprod);

        /**
         * @brief Write the optimization model to an output stream
         * @param out Output stream to write to
         * @param cost_one If true, use uniform costs of 1 for all reactions
         */
        void write_model (std::ostream& out, bool cost_one);
        /**
         * @brief Write metabolite network in DOT graph format
         * @param out Output stream to write DOT graph to
         */
        void write_met_dot (std::ostream& out);
        /**
         * @brief Write metabolic balance information to an output stream
         * @param out Output stream to write to
         */
        void write_metbal (std::ostream& out);
        /**
         * @brief Write metabolite usage information to an output stream
         * @param out Output stream to write to
         */
        void write_met_use (std::ostream& out);
        /**
         * @brief Write carbon source network in DOT graph format
         * @param out Output stream to write DOT graph to
         */
        void write_carbon_source_dot (std::ostream& out);
        /**
         * @brief Write supply-demand network in DOT graph format
         * @param out Output stream to write DOT graph to
         */
        void write_supply_demand_dot (std::ostream& out);

        /**
         * @brief Draw solution visualization using Dig graphics interface
         * @param dig Pointer to Dig graphics object
         */
        void draw (lime::Dig* dig);
        /**
         * @brief Draw metabolic balance visualization using Dig graphics interface
         * @param dig Pointer to Dig graphics object
         */
        void draw_metbal (lime::Dig* dig);
        /**
         * @brief Draw reduced cost visualization using Dig graphics interface
         * @param dig Pointer to Dig graphics object
         */
        void draw_reduced_cost (lime::Dig* dig);

    protected:
        /**
         * @brief Plot metabolic pathways between source and target metabolite sets
         * @param out Output stream to write pathway information to
         * @param from_mets Set of source metabolites
         * @param to_mets Set of target metabolites
         */
        void plot_paths (
            std::ostream& out,
            std::set<const Metabolite*>& from_mets,
            std::set<const Metabolite*>& to_mets
        );
        
        const Scenario* scenario_; /**< Pointer to the metabolic scenario */
        const Params* params_; /**< Pointer to algorithm parameters */
        std::vector<double> flux_; /**< Vector storing flux values for all reactions */
    };

    /**
     * @typedef SolutionPtr
     * @brief Shared pointer type for Solution objects
     */
    using SolutionPtr = std::shared_ptr<Solution>;
}
