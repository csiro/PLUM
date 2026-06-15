/**
 * @file scenario.h
 * @brief Core metabolic network scenario definition for gap-filling and flux balance analysis
 *
 * This file defines the Scenario class which represents a complete metabolic network model
 * including metabolites, reactions, experimental conditions, and biomass objectives.
 * It provides functionality for metabolic gap-filling, reachability analysis, and
 * flux balance analysis (FBA) in the context of genome-scale metabolic reconstruction.
 */
#pragma once

#include <vector>
#include <map>
#include <list>
#include <memory>

#include "lime/dig.h"
#include "lime/timekeeper.h"
#include "lime/numutil.h"
#include "lime/constants.h"
#include "lime/error.h"
#include "lime/dijkstra.h"

#include "mosh/constants.h"
#include "mosh/reaction.h"
#include "mosh/metabolite.h"
#include "mosh/experiment.h"

/**
 * @namespace mosh
 * @brief Metabolic optimization and scenario handling namespace
 */
namespace mosh
{
    class Solution;
    class MultiSol;
    class Params;

    /** @brief Map from string names to indices for fast lookup */
    using StrIdxMap = std::map<std::string,size_t>;
    /** @brief Shared pointer to MultiSol for managing multiple solutions */
    using MultiSolPtr = std::shared_ptr<MultiSol>;

    /**
     * @enum ReactCostPolicy
     * @brief Policy for handling reaction cost assignment when reading cost files
     *
     * @var ReactCostPolicy::USEMIN
     * Use minimum cost when multiple values are available
     * @var ReactCostPolicy::USEMAX
     * Use maximum cost when multiple values are available
     * @var ReactCostPolicy::REPLACE
     * Replace existing cost values with new ones
     */
    enum ReactCostPolicy {USEMIN, USEMAX, REPLACE};

    /**
     * @class Scenario
     * @brief Represents a complete metabolic network scenario for gap-filling analysis
     *
     * The Scenario class encapsulates all components of a metabolic network model including:
     * - Metabolites and their properties
     * - Reactions with stoichiometry and bounds
     * - Experimental conditions (media, supplies, demands)
     * - Biomass reaction and objective coefficients
     * - Dummy reactions for gap-filling
     *
     * This class provides methods for:
     * - Loading network data from files
     * - Computing metabolite reachability
     * - Managing reaction costs for optimization
     * - Adding dummy exchange and demand reactions
     * - Analyzing solution depth and network connectivity
     */
    class Scenario 
    {
    public:
        /**
         * @brief Default constructor initializing an empty scenario
         */
        Scenario () :
            metabolite_(),
            reaction_(),
            experiment_(),
            biomass_react_(nullptr),
            biolog_rank_(),
            dummy_ex_for_(),
            dummy_dm_for_(),
            orig_num_reactions_(0),
            max_react_cost_(1.0f),
            has_base_flux_(false),
            met_map_(),
            react_map_()
        {
        }

        /**
         * @brief Get the total number of metabolites in the scenario
         * @return Number of metabolites
         */
        size_t num_metabolites() const {return metabolite_.size();}
        /**
         * @brief Access a specific metabolite by index
         * @param k Index of the metabolite
         * @return Pointer to the metabolite (const)
         */
        const Metabolite* metabolite(size_t k) const {
            return metabolite_[k].get();
        }
        /**
         * @brief Add a metabolite to the scenario
         * @param met Shared pointer to the metabolite to add
         *
         * Sets the metabolite's index and updates the name-to-index map.
         */
        void add_metabolite (MetabolitePtr met)
        {
            met->set_index (metabolite_.size());
            metabolite_.push_back (met);
            met_map_[met->name()] = met->index();
        }
        /**
         * @brief Get all metabolites in the scenario
         * @return Const reference to vector of metabolite pointers
         */
        const std::vector<MetabolitePtr>& metabolites() const {
            return metabolite_;
        }
        /**
         * @brief Calculate which metabolites are sources (supplied) in an experiment
         * @param exp_id Experiment ID to check (default: 0)
         *
         * Marks metabolites as sources if they are supplied in the specified experiment.
         */
        void calc_sources (int exp_id = 0)
        {
            for (auto& met : metabolite_) {
                if (experiment(exp_id)->is_supply(met->index())) {
                    DEBUG ('f', "  " << *met << " is source (supply)");
                    met->set_is_source(true);
                }
            }
        }

        /**
         * @brief Get the total number of reactions in the scenario
         * @return Number of reactions (including dummy reactions if added)
         */
        size_t num_reactions() const {return reaction_.size();}
        /**
         * @brief Get all reactions in the scenario (non-const)
         * @return Reference to vector of reaction pointers
         */
        const std::vector<ReactionPtr>& reactions() {return reaction_;}
        /**
         * @brief Access a specific reaction by index
         * @param k Index of the reaction
         * @return Pointer to the reaction
         */
        Reaction* reaction(size_t k) const {return reaction_[k].get();}
        /**
         * @brief Get all reactions in the scenario (const)
         * @return Const reference to vector of reaction pointers
         */
        const std::vector<ReactionPtr>& reactions() const {return reaction_;}
        /**
         * @brief Add a reaction to the scenario
         * @param react Shared pointer to the reaction to add
         *
         * Sets the reaction's index, updates the name-to-index map, identifies biomass
         * reactions, and tracks maximum reaction cost.
         */
        void add_reaction (ReactionPtr react)
        {
            react->set_index (reaction_.size());
            reaction_.push_back (react);
            react_map_[react->name()] = react->index();
            if (react->is_biomass()) {
                if (biomass_react_ != nullptr) {
                    limeCrash (
                        "Too many biomass reactions adding " << react->name() <<
                        " already have " << biomass_react_->name()
                    );
                }
                biomass_react_ = react.get();
            }
            if (react->obj_coeff() > max_react_cost_)
                max_react_cost_ = react->obj_coeff();
        }
        /**
         * @brief Get the original number of reactions before adding dummy reactions
         * @return Number of original reactions (excluding dummies)
         */
        size_t orig_num_reactions() const {return orig_num_reactions_;}

        /**
         * @brief Get the total number of experimental conditions
         * @return Number of experiments
         */
        size_t num_experiments() const {return experiment_.size();}
        /**
         * @brief Check if multiple experimental conditions are defined
         * @return True if more than one experiment exists
         */
        bool is_multi_exp() const {return experiment_.size() > 1;}
        /**
         * @brief Access a specific experiment by index
         * @param k Index of the experiment
         * @return Pointer to the experiment
         */
        Experiment* experiment(size_t k) const
        {
            return experiment_[k].get();
        }
        /**
         * @brief Get all experimental conditions
         * @return Const reference to vector of experiment pointers
         */
        const std::vector<ExperimentPtr>& experiments() const {
            return experiment_;
        }

        /**
         * @brief Get the biomass reaction (growth objective)
         * @return Pointer to the biomass reaction (const)
         */
        const Reaction* biomass_react() const {return biomass_react_;}

        /**
         * @brief Finalize scenario setup after loading all data
         * @param params Pointer to parameter configuration
         *
         * Performs final initialization steps including ranking experiments by biolog score.
         */
        void finalise (Params* params);

        /**
         * @brief Get the index of the highest-ranked biolog experiment
         * @return Index of the top-ranked experiment
         */
        size_t biolog_rank0() const {return biolog_rank_[0]->index();}
        /**
         * @brief Get a biolog-ranked experiment by rank position
         * @param k Rank position (0 = highest biolog score)
         * @return Pointer to the experiment at that rank
         */
        const Experiment* biolog_rank(size_t k) const
        {
            return biolog_rank_[k];
        }

        /**
         * @brief Get the dummy exchange metabolite for a given metabolite
         * @param met Pointer to the metabolite
         * @return Pointer to the corresponding dummy exchange metabolite
         */
        const Metabolite* dummy_ex_for (Metabolite* met)
        {
            return dummy_ex_for_[met];
        }
        /**
         * @brief Get the dummy demand metabolite for a given metabolite
         * @param met Pointer to the metabolite
         * @return Pointer to the corresponding dummy demand metabolite
         */
        const Metabolite* dummy_dm_for (Metabolite* met)
        {
            return dummy_dm_for_[met];
        }


        /**
         * @brief Get the maximum reaction cost in the scenario
         * @return Maximum objective coefficient across all reactions
         */
        double max_react_cost() const {return max_react_cost_;}
        /**
         * @brief Check if base flux data has been loaded
         * @return True if base flux distribution is available
         */
        bool has_base_flux() const {return has_base_flux_;}
        /**
         * @brief Calculate target flux distribution for optimization
         * @param target Output vector to store target flux values
         * @param rank0_flux Flux value for the top-ranked experiment
         */
        void calc_target_flux (std::vector<double>& target, double rank0_flux);

        /**
         * @brief Read metabolic network data from file
         * @param data_fn Filename containing network definition (reactions, metabolites)
         */
        void read_data (std::string data_fn);
        /**
         * @brief Read flux distribution data from file
         * @param flux_fn Filename containing flux values
         */
        void read_flux (std::string flux_fn);
        /**
         * @brief Read reaction cost (objective coefficient) data from file
         * @param react_cost_fn Filename containing reaction costs
         * @param params Pointer to parameter configuration
         * @param policy Policy for handling multiple cost values (USEMIN, USEMAX, REPLACE)
         */
        void read_react_cost (
            std::string react_cost_fn, const Params* params,
            ReactCostPolicy policy
        );
        /**
         * @brief Read list of reaction cost files to process
         * @param react_cost_list_fn Filename containing list of cost files
         */
        void read_react_cost_list (std::string react_cost_list_fn);
        /**
         * @brief Set uniform cost for all non-biomass reactions
         * @param cost Cost value to assign
         *
         * Updates objective coefficients for all reactions except biomass.
         */
        void set_react_cost (double cost)
        {
            int count = 0;
            for (auto react : reactions()) {
                if (
                    !react->is_biomass() &&
                    !limeDblEqual (react->obj_coeff(), cost)
                ) {
                    react->set_obj_coeff (cost);
                    count++;
                }
            }
            std::cout << "Set cost on " << count <<
                " reactions to " << cost << std::endl;
        }
        /**
         * @brief Deselect reactions with cost above threshold
         * @param max_cost Maximum allowed cost
         *
         * Marks reactions exceeding max_cost as unselected for optimization.
         */
        void unselect_above_cost (double max_cost)
        {
            DEBUG ('A', "Set max react cost to " << max_cost);
            int count = 0;
            for (auto react : reactions()) {
                if (
                    !react->is_dummy() &&
                    react->is_selected() &&
                    react->obj_coeff() > max_cost
                ) {
                    react->set_selected(false);
                    count++;
                }
            }
            std::cout << "Deselected " << count <<
                " reactions with cost > " << max_cost << std::endl;
        }
        /**
         * @brief Read carbon source definitions from file
         * @param c_source_fn Filename containing carbon source metabolites
         */
        void read_c_sources (std::string c_source_fn);
        /**
         * @brief Read metabolites involved in metabolic cycles
         * @param cycle_met_fn Filename containing cycle metabolite names
         */
        void read_cycle_mets (std::string cycle_met_fn);
        /**
         * @brief Read supply and demand constraints for experimental conditions
         * @param supply_demand_fn Filename containing supply/demand specifications
         * @param params Pointer to parameter configuration
         */
        void read_supply_demand (
            std::string supply_demand_fn, const Params* params
        );
        /**
         * @brief Create an empty experiment with no supply/demand constraints
         *
         * Used when no experimental conditions are specified.
         */
        void no_supply_demand ()
        {
            // Create an empty experiment
            ExperimentPtr experiment =
                std::make_shared<Experiment> ("None", 0, num_metabolites());
            experiment_.push_back (experiment);
        }
        /**
         * @brief Read base flux distribution from file
         * @param base_flux_fn Filename containing base flux values
         *
         * Sets has_base_flux_ flag upon successful load.
         */
        void read_base_flux (std::string base_flux_fn);
        /**
         * @brief Add dummy exchange and demand reactions for gap-filling
         * @param dummy_cost Cost assigned to dummy reactions
         *
         * Creates artificial reactions to enable metabolite production/consumption.
         */
        void add_dummy_reactions (double dummy_cost);
        /**
         * @brief Add dummy demand reactions for biomass components
         * @param params Pointer to parameter configuration
         */
        void add_dummy_biomass_dm (Params* params);
        /**
         * @brief Mark dummy reactions for preservation in solutions
         *
         * Ensures dummy reactions are maintained during optimization.
         */
        void preserve_dummies ();

        /**
         * @brief Calculate overall metabolic network depth
         * @return Maximum depth across all metabolites and reactions
         *
         * Depth represents the longest pathway from sources to products.
         */
        int depth (); // Overall problem depth
        /**
         * @brief Calculate network depth for a specific solution
         * @param sol Pointer to the solution to analyze
         * @return Depth of the given solution
         */
        int sol_depth (Solution* sol); // Depth for a given sol

        /**
         * @brief Count reachable reactions, metabolites, and residuals
         * @param experiment Index of the experiment to analyze
         * @param reactions Output: number of reachable reactions
         * @param mets Output: number of reachable metabolites
         * @param residuals Output: number of unreachable metabolites
         */
        void count_reachability (
            size_t experiment, int& reactions, int& mets, int& residuals
        );
        /**
         * @brief Calculate detailed metabolite reachability analysis
         * @param experiment Index of the experiment to analyze
         * @param react_depth Output: depth at which each reaction becomes active
         * @param avail Output: availability status for each metabolite
         * @param enabled_by Output: which reaction enables each metabolite
         * @param residual_depth Output: maximum depth of unreachable metabolites
         * @param never Output: count of never-reachable metabolites
         * @param flux Optional flux distribution to constrain analysis
         * @return Maximum reachability depth
         */
        int calc_reachability (
            size_t experiment, 
            std::vector<int>& react_depth, std::vector<bool>& avail,
            std::vector<int>& enabled_by,
            int& residual_depth, int& never, Solution* flux = nullptr
        );

        /**
         * @brief Select all active reactions for optimization
         *
         * Marks all active reactions as selected.
         */
        void select_all()
        {
            for (auto react : reaction_) {
                react->set_selected(react->is_active());
            }
        }
        /**
         * @brief Count the number of selected reactions
         * @return Number of reactions marked as selected
         */
        size_t num_selected()
        {
            size_t count = 0;
            for (size_t k = 0; k < num_reactions(); k++)
                if (reaction(k)->is_selected())
                    count++;
            return count;
        }

        /**
         * @brief Find a metabolite by name
         * @param name Name of the metabolite to find
         * @return Pointer to the metabolite, or nullptr if not found
         */
        const Metabolite* find_metabolite (std::string name) 
        {
            if (met_map_.find (name) == met_map_.end())
                return nullptr;
            return metabolite(met_map_[name]);
        }
        /**
         * @brief Find a reaction by name
         * @param name Name of the reaction to find
         * @return Pointer to the reaction, or nullptr if not found
         */
        Reaction* find_reaction (std::string name)
        {
            if (react_map_.find (name) == react_map_.end())
                return nullptr;
            return reaction(react_map_[name]);
        }

        /**
         * @brief Find an experiment by name
         * @param name Name of the experiment to find
         * @return Pointer to the experiment, or nullptr if not found
         */
        Experiment* find_experiment (std::string name)
        {
            for (auto exp : experiment_) {
                if (exp->name().compare (name) == 0)
                    return exp.get();
            }
            return nullptr;
        }

        /**
         * @brief Construct a Dijkstra graph for shortest path analysis
         * @param graph Output Dijkstra graph structure
         * @param params Pointer to parameter configuration
         * @return Index of the dummy source node
         *
         * Creates a graph representation for computing metabolite distances.
         */
        size_t make_dijkstra (
            lime::Dijkstra<double>& graph, const Params* params
        );
        /**
         * @brief Construct a Dijkstra graph with initial availability constraints
         * @param graph Output Dijkstra graph structure
         * @param params Pointer to parameter configuration
         * @param avail Initial metabolite availability flags
         * @param supplied_by Reactions that supply each metabolite
         * @return Index of the dummy source node
         */
        size_t make_dijkstra (
            lime::Dijkstra<double>& graph, const Params* params,
            std::vector<bool>& avail, std::vector<const Reaction*>& supplied_by
        );

        /**
         * @brief Generate visual representation of metabolite reachability
         * @param dig Pointer to graphical output handler
         * @param sol Optional solution to visualize (default: nullptr)
         */
        void draw_reachability (lime::Dig* dig, Solution* sol = nullptr);

    private:
        std::vector<MetabolitePtr> metabolite_; /**< Collection of all metabolites in the network */
        std::vector<ReactionPtr> reaction_; /**< Collection of all reactions in the network */
        std::vector<ExperimentPtr> experiment_; /**< Collection of experimental conditions */

        Reaction* biomass_react_; /**< Pointer to the biomass/growth objective reaction */
        // The biolog order of experiments (descending)
        // biolog_rank_[0] is the experiment with the highest biolog score
        std::vector<const Experiment*> biolog_rank_; /**< Experiments ranked by biolog score (descending) */

        std::map<const Metabolite*,const Metabolite*> dummy_ex_for_; /**< Map from metabolites to their dummy exchange metabolites */
        std::map<const Metabolite*,const Metabolite*> dummy_dm_for_; /**< Map from metabolites to their dummy demand metabolites */

        size_t orig_num_reactions_; /**< Original reaction count before adding dummies */
        double max_react_cost_; /**< Maximum objective coefficient across all reactions */
        bool has_base_flux_; /**< Flag indicating whether base flux data has been loaded */

        StrIdxMap met_map_; /**< Map from metabolite names to indices for fast lookup */
        StrIdxMap react_map_; /**< Map from reaction names to indices for fast lookup */"
    };

    /** @brief Shared pointer to Scenario for memory management */
    using ScenarioPtr = std::shared_ptr<Scenario>;

    /** @brief Shared pointer to output stream for flexible output handling */
    using ostream_ptr = std::shared_ptr<std::ostream>;
}
