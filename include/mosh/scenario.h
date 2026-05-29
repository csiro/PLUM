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

namespace mosh
{
    class Solution;
    class MultiSol;
    class Params;

    using StrIdxMap = std::map<std::string,size_t>;
    using MultiSolPtr = std::shared_ptr<MultiSol>;

    enum ReactCostPolicy {USEMIN, USEMAX, REPLACE};
    
    class Scenario 
    {
    public:
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


        size_t num_metabolites() const {return metabolite_.size();}
        const Metabolite* metabolite(size_t k) const {
            return metabolite_[k].get();
        }
        void add_metabolite (MetabolitePtr met)
        {
            met->set_index (metabolite_.size());
            metabolite_.push_back (met);
            met_map_[met->name()] = met->index();
        }
        const std::vector<MetabolitePtr>& metabolites() const {
            return metabolite_;
        }
        void calc_sources (int exp_id = 0)
        {
            for (auto& met : metabolite_) {
                if (experiment(exp_id)->is_supply(met->index())) {
                    DEBUG ('f', "  " << *met << " is source (supply)");
                    met->set_is_source(true);
                }
            }
        }
        
        size_t num_reactions() const {return reaction_.size();}
        const std::vector<ReactionPtr>& reactions() {return reaction_;}
        Reaction* reaction(size_t k) const {return reaction_[k].get();}
        const std::vector<ReactionPtr>& reactions() const {return reaction_;}
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
        // Num reactions before adding dummies
        size_t orig_num_reactions() const {return orig_num_reactions_;}

        size_t num_experiments() const {return experiment_.size();}
        // Do we have multiple experiments?
        bool is_multi_exp() const {return experiment_.size() > 1;}
        Experiment* experiment(size_t k) const
        {
            return experiment_[k].get();
        }
        const std::vector<ExperimentPtr>& experiments() const {
            return experiment_;
        }

        const Reaction* biomass_react() const {return biomass_react_;}

        void finalise (Params* params);

        size_t biolog_rank0() const {return biolog_rank_[0]->index();}
        const Experiment* biolog_rank(size_t k) const
        {
            return biolog_rank_[k];
        }

        const Metabolite* dummy_ex_for (Metabolite* met)
        {
            return dummy_ex_for_[met];
        }
        const Metabolite* dummy_dm_for (Metabolite* met)
        {
            return dummy_dm_for_[met];
        }


        double max_react_cost() const {return max_react_cost_;}
        bool has_base_flux() const {return has_base_flux_;}
        void calc_target_flux (std::vector<double>& target, double rank0_flux);
        
        void read_data (std::string data_fn);
        void read_flux (std::string flux_fn);
        void read_react_cost (
            std::string react_cost_fn, const Params* params,
            ReactCostPolicy policy
        );
        void read_react_cost_list (std::string react_cost_list_fn);
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
        void read_c_sources (std::string c_source_fn);
        void read_cycle_mets (std::string cycle_met_fn);
        void read_supply_demand (
            std::string supply_demand_fn, const Params* params
        );
        void no_supply_demand ()
        {
            // Create an empty experiment
            ExperimentPtr experiment =
                std::make_shared<Experiment> ("None", 0, num_metabolites());
            experiment_.push_back (experiment);
        }
        void read_base_flux (std::string base_flux_fn);
        void add_dummy_reactions (double dummy_cost);
        void add_dummy_biomass_dm (Params* params);
        void preserve_dummies ();

        int depth (); // Overall problem depth
        int sol_depth (Solution* sol); // Depth for a given sol

        void count_reachability (
            size_t experiment, int& reactions, int& mets, int& residuals
        );
        int calc_reachability (
            size_t experiment, 
            std::vector<int>& react_depth, std::vector<bool>& avail,
            std::vector<int>& enabled_by,
            int& residual_depth, int& never, Solution* flux = nullptr
        );

        void select_all()
        {
            for (auto react : reaction_) {
                react->set_selected(react->is_active());
            }
        }
        size_t num_selected()
        {
            size_t count = 0;
            for (size_t k = 0; k < num_reactions(); k++)
                if (reaction(k)->is_selected())
                    count++;
            return count;
        }
        
        const Metabolite* find_metabolite (std::string name) 
        {
            if (met_map_.find (name) == met_map_.end())
                return nullptr;
            return metabolite(met_map_[name]);
        }
        Reaction* find_reaction (std::string name)
        {
            if (react_map_.find (name) == react_map_.end())
                return nullptr;
            return reaction(react_map_[name]);
        }
        
        Experiment* find_experiment (std::string name)
        {
            for (auto exp : experiment_) {
                if (exp->name().compare (name) == 0)
                    return exp.get();
            }
            return nullptr;
        }

        // Make disjkstra graph, and return idx of dummy source
        size_t make_dijkstra (
            lime::Dijkstra<double>& graph, const Params* params
        );
        // Supply an initial availability graph
        size_t make_dijkstra (
            lime::Dijkstra<double>& graph, const Params* params,
            std::vector<bool>& avail, std::vector<const Reaction*>& supplied_by
        );

        void draw_reachability (lime::Dig* dig, Solution* sol = nullptr);

    private:
        std::vector<MetabolitePtr> metabolite_;
        std::vector<ReactionPtr> reaction_;
        std::vector<ExperimentPtr> experiment_;

        Reaction* biomass_react_;
        // The biolog order of experiments (descending)
        // biolog_rank_[0] is the experiment with the highest biolog score
        std::vector<const Experiment*> biolog_rank_;

        std::map<const Metabolite*,const Metabolite*> dummy_ex_for_;
        std::map<const Metabolite*,const Metabolite*> dummy_dm_for_;

        size_t orig_num_reactions_;
        double max_react_cost_;
        bool has_base_flux_;

        StrIdxMap met_map_;
        StrIdxMap react_map_;
    };

    using ScenarioPtr = std::shared_ptr<Scenario>;

    using ostream_ptr = std::shared_ptr<std::ostream>;
}
