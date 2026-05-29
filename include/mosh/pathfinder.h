#pragma once

#include <list>
#include <set>
#include <memory>


#include "lime/dijkstra.h"

#include "mosh/scenario.h"


namespace mosh
{
    using ofstream_ptr = std::shared_ptr<std::ofstream>;
    
    class PathFinder {
        
    public:
        PathFinder (Scenario* scenario, const Params* params) :
            scenario_(scenario),
            params_(params),
            source_(
                std::make_shared<Reaction> (
                    "Source", 0.0f, 0.0f, false, false, true
                )
            ),
            supplied_by_(
                scenario->num_metabolites(), nullptr
            ),
            graph_(scenario->num_metabolites() + 1),
            source_idx_(0),
            in_loop_(scenario->num_metabolites(), false),
            dot_mets_(),
            dot_reacts_(),
            dot_edges_(),
            max_flux_(1.0f),
            out_(nullptr)
        {
            std::vector<bool> avail(scenario->num_metabolites(), true);
            
            source_idx_ =
                scenario->make_dijkstra (graph_, params, avail, supplied_by_);

            // Now redo graph with only those that are supplied as avail
            for (size_t k = 0 ; k < scenario_->num_metabolites(); k++) {
                avail[k] = (supplied_by_[k] != nullptr);
            }

            std::fill (supplied_by_.begin(), supplied_by_.end(), nullptr);
            source_idx_ =
                scenario->make_dijkstra (graph_, params, avail, supplied_by_);
            
            for (auto& react : scenario_->reactions())
                if (react->known_flux() > max_flux_)
                    max_flux_ = react->known_flux();
        }
            

        void set_outfile (ofstream_ptr out) {out_ = out;}
        
        std::string penwidth (const Reaction* react) const {
            auto val = 1.0 + 10.0 * (react->known_flux() / max_flux_);
            return " [penwidth=" + std::to_string(val) + "]";
        }
        void check_biomass (int& count, int& count_ok);
        void full_check_biomass (int& count, int& count_ok);
        void carbon_source_graph (
            const Metabolite* met, int carbon_depth, std::ostream& dot
        );
        void find_loops (double min_flux, int& num_loops, int& num_mets);
        bool find_loop (const Metabolite* met, double min_flux, int& num_mets);
        size_t num_edges() const {return graph_.num_edges();}
        bool find_sp_to (const Metabolite* met);
        const Metabolite* parent (const Metabolite* met) const {
            size_t idx = graph_.parent (met->index());
            if (idx == source_idx_)
                return nullptr;
            if (idx == graph_.num_nodes())
                return nullptr;
            return scenario_->metabolite(idx);
        }
        size_t react_idx_for (const Metabolite* met) const {
            return graph_.parent_edge (met->index());
        }
        const Reaction* react_for (const Metabolite* met) const {
            size_t idx = graph_.parent_edge (met->index());
            return scenario_->reaction(idx);
        }

        bool is_supplied (const Metabolite* met) const
        {
            return supplied_by_[met->index()] != nullptr;
        }

        void add_dot_edge (
            std::string a, std::string b, const Reaction* react = nullptr
        )
        {
            std::string edge_str = "\"" + a + "\" -> \"" + b + "\"";
            if (react != nullptr)
                edge_str += penwidth (react);
            
            dot_edges_.insert (edge_str);
        }
        void write_dot (
            std::ostream& dot, bool explode_react, bool explode_met
        );

        
    private:
        void flood();
        double calc_flux (size_t start_idx);
        
        Scenario* scenario_;
        const Params* params_;
        // Dummy sorce reaction
        ReactionPtr source_;
        // The (first) reaction supplying each metabolite
        std::vector<const Reaction*> supplied_by_;

        lime::Dijkstra<double> graph_;
        size_t source_idx_;

        std::vector<bool> in_loop_;
        std::set<const Metabolite*> dot_mets_;
        std::set<const Reaction*> dot_reacts_;
        std::set<std::string> dot_edges_;
        double max_flux_;

        ofstream_ptr out_;
    };
    
    using PathFinderPtr = std::shared_ptr<PathFinder>;
}
