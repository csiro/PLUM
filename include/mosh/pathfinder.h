/**
 * @file pathfinder.h
 * @brief Pathfinding and metabolic network analysis for gap-filling
 *
 * This file defines the PathFinder class which implements shortest path algorithms
 * for metabolic network analysis, carbon source tracing, and loop detection in the
 * context of flux balance analysis and metabolic gap-filling.
 */

#pragma once

#include <list>
#include <set>
#include <memory>


#include "lime/dijkstra.h"

#include "mosh/scenario.h"


/**
 * @brief Namespace for metabolic optimization and shortest-path heuristics
 */
namespace mosh
{
    /** @brief Shared pointer type for output file streams */
    using ofstream_ptr = std::shared_ptr<std::ofstream>;

    /**
     * @class PathFinder
     * @brief Analyzes metabolic pathways using shortest-path algorithms
     *
     * PathFinder constructs a directed graph representation of metabolic networks
     * and provides methods for finding shortest paths from carbon sources to biomass
     * components, detecting metabolic loops, and generating visualization output in
     * DOT format. It uses Dijkstra's algorithm to compute optimal pathways based on
     * reaction costs and flux constraints.
     */
    class PathFinder {
        
    public:
        /**
         * @brief Constructs a PathFinder and initializes the metabolic network graph
         *
         * Builds a Dijkstra graph representation of the metabolic network by calling
         * make_dijkstra twice: first to identify all supplied metabolites, then to
         * rebuild the graph using only supplied metabolites as available nodes.
         * Also computes the maximum flux across all reactions for visualization scaling.
         *
         * @param scenario Pointer to the metabolic scenario containing reactions and metabolites
         * @param params Pointer to the parameter configuration for gap-filling
         */
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
            

        /**
         * @brief Sets the output file stream for logging and debugging
         * @param out Shared pointer to an output file stream
         */
        void set_outfile (ofstream_ptr out) {out_ = out;}
        
        /**
         * @brief Computes DOT penwidth attribute based on reaction flux
         *
         * Calculates a penwidth value scaled by the reaction's known flux relative
         * to the maximum flux in the network, for visual representation of flux magnitude.
         *
         * @param react Pointer to the reaction whose penwidth is being computed
         * @return String containing the DOT penwidth attribute specification
         */
        std::string penwidth (const Reaction* react) const {
            auto val = 1.0 + 10.0 * (react->known_flux() / max_flux_);
            return " [penwidth=" + std::to_string(val) + "]";
        }
        /**
         * @brief Checks biomass metabolite connectivity in the network
         * @param count Output parameter for total number of biomass metabolites checked
         * @param count_ok Output parameter for number of biomass metabolites with valid paths
         */
        void check_biomass (int& count, int& count_ok);
        /**
         * @brief Performs comprehensive biomass metabolite connectivity analysis
         * @param count Output parameter for total number of biomass metabolites checked
         * @param count_ok Output parameter for number of biomass metabolites with valid paths
         */
        void full_check_biomass (int& count, int& count_ok);
        /**
         * @brief Generates a carbon source tracing graph in DOT format
         *
         * Traces metabolic pathways from a specified metabolite back to carbon sources,
         * following reactions up to a specified depth based on carbon atom flow.
         *
         * @param met Pointer to the starting metabolite for carbon tracing
         * @param carbon_depth Maximum depth for carbon source tracing
         * @param dot Output stream for writing DOT graph representation
         */
        void carbon_source_graph (
            const Metabolite* met, int carbon_depth, std::ostream& dot
        );
        /**
         * @brief Detects metabolic loops in the network with flux above threshold
         * @param min_flux Minimum flux threshold for loop detection
         * @param num_loops Output parameter for number of loops found
         * @param num_mets Output parameter for number of metabolites involved in loops
         */
        void find_loops (double min_flux, int& num_loops, int& num_mets);
        /**
         * @brief Searches for a metabolic loop starting from a specific metabolite
         * @param met Pointer to the starting metabolite for loop detection
         * @param min_flux Minimum flux threshold for considering reactions in the loop
         * @param num_mets Output parameter for number of metabolites in the detected loop
         * @return True if a loop was found, false otherwise
         */
        bool find_loop (const Metabolite* met, double min_flux, int& num_mets);
        /**
         * @brief Returns the number of edges in the metabolic network graph
         * @return Number of edges in the Dijkstra graph
         */
        size_t num_edges() const {return graph_.num_edges();}
        /**
         * @brief Finds shortest path to a specified metabolite
         * @param met Pointer to the target metabolite
         * @return True if a path was found, false otherwise
         */
        bool find_sp_to (const Metabolite* met);
        /**
         * @brief Returns the parent metabolite in the shortest path tree
         *
         * Retrieves the parent metabolite for a given metabolite in the shortest path
         * tree rooted at the carbon source. Returns nullptr if the metabolite has no
         * parent or is directly connected to the source.
         *
         * @param met Pointer to the metabolite whose parent is being queried
         * @return Pointer to the parent metabolite, or nullptr if none exists
         */
        const Metabolite* parent (const Metabolite* met) const {
            size_t idx = graph_.parent (met->index());
            if (idx == source_idx_)
                return nullptr;
            if (idx == graph_.num_nodes())
                return nullptr;
            return scenario_->metabolite(idx);
        }
        /**
         * @brief Returns the reaction index for the edge leading to a metabolite
         * @param met Pointer to the metabolite
         * @return Index of the reaction that produces this metabolite in the shortest path
         */
        size_t react_idx_for (const Metabolite* met) const {
            return graph_.parent_edge (met->index());
        }
        /**
         * @brief Returns the reaction that produces a metabolite in the shortest path
         * @param met Pointer to the metabolite
         * @return Pointer to the reaction that produces this metabolite in the shortest path
         */
        const Reaction* react_for (const Metabolite* met) const {
            size_t idx = graph_.parent_edge (met->index());
            return scenario_->reaction(idx);
        }

        /**
         * @brief Checks if a metabolite is supplied by any reaction
         * @param met Pointer to the metabolite to check
         * @return True if the metabolite is supplied by at least one reaction
         */
        bool is_supplied (const Metabolite* met) const
        {
            return supplied_by_[met->index()] != nullptr;
        }

        /**
         * @brief Adds an edge to the DOT graph representation
         *
         * Constructs a DOT edge string between two nodes and optionally applies
         * penwidth styling based on reaction flux.
         *
         * @param a Source node name
         * @param b Target node name
         * @param react Optional pointer to reaction for flux-based styling (default: nullptr)
         */
        void add_dot_edge (
            std::string a, std::string b, const Reaction* react = nullptr
        )
        {
            std::string edge_str = "\"" + a + "\" -> \"" + b + "\"";
            if (react != nullptr)
                edge_str += penwidth (react);

            dot_edges_.insert (edge_str);
        }
        /**
         * @brief Writes the complete DOT graph representation to an output stream
         * @param dot Output stream for DOT format graph
         * @param explode_react If true, expand reactions into explicit nodes
         * @param explode_met If true, expand metabolites into explicit nodes
         */
        void write_dot (
            std::ostream& dot, bool explode_react, bool explode_met
        );

        
    private:
        /**
         * @brief Performs flood-fill traversal of the metabolic network
         */
        void flood();
        /**
         * @brief Calculates flux starting from a specified node index
         * @param start_idx Index of the starting node in the graph
         * @return Calculated flux value
         */
        double calc_flux (size_t start_idx);
        
        /** @brief Pointer to the metabolic scenario containing reactions and metabolites */
        Scenario* scenario_;
        /** @brief Pointer to the parameter configuration for gap-filling */
        const Params* params_;
        /** @brief Dummy source reaction representing the carbon source entry point */
        ReactionPtr source_;
        /** @brief Vector storing the first reaction that supplies each metabolite */
        std::vector<const Reaction*> supplied_by_;

        /** @brief Dijkstra graph structure for shortest path computation */
        lime::Dijkstra<double> graph_;
        /** @brief Index of the source node in the Dijkstra graph */
        size_t source_idx_;

        /** @brief Vector tracking which metabolites are involved in metabolic loops */
        std::vector<bool> in_loop_;
        /** @brief Set of metabolites to include in DOT graph output */
        std::set<const Metabolite*> dot_mets_;
        /** @brief Set of reactions to include in DOT graph output */
        std::set<const Reaction*> dot_reacts_;
        /** @brief Set of edge strings for DOT graph representation */
        std::set<std::string> dot_edges_;
        /** @brief Maximum flux value across all reactions for visualization scaling */
        double max_flux_;

        /** @brief Shared pointer to output file stream for logging and debugging */
        ofstream_ptr out_;
    };
    
    /** @brief Shared pointer type for PathFinder objects */
    using PathFinderPtr = std::shared_ptr<PathFinder>;
}
