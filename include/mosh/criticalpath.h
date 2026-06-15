/**
 * @file criticalpath.h
 * @brief Critical path analysis for metabolic networks
 *
 * This file defines the CriticalPath class which identifies critical metabolic
 * pathways from a carbon source to target metabolites in gap-filled metabolic
 * networks. Uses Dijkstra's algorithm to find shortest paths through the
 * metabolic network based on flux balance analysis solutions.
 */

#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <memory>

#include "lime/dijkstra.h"
#include "lime/numutil.h"
#include "lime/dijkstra.h"

#include "mosh/scenario.h"
#include "mosh/solution.h"
#include "mosh/reaction.h"

/**
 * @brief Metabolic gap-filling namespace
 *
 * Contains classes and utilities for metabolic network gap-filling,
 * flux balance analysis, and pathway analysis.
 */
namespace mosh
{
    /**
     * @class CriticalPath
     * @brief Identifies critical metabolic pathways in gap-filled networks
     *
     * This class analyzes flux balance analysis solutions to identify the most
     * critical metabolic pathways from a carbon source to target metabolites.
     * It constructs a directed graph from active reactions in the solution and
     * uses Dijkstra's algorithm to find shortest paths through the metabolic network.
     *
     * Critical paths help identify essential reactions and metabolites that are
     * necessary for biomass production or other metabolic objectives.
     */
    class CriticalPath 
    {
    public:
        /**
         * @brief Construct a CriticalPath analyzer
         *
         * Initializes the critical path analysis for a given FBA solution and
         * carbon source metabolite. Constructs the metabolic network graph from
         * active reactions in the solution.
         *
         * @param sol Shared pointer to the FBA solution to analyze
         * @param carbon_source Pointer to the carbon source metabolite serving as the
         *                      starting point for path analysis
         */
        CriticalPath (
            SolutionPtr sol,
            const Metabolite* carbon_source
        ) :
            sol_(sol),
            scenario_(sol->scenario()),
            carbon_source_(carbon_source),
            graph_(scenario_->num_metabolites())
        {
            init_graph();
        }

        /**
         * @brief Find and output critical paths to target metabolites
         *
         * Computes shortest paths from the carbon source to all reachable metabolites
         * in the network using Dijkstra's algorithm. Outputs the identified critical
         * paths including path length and intermediate metabolites.
         *
         * @param out Output stream to write the discovered paths and analysis results
         */
        void find_paths(std::ostream& out);

    private:
        /**
         * @brief Initialize the metabolic network graph
         *
         * Constructs a directed graph representation of the metabolic network from
         * the active reactions in the FBA solution. Each metabolite becomes a node,
         * and each reaction with non-zero flux creates directed edges between
         * substrate and product metabolites.
         */
        void init_graph();
        
        SolutionPtr sol_; /**< Shared pointer to the FBA solution being analyzed */
        const Scenario* scenario_; /**< Pointer to the metabolic scenario containing network structure */
        const Metabolite* carbon_source_; /**< Pointer to the carbon source metabolite (starting point) */
        lime::Dijkstra<double> graph_; /**< Dijkstra graph for shortest path computation through metabolic network */
    };
}
