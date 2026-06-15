/**
 * @file pathsol.h
 * @brief Path solution finder for metabolic network analysis
 *
 * This file defines the PathSol class which computes shortest paths between
 * metabolites in a metabolic network using Dijkstra's algorithm. It is used
 * in gap-filling and metabolic reconstruction to identify missing reactions
 * needed to connect source and target metabolites.
 */
#pragma once

#include <vector>
#include <memory>
#include <set>

#include "mosh/scenario.h"
#include "mosh/reaction.h"
#include "mosh/path.h"

#include "lime/dijkstra.h"

/**
 * @namespace mosh
 * @brief Metabolic optimization and stoichiometric handling namespace
 */
namespace mosh
{
    /**
     * @class PathSol
     * @brief Finds shortest paths between metabolites in a metabolic network
     *
     * PathSol implements path-finding capabilities for metabolic networks using
     * Dijkstra's algorithm. It operates on a Scenario (metabolic network state)
     * to compute paths between source and target metabolites, which is essential
     * for gap-filling analysis to identify missing reactions or pathways.
     *
     * The class maintains a reference to a Scenario and reuses a Dijkstra solver
     * instance for efficient repeated path computations.
     */
    class PathSol 
    {
    public:
        /**
         * @brief Constructs a PathSol instance for a given metabolic scenario
         * @param scenario Pointer to the metabolic network scenario to analyze
         */
        PathSol(const Scenario* scenario) :
            scenario_(scenario)
        {
        }
        /**
         * @brief Copy constructor from another PathSol instance
         * @param other Pointer to the PathSol instance to copy from
         */
        PathSol(PathSol* other) :
            scenario_(other->scenario_)
        {
        }

        /**
         * @brief Gets the associated metabolic scenario
         * @return Const pointer to the Scenario being analyzed
         */
        const Scenario* scenario() const {return scenario_;}

        /**
         * @brief Computes the shortest path between two metabolites
         * @param from Source metabolite to start the path from
         * @param to Target metabolite to reach
         * @param path Output parameter that will contain the computed path
         * @return true if a valid path was found, false otherwise
         *
         * Uses Dijkstra's algorithm to find the shortest path through the metabolic
         * network from the source to target metabolite. The resulting path contains
         * the sequence of reactions connecting the metabolites.
         */
        bool calc_path (Metabolite* from, Metabolite* to, Path& path);

    private:
        lime::Dijkstra<int> dijkstra_; /**< Dijkstra's algorithm solver instance for path computation */
        
        const Scenario* scenario_; /**< Pointer to the metabolic scenario being analyzed */
        Path path_; /**< Internal path storage for computation results */
    };

    /**
     * @typedef PathSolPtr
     * @brief Shared pointer type for PathSol instances
     */
    using PathSolPtr = std::shared_ptr<PathSol>;
}
